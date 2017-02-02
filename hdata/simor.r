## For a vector of n binary outcomes y, simulates p binary features
## x that have a p-vector of fixed prevalences | y=0 of prev and are connected
## to y by a p-vector of true population odds ratios ors.
## Estimates the p odds ratios against the simulated outcomes and
## returns a data frame summarizing the information
##
## Note: the odds ratio for x predicting y is the same as the odds ratio
## for y predicting x.  y is simulated first so that all features will
## be analyzed against the same outcomes

sim <- function(y, prev, or) {
  n <- length(y)
  p <- length(prev)
  if(p != length(or)) stop('prev and or must have the same length')

  ## prev = Pr(x=1 | y=0); let the odds for this be oprev = prev / (1-prev)
  ## or = odds(x=1 | y=1) / oprev
  ## Pr(x=1 | y=1) = oprev / ((1 / or) + oprev)

  oprev <- prev / (1 - prev)
  p1 <- oprev / ((1 / or) + oprev)
  n0 <- sum(y == 0)
  n1 <- sum(y == 1)
  ## For n0 observations sample x so that Pr(x=0 | y=0) = prev
  nxy0 <- rbinom(p, size=n0, prob=prev)
  nxy1 <- rbinom(p, size=n1, prob=p1)

  ## Compute p sample odds ratios
  sor <- (n0 - nxy0) * nxy1 / (nxy0 * (n1 - nxy1))
  g <- function(x) ifelse(x >= 1, x, 1 / x)
  r1 <- rank(sor)[which.max(or) ] / p
  r2 <- rank(or) [which.max(sor)] / p
  data.frame(prev, or, nx=nxy0 / n0, obsprev0=(nxy0 + nxy1) / n,
             obsprev=nxy1 / (nxy0 + nxy1), obsor=sor, n=n,
             N       =paste('n', n, sep=':'),
             Features=paste('Features', p, sep=':'),
             mmoe    =quantile(g(sor / or), 0.99, na.rm=TRUE),
             obsranktrue=r1, truerankobs=r2,
             rho=cor(sor, or, method='spearman', use='pair'))
}
U <- NULL
set.seed(1)
for(n in c(50, 100, 250, 500, 1000, 2000)) {
  for(p in c(10, 50, 500, 1000, 2000)) {
    for(yprev in c(.1, .3)) {
      y <- rbinom(n, 1, yprev)
      prev <- runif(p, .05, .5)
      or <- exp(rnorm(p, 0, .25))
      u <- cbind(sim(y, prev, or),
                 Yprev=paste('Prevalence of Outcome', yprev, sep=':'))
      U <- rbind(U, u)
    }
  }
}

require(ggplot2)
pl <- function(yprev) {
  br <- c(.01, .1, .5, 1, 2.5, 5, 25, 100)
  ggplot(subset(U, Yprev=yprev),
         aes(x=or, y=obsor)) + geom_point() + facet_grid(Features ~ N) +
         ggtitle(paste('Prevalence of Outcome', yprev, sep=':')) +
         xlab('True ORs') + ylab('Estimated ORs') +
         scale_x_log10(breaks=br) + scale_y_log10(breaks=br) +
         theme(axis.text.x = element_text(size = rel(0.8), angle=-45,
                                    hjust=0, vjust=1)) +
         geom_abline(col='red')
}
pl(0.1)
pdf('simor-ors.pdf', width=8, height=8)
pl(0.3)
dev.off()

ggplot(U, aes(x=n, y=mmoe)) + geom_point() + facet_grid(Features ~ Yprev) +
  ylim(1, 10) +
  ylab('0.9 Quantile of Multiplicative Margin of Error in OR Across Features')

ggplot(U, aes(x=n, y=rho)) + geom_point() +
  facet_grid(Features ~ Yprev) +
  ylab(expression(paste('Spearman ', rho, ' Rank Correlation Between ',
                        OR, ' and ', hat(OR), ' Across Features')))

