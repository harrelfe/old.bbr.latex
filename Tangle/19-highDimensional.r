## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('hdata', cache=TRUE)

## ----simor-sim-----------------------------------------------------------
# For a vector of n binary outcomes y, simulates p binary features
# x that have a p-vector of fixed prevalences | y=0 of prev and are connected
# to y by a p-vector of true population odds ratios ors.
# Estimates the p odds ratios against the simulated outcomes and
# returns a data frame summarizing the information
#
# Note: the odds ratio for x predicting y is the same as the odds ratio
# for y predicting x.  y is simulated first so that all features will
# be analyzed against the same outcomes

sim <- function(y, prev, or) {
  n <- length(y)
  p <- length(prev)
  if(p != length(or)) stop('prev and or must have the same length')

  # prev = Pr(x=1 | y=0); let the odds for this be oprev = prev / (1-prev)
  # or = odds(x=1 | y=1) / oprev
  # Pr(x=1 | y=1) = oprev / ((1 / or) + oprev)

  oprev <- prev / (1 - prev)
  p1 <- oprev / ((1 / or) + oprev)
  n0 <- sum(y == 0)
  n1 <- sum(y == 1)
  # For n0 observations sample x so that Pr(x=0 | y=0) = prev
  nxy0 <- rbinom(p, size=n0, prob=prev)
  nxy1 <- rbinom(p, size=n1, prob=p1)

  # Compute p sample odds ratios
  sor <- (n0 - nxy0) * nxy1 / (nxy0 * (n1 - nxy1))
  g <- function(x) ifelse(x >= 1, x, 1 / x)
  r1 <- rank(sor)[which.max(or) ] / p
  r2 <- rank(or) [which.max(sor)] / p
  data.frame(prev, or, nx=nxy0 / n0, obsprev0=(nxy0 + nxy1) / n,
             obsprev=nxy1 / (nxy0 + nxy1), obsor=sor, n=n,
             N       =paste('n', n, sep=':'),
             Features=paste('Features', p, sep=':'),
             mmoe    =quantile(g(sor / or), 0.90, na.rm=TRUE),
             obsranktrue=r1, truerankobs=r2,
             rho=cor(sor, or, method='spearman', use='pair'))
}

## ----simor-dosim---------------------------------------------------------
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

## ----simor-plot,w=7,h=7--------------------------------------------------
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

## ----simor-plotb,w=7,h=7-------------------------------------------------
pl(0.3)

## ----simor-mmoe,w=5.5,h=5.5----------------------------------------------
ggplot(U, aes(x=n, y=mmoe)) + geom_point() + facet_grid(Features ~ Yprev) +
  geom_hline(aes(yintercept=1.5, col='red')) +
  ylim(1, 10) +
  ylab('0.9 Quantile of Multiplicative Margin of Error in OR Across Features')

## ----simor-rho,w=5.5,h=5.5-----------------------------------------------
ggplot(U, aes(x=n, y=rho)) + geom_point() +
  facet_grid(Features ~ Yprev) +
  ylab(expression(paste('Spearman ', rho, ' Rank Correlation Between ',
                        OR, ' and ', hat(OR), ' Across Features')))

## ----simbf---------------------------------------------------------------
# Function to simulate the raw data
# prev is the vector of prevalences of x when y=0 as before
# yprev is the overall prevalence of y
# n is the sample size
# or is the vector of true odds ratios
sim <- function(n, yprev, prev, or) {
  y <- rbinom(n, 1, yprev)
  p <- length(prev)
  if(p != length(or)) stop('prev and or must have the same length')

  # prev = Pr(x=1 | y=0); let the odds for this be oprev = prev / (1-prev)
  # or = odds(x=1 | y=1) / oprev
  # Pr(x=1 | y=1) = oprev / ((1 / or) + oprev)

  oprev <- prev / (1 - prev)
  p1 <- oprev / ((1 / or) + oprev)
  x <- matrix(NA, nrow=n, ncol=p)
  for(j in 1 : p) 
  	x[, j] <- ifelse(y == 1, rbinom(n, 1, prob = p1[j]  ), 
                             rbinom(n, 1, prob = prev[j]))
  list(x=x, y=y)
}
  
# Function to compute the sample odds ratios given x matrix and y vector
ors <- function(x, y) {
  p <- ncol(x)
  or <- numeric(p)
  for(j in 1 : p) {
    f <- table(x[, j], y)
    or[j] <- f[2, 2] * f[1, 1] / (f[1, 2] * f[2, 1])
  }
  or
}

## ----simbd---------------------------------------------------------------
# Generate sample of size 600 with 300 features
# Log odds ratios have a normal distribution with mean 0 SD 0.3
# x have a random prevalence uniform [0.05, 0.5]
# y has prevalence 0.3

set.seed(188)	
n <- 600; p <- 300
prev <- runif(p, .05, .5)
or   <- exp(rnorm(p, 0, .3))
z    <- sim(n, 0.3, prev, or)

# Compute estimated ORs
x <- z$x;   y <- z$y
sor <- ors(x, y)
# Show how estimates related to true ORs
ggplot(data.frame(or, sor), aes(x=or, y=sor)) + geom_point() +
  xlab('True OR') + ylab('Estimated OR')

# Print the largest estimated OR and its column number,
# and corresponding true OR, and similarly for the smallest.
largest   <- max(sor)
imax      <- which.max(sor)
true.imax <- or[imax]
mmoe.imax <- largest / true.imax
smallest  <- min(sor)
imin      <- which.min(sor)
true.imin <- or[imin]
mmoe.imin <- smallest / true.imin
cat('\nLargest observed OR\n')
cat('OR:', round(largest, 2), '  Feature #', imax, '  True OR:',
    round(true.imax, 2), '  MMOE:', round(mmoe.imax, 2), '\n')
cat('Rank of winning feature among true ORs:', sum(or <= or[imax]), '\n\n')
cat('Smallest observed OR\n')
cat('OR:', round(smallest, 2), '  Feature #', imin, '  True OR:',
    round(true.imin, 2), '  MMOE:', round(mmoe.imin, 2), '\n')

## ----bootor--------------------------------------------------------------
set.seed(11)
B <- 1000
ranksS <- ranksL <- mmoeS <- mmoeL <- numeric(B)

for(k in 1 : B) {
  # Draw a sample of size n with replacement
  i <- sample(1 : n, n, replace=TRUE)
  # Compute sample ORs on the new sample
  bor      <- ors(x[i, ], y[i])
  blargest <- max(bor)
  bmax     <- which.max(bor)
  ranksL[k] <- sum(bor <= largest)
  mmoeL[k]  <- blargest / sor[bmax]
  bsmallest <- min(bor)
  bmin      <- which.min(bor)
  ranksS[k] <- sum(bor <= smallest)
  mmoeS[k]  <- bsmallest / sor[bmin]
}

## ----summarizeb----------------------------------------------------------
pr <- function(which, ranks, mmoe, mmoe.true, estor, or.true) {
  gm <- exp(mean(log(mmoe)))
  cat(which, 'OR\n')
  cat('CL for rank:', quantile(ranks, c(0.025, 0.975)),
      '  Median MMOE:', round(median(mmoe), 2),
      '  Geometric mean MMOE:', round(gm, 2),
      '\nTrue MMOE:', round(mmoe.true, 2), '\n')
  bmmoe <- if(which == 'Largest') gm else median(mmoe)
  cat('Bootstrap bias-corrected', tolower(which), 'OR:',
      round(estor / bmmoe, 2),
      '  Original OR:', round(estor, 2),
      '  True OR:', round(or.true, 2),
      '\n\n')
}
pr('Largest',  ranksL, mmoeL, mmoe.imax, largest,  true.imax)
pr('Smallest', ranksS, mmoeS, mmoe.imin, smallest, true.imin)

