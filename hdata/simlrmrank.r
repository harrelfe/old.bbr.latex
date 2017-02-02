## Simulate data using a binary logistic model with known p-vector of
## betas.  Among other things show the variation in which coefficient is
## largest in absolute value
## X is the covariate matrix, bint is the intercept
siml <- function(beta, bint, X) {
  n <- nrow(X)
  L <- bint + X %*% beta
  y <- ifelse(runif(n) <= plogis(L), 1, 0)
  f <- lrm.fit(X, y, tol=1e-10)
  if(f$fail) return(c(n=n, prevy=mean(y), rankbest=NA))
  b <- f$coefficients[-1]
  bhi <- b[which.max(beta)]
  ## Get observed rank of estimated beta that has the true largest beta
  ## (1=highest, 2=next highest, etc.)
  r <- sum(b >= bhi)
  c(n=n, prevy=mean(y), rankbest=r)
}

## Simulate binary predictor matrix having true prevalences ranging
## uniformly from low to high.  Xs are centered by subtracting the mean.
simX <- function(n, p, low=0.05, hi=0.5) {
  prev <- runif(p, low, hi)
  X <- matrix(NA, nrow=n, ncol=p)
  for(j in 1 : p) {
    x <- ifelse(runif(n) <= prev[j], 1, 0)
    X[, j] <- x - mean(x)
  }
  X
}

require(rms)

## Do B simulations for a given n and p where 3/4 of the betas
## are truly zero, all but one of the remaining are normal with SD 0.15,
## and the last true beta is log(2).  Intercept is -1 (Xs are centered).
## Result is B ranks of the truly largest beta among the estimated betas
simr <- function(n, p, B) {
  nzeros <- floor(3 * p / 4)
  ngauss <- p - nzeros - 1
  beta <- c(rep(0, nzeros), rnorm(ngauss, sd=0.15), log(2))
  X <- simX(n, p)
  r <- numeric(B)
  for(i in 1 : B) {
    s <- siml(beta, -1, X)
    rb <- s['rankbest']
    yprev <- s['prevy']
    if(i < 6) cat('Simulated prevalance of Y:', round(yprev, 2),
                  '  Rank best:', rb, '\n')
    else if(i %% 10 == 0) cat(i, '')
    r[i] <- rb
  }
  cat('\n')
  data.frame(n, p, yprev, rank=r)
}

set.seed(1)
U <- NULL
for(n in c(50, 100, 250, 500, 1000, 2000)) {
  for(p in c(10, 50, 100, 500)) {
    if(p >= n / 3) next
    cat('n:', n, '  p:', p, '\n')
    u <- simr(n=n, p=p, B=1000)
    cat('Number of simulations where model could not be fit:',
        sum(is.na(u$rank)), '\n')
    U <- rbind(U, u)
  }
  simrank <- U
  Save(simrank)   # save every time n changes so can monitor results
}

graph <- FALSE
if(graph) {
system('scp -p $vu:doc/teaching/hes/CI2/hdata/simrank.rda .')
Load(simrank)
dim(simrank)

s <- transform(simrank,
               p = factor(p), n = factor(n))
levels(s$p) <- paste('Features', levels(s$p), sep=':')
levels(s$n) <- paste('n', levels(s$n), sep=':')
pdf('simlrmrank.pdf', width=6, height=4.5)
ggplot(s, aes(x=rank)) + geom_histogram() +
  facet_grid(p ~ n, scales='free_x') +
  xlab(expression(paste('Rank of True Strongest Feature Among ',
       hat(beta), 's'))) +
  ylab('Frequency') +
  theme(axis.text.x = element_text(size=rel(0.8), angle=-45, hjust=0, vjust=1))
i <- round(mean(s$yprev), 2)
prfootnote(sprintf('1000 simulations, mean outcome incidence %s', i))
dev.off()
}
