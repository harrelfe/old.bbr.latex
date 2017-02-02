## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('htest')

## ----tpdfs,w=5,h=3.75,cap='Comparison of probability densities for $t_2$, $t_5$, $t_{50}$, and normal distributions',scap='$t$ distribution for varying d.f.'----
x <- seq(-3, 3, length=200)   # Fig. (*\ref{fig:htest-tpdfs}*)
w <- list(Normal     = list(x=x, y=dnorm(x)),
          't (50 df)'= list(x=x, y=dt(x, 50)),
          't (5 df)' = list(x=x, y=dt(x, 5)),
          't (2 df)' = list(x=x, y=dt(x,2)))
require(Hmisc)
labcurve(w, pl=TRUE, keys='lines', col=c(1,2,4,5), lty=c(2,1,1,1),
         xlab=expression(x), ylab='')

## ----tone----------------------------------------------------------------
xbar  <- 181.52
s     <- 40
n     <- 100
mu0   <- 190
tstat <- (xbar - mu0) / (s / sqrt(n))
pval  <- 2 * (1 - pt(abs(tstat), n - 1))
c(tstat=tstat, pval=pval)

## ----powerone------------------------------------------------------------
sigma <- 1
mu    <- 2.5
mu0   <- 3
n     <- 10.51 * (1 / (mu - mu0)) ^ 2
# General formula for approximate power of 1-sample t-test
# Approximate because it uses the normal distribution throughout,
# not the t distribution
alpha <- 0.05
power <- 0.9
delta <- mu - mu0
za    <- qnorm(1 - alpha / 2)
zb    <- qnorm(power)
n     <- ((za + zb) * sigma / delta) ^ 2
c(alpha=alpha, power=power, delta=delta, za=za, zb=zb, n=n)

## ----pwrpacktone---------------------------------------------------------
# Make sure pwr package is installed
require(pwr)
pwr.t.test(d = delta / sigma, power = 0.9, sig.level = 0.05, type='one.sample')

## ----precmean------------------------------------------------------------
sigma <- 10
delta <- 1
3.84 * (sigma / delta) ^ 2

## ----zprop---------------------------------------------------------------
p  <- 0.8
p0 <- 0.5
n  <- 10
z  <- (p - p0) / sqrt(p0 * (1 - p0) / n)
c(z=z, Pvalue=2 * pnorm(-abs(z)))

## ----exactprob-----------------------------------------------------------
# Pr(X >= 8) = 1 - Pr(X < 8) = 1 - Pr(X <= 7)
pbinom(2, 10, 0.5) + 1 - pbinom(7, 10, 0.5)
# Also compute as the probability of getting 0, 1, 2, 8, 9, 10 heads
sum(dbinom(c(0, 1, 2, 8, 9, 10), 10, 0.5))

## ----nmoep---------------------------------------------------------------
delta <- 0.1
n <- 0.25 * (qnorm(0.975) / delta) ^ 2
n

## ----ptwins--------------------------------------------------------------
xbar  <- -5
se    <- 2
n     <- 41
mu0   <- 0
tstat <- (xbar - mu0) /se
pval  <- 2 * (1 - pt(abs(tstat), n - 1))
c(tstat=tstat, Pvalue=pval)

## ----ptwot---------------------------------------------------------------
n1    <- 8;         n2 <- 21
xbar1 <- 132.86; xbar2 <- 127.44
s1    <- 15.34;     s2 <- 18.23
s     <- sqrt(((n1 - 1) * s1 ^ 2 + (n2 - 1) * s2 ^ 2) / (n1 + n2 - 2))
se    <- s * sqrt(1 / n1 + 1 / n2)
tstat <- (xbar1 - xbar2) / se
pval  <- 2 * (pt(- abs(tstat), n1 + n2 - 2))
c(s=s, se=se, tstat=tstat, Pvalue=pval)

## ----tuneqpow------------------------------------------------------------
delta <- 5
require(pwr)
pwr.t2n.test(n1=100, n2=100, d=delta / s, sig.level = 0.05)

## ----ssizet--------------------------------------------------------------
pwr.t.test(d = delta / s, sig.level = 0.05, power = 0.8)

## ----powtuneqn-----------------------------------------------------------
pwr.t2n.test(n1 = 129, n2 = 259, d = delta / s, sig.level = 0.05)

## ----smmoe,cap='Multiplicative margin of error in estimating $\\sigma$ as a function of sample size, with 0.95 confidence',scap='Margin of error in estimating $\\sigma$'----
n    <- 10:300
low  <- sqrt((n - 1) / qchisq(.975, n - 1))
hi   <- sqrt((n - 1) / qchisq(.025, n - 1))
m    <- pmax(1 / low, hi)
ggplot(data.frame(n, m), aes(x=n, y=m)) + geom_line() +
  ylab('MMOE for s')
nmin <- min(n[m <= 1.2])

## ----sleepa,w=5,h=3,cap='Data for two-sample RCT.  Measurements are differences from a control period while subjects were on placebo.  Control period data were not used in the analysis except for normalization.  Diamonds depict means.',scap='Two-sample parallel group RCT'----
require(ggplot2)   # Fig (*\ref{fig:htest-sleepa}*)
sleep <- data.frame(drug=c(rep('Drug 1', 10), rep('Drug 2', 10)),
                    extra=c(.7, -1.6, -.2, -1.2, -.1, 3.4, 3.7, .8, 0, 2,
                      1.9, .8, 1.1, .1, -.1, 4.4, 5.5, 1.6, 4.6, 3.4))
ggplot(sleep, aes(x=drug, y=extra, group=drug)) +
  geom_boxplot(col='lightyellow1', alpha=.3, width=.3) + 
  geom_dotplot(binaxis='y', stackdir='center', position='dodge') +
  stat_summary(fun.y=mean, geom="point", col='red', shape=5, size=3) +
  xlab('') + ylab('Extra Hours of Sleep') + coord_flip() 

## ------------------------------------------------------------------------
t.test(extra ~ drug, data=sleep)

## ----tplot,w=5,h=3,cap='Raw data and box plots for paired data and their paired differences, with lines connecting points from the same subject.  Diamonds depict means.',scap='Data and box plots for paired data'----
drug1 <- c(.7, -1.6, -.2, -1.2, -.1, 3.4, 3.7, .8, 0, 2)
drug2 <- c(1.9, .8, 1.1, .1, -.1, 4.4, 5.5, 1.6, 4.6, 3.4)
d <- data.frame(Drug=c(rep('Drug 1', 10), rep('Drug 2', 10),
                  rep('Difference', 10)),
                extra=c(drug1, drug2, drug2 - drug1))

ggplot(d, aes(x=Drug, y=extra)) +   # Fig. (*\ref{fig:htest-tplot}*)
  geom_boxplot(col='lightyellow1', alpha=.3, width=.5) + 
  geom_dotplot(binaxis='y', stackdir='center', position='dodge') +
  stat_summary(fun.y=mean, geom="point", col='red', shape=18, size=5) +
  geom_segment(aes(x='Drug 1', xend='Drug 2', y=drug1, yend=drug2),
               col=gray(.8)) +
  geom_segment(aes(x='Drug 1', xend='Difference', y=drug1, yend=drug2 - drug1),
               col=gray(.8)) +
  xlab('') + ylab('Extra Hours of Sleep') + coord_flip() 

## ----pairedt-------------------------------------------------------------
with(d, t.test(drug1, drug2, paired=TRUE))

## ----xover---------------------------------------------------------------
drug1 <- c(87, 64, 78, 68, 79, 114, 117, 88, 80, 100)/10
drug2 <- c(99, 88, 91, 81, 79, 124, 135, 96, 126, 114)/10
t.test(drug1, drug2, paired=TRUE)

## ----xovercarry----------------------------------------------------------
# Unpaired t-test
t.test((drug2 - drug1)[1:5], (drug2 - drug1)[6:10], var.equal=TRUE)

