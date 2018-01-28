## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('corr')

## ----calcr---------------------------------------------------------------
n <- 13
r <- 0.7283
z.transform <- 0.5 * log((1 + r) / (1 - r))
clz <- z.transform + c(-1, 1) * qnorm(0.975) / sqrt(n - 3)
clr <- (exp(2 * clz) - 1) / (exp(2 * clz) + 1)
round(c(z.transform, clz, clr), 4)

## ----corrplota,w=6.5,h=5,cap='Samples of size $n=50$ for X and Y are drawn from bivariate normal populations with true correlations ranging from 0.0 to 0.9. Pearson and Spearman sample correlations are shown for samples of size 50.  Besides the population correlation coefficient, each panel is labeled with the estimated Pearson $r$, its $t$ statistic, the estimated Spearman $\\rho$, and its $t$ statistic',scap='Example correlation coefficients'----
# Generate 50 data points with Popluation correlations of 0, .2, .4, .6,
# .8, and .9 and plot results
require(ggplot2)
n <- 50
set.seed(123)
x <- rnorm(n, 5, 1)
d <- expand.grid(x=x, R=c(0, .2, .4, .6, .8, .9))
d <- transform(d, y = x + rnorm(nrow(d), 0,
                      ifelse(R == 0, 5, sqrt(R ^ -2 - 1))))
sfun <- function(i) {
  x <- d$x[i]; y <- d$y[i]; R <- d$R[i][1]
  r <- cor(x, y)
  tr <- r * sqrt(n - 2) / sqrt(1 - r^2)
  rho <- cor(rank(x), rank(y))
  trho <- rho * sqrt(n - 2) / sqrt(1 - rho^2)
  label <- paste('True r:', R[1], '  r:', round(r,2), '  t:', round(tr,2),
        '  rho:', round(rho,2), '  t:', round(trho,2), sep='')
  names(label) <- R
  label
  }
stats <- tapply(1 : nrow(d), d$R, sfun)
d$stats <- factor(stats[as.character(d$R)], unique(stats))
   
ggplot(d, aes(x=x, y=y)) + geom_point() + facet_wrap(~ stats) +
  theme(strip.text.x = element_text(size=7))   # Fig. (*\ref{fig:corr-corrplota}*)

## ----corrplotb,mfrow=c(3,2),mar=c(3,1,.5,1),w=6,h=6,cap='Different observed datasets that have the same correlation.  All six plots have a sample Pearson correlation of $0.7$.',scap='Multiple datasets having same Pearson $r$'----
# Different scenarios that can lead to a correlation of 0.7

set.seed(123)   # Fig. (*\ref{fig:corr-corrplotb}*)
rho <- 0.7; n <- 50
var.eps <- rho^-2 - 1
x <- rnorm(n, 5, 1)
y <- x + rnorm(n, 0, sqrt(var.eps))
cor(x,y)
plot(x,y,xlab='',ylab='')

x <- c(1:20,30)
y <- c(1:20,6.2)
cor(x,y)
plot(x,y,xlab='',ylab='')

set.seed(123)
x <- rnorm(40)
y <- rnorm(40)
x[21] <- y[21] <- 8.5
cor(x,y)
plot(x,y,xlab='',ylab='')

x <- rep(0:19,2)
y <- c(rep(.62,20),rep(2,20)) * x
cor(x,y)
plot(x,y,xlab='',ylab='')

x <- -7:12
y <- x^2
cor(x,y)
plot(x,y,xlab='',ylab='')

set.seed(123)
tmp <- 1:20 / 2
x <- c(rnorm(20, tmp, 1), tmp + rnorm(20,14.5,1))
y <- c(rnorm(20, -tmp, 1), -tmp + rnorm(20,14.5,1))
cor(x,y)
plot(x,y,xlab='',ylab='')

## ----baplot,w=5.5,h=4.5,cap='Bland-Altman plot for the oroesophageal and conventional pH measurements, using hexagonal binning because of the large sample size.  The difference in pH mesaurements (oro.\\ -conventional) is presented on the $y$-axis and the average of the two devices on the $x$-axis.  We see poor agreement around pH values of 4-5',scap='Bland-Altman plot for 2 pH measurements'----
require(Hmisc)
getHdata(esopH)
esopH$diff <- with(esopH, orophar - conv)
ggplot(esopH, aes(x=(conv + orophar)/2, y=diff)) +   # Fig. (*\ref{fig:corr-baplot}*)
  stat_binhex(aes(alpha=..count.., color=Hmisc::cut2(..count.., g=20)),
              bins=80) +
  stat_smooth() +
  geom_hline(yintercept = mean(esopH$diff, na.rm=TRUE) +
   c(-1.96, 0, 1.96) * sd(esopH$diff, na.rm=TRUE),
   linetype=c(2,1,2), color='brown') +
  xlab('Average of Conventional and Oropharyngeal pH') +
  ylab('Oropharyngeal Minus Conventional pH') +
  guides(alpha=FALSE, fill=FALSE, color=guide_legend(title='Frequency'))

## ----phtimediff,w=5,h=4,cap='Difference in pH measurements (oro.\ - conventional) by time of day along with a loess smoother and pointwise 0.95 confidence bands.  Is the difference modified by a subject being in a supine position rather than being upright?',scap='Difference in pH by time of day'----
getHdata(esopH2)
ggplot(esopH2, aes(x=time, y=diffpH)) +    # Fig. (*\ref{fig:corr-phtimediff}*)
       geom_point(pch='.') + stat_smooth() +
       geom_hline(yintercept = 0, col='gray60') +
       scale_x_continuous(breaks=seq(16, 38, by=4),
                          labels=c("4 PM", "8 PM", "12 AM",
                            "4 AM", "8AM", "12 PM"),
                          limits=c(14, 14+24)) +
       ylab('Average of Oropharyngeal Minus Conventional pH') +
       xlab('Time of Day')

## ----moe,w=7,h=5.5,cap='Margin for error (length of longer side of asymmetric 0.95 confidence interval) for $r$ in estimating $\\rho$, when $\\rho=0, 0.25, 0.5, 0.75$.  Calculations are based on Fisher $z$ transformation of $r$.',scap='Margin of error for estimating correlation coefficient'----
require(Hmisc)
plotCorrPrecision(rho=c(0, .25, .5, .75),
                  n=seq(10, 1000, length=100),
                  ylim=c(0, .4), col=1:4, opts=list(keys='lines'))
abline(h=seq(0, .4, by=0.025),
       v=seq(25, 975, by=25), col=gray(.9))

