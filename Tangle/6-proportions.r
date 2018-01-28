## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('prop')

## ----ptwo----------------------------------------------------------------
n1 <- 3220;     n2 <- 10245
p1 <- 683 / n1; p2 <- 1498 / n2
pp <- (n1 * p1 + n2 * p2) / (n1 + n2)
se <- sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
z  <- (p1 - p2) / se
pval <- 2 * (1 - pnorm(abs(z)))
round(c(p1=p1, p2=p2, pooled=pp, se=se, z=z, pval=pval), 4)

## ----pearsonchi----------------------------------------------------------
x <- matrix(c(2537, 8747, 683, 1498), nrow=2, byrow=TRUE)
x
chisq.test(x, correct=FALSE)
# Also compute more accurate P-value based on 1M Monte-Carlo simulations
chisq.test(x, correct=FALSE, simulate.p.value=TRUE, B=1e6)

## ----bpower--------------------------------------------------------------
require(Hmisc)
bpower(.5, .7, n1=100, n2=100)

## ----bsamsize------------------------------------------------------------
bsamsize(.0015, 0.8 * .0015, alpha=0.05, power=0.8)

## ----moeprop-------------------------------------------------------------
diff <- .05
qnorm(.975)^2 / 2 / (diff ^ 2)

## ----mmeor,w=7,h=6,cap='Multiplicative margin of error related to 0.95 confidence limits of an odds ratio, for varying $n$ and $p$ (different curves), assuming the unknown true probability in each group is no lower than $p$',scap='Multiplicative margin of error for odds ratios'----
require(ggplot2)
d <- expand.grid(n=c(seq(10, 1000, by=10), seq(1100, 50000, by=100)),
                 p=c(.02, .05, .075, .1, .15, .2, .25, .3, .4, .5))
d$selor <- with(d, sqrt(2 / (p * (1 - p) * n)))
d$mmoe  <- with(d, exp(qnorm(0.975) * selor))
mb <- c(1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 10, 20, 30, 40, 50, 100, 400)
ggplot(aes(x=n, y=mmoe, color=factor(p)), data=d) +   # Fig. (*\ref{fig:prop-mmeor}*)
  geom_line() +
  scale_x_log10(breaks=c(10,20,30,50,100,200,500,1000,2000,5000,10000,
                  20000,50000)) +
  scale_y_log10(breaks=mb, labels=as.character(mb)) +
  xlab(expression(n)) + ylab('Multiplicative Margin of Error for OR') +
  guides(color=guide_legend(title=expression(p)))

## ----surgss--------------------------------------------------------------
round(bsamsize(.3, .1, fraction=1/11, power=.9))
round(bsamsize(.4, .2, fraction=1/11, power=.9))
round(bsamsize(.7, .9, fraction=1/11, power=.9))

## ----surgdata------------------------------------------------------------
n1 <- 25;     n2 <- 111
p1 <- 6 / n1; p2 <- 11 / n2
or <- p1 / (1 - p1) / (p2 / (1 - p2))
or
# Standard error of log odds ratio:
selor <- sqrt(1 / (n1 * p1 * (1 - p1)) + 1 / (n2 * p2 * (1 - p2)))
# Get 0.95 confidence limits
cls <- exp(log(or) + c(-1, 1) * qnorm(0.975) * selor)
cls
tcls <- paste0(round(or, 2), ' (0.95 CI: [', round(cls[1], 2),
               ', ', round(cls[2], 2), '])')
# Multiplying a constant by the vector -1, 1 does +/-
x <- matrix(c(6, 19, 11, 100), nrow=2, byrow=TRUE)
x
chisq.test(x, correct=FALSE)

## ----fet-----------------------------------------------------------------
fisher.test(x)

