## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('ancova', width=80)
options(prType='latex')

## ----gustohistrisk,cap='Distribution of baseline risk in GUSTO-I.  Kernel density estimate of risk distribution for SK treatment.  Average risk is 0.07.  See also \\cite{ioa97imp}.',scap='Distribution of baseline risk in GUSTO-I'----
load('gustomin.rda')
with(gustomin,
     plot(density(p.sk), xlim=c(0, .4), xlab='Baseline Expected Risk',
          ylab='Probability Density', main='') )    # Fig. (*\ref{fig:ancova-gustohistrisk}*)

## ----htesim,results='asis'-----------------------------------------------
require(rms)
options(prType='latex')   # for cph print, anova
set.seed(1)
n <- 3000    # total of 3000 subjects
age <- rnorm(n, 60, 12)
label(age) <- 'Age'
sex   <- factor(sample(c('Male', 'Female'), n, rep=TRUE))
treat <- factor(sample(c('A', 'B'), n, rep=TRUE))
cens  <- 15 * runif(n)     # censoring time
h <- 0.02 * exp(0.04 * (age - 60) + 0.4 * (sex == 'Female') -  
                0.04 * (treat == 'B') * pmax(age - 60, 0))
dt <- -log(runif(n)) / h
label(dt) <- 'Time Until Death or Censoring'
e <- ifelse(dt <= cens, 1, 0)
dt <- pmin(dt, cens)
units(dt) <- 'Year'
dd <- datadist(age, sex, treat); options(datadist='dd')
S <- Surv(dt, e)
f <- cph(S ~ sex + rcs(age, 4) * treat)
f
anova(f)

## ----hteplot-------------------------------------------------------------
ggplot(Predict(f, age, treat), rdata=data.frame(age, treat))
ages <- seq(30, 87, length=200)
k <- contrast(f, list(treat='B', age=ages), list(treat='A', age=ages))
class(k) <- 'data.frame'
ggplot(k, aes(x=age, y=exp(Contrast))) + scale_y_log10(minor_breaks=seq(.2, .9, by=.1)) + 
  geom_ribbon(aes(ymin=exp(Lower), ymax=exp(Upper)), fill='gray80') +
  geom_line() +
  ylab('B:A Hazard Ratio') + xlab('Age')

## ----hteplot2,results='asis'---------------------------------------------
g <- cph(S ~ sex * treat + rcs(age, 4))
g
anova(g)
ggplot(Predict(g, sex, treat))
k <- contrast(g, list(treat='B', sex=levels(sex)), list(treat='A', sex=levels(sex)))
class(k) <- 'data.frame'
ggplot(k, aes(y=exp(Contrast), x=sex)) + geom_point() + scale_y_log10(breaks=c(.5, .6, .7, .8, .9, 1, 1.1)) +
  geom_linerange(aes(ymin=exp(Lower), ymax=exp(Upper))) + 
  xlab('Sex') + ylab('B:A Hazard Ratio') + coord_flip()

## ----htejoint,results='asis'---------------------------------------------
h <- cph(S ~ treat * (rcs(age, 4) + sex))
anova(h)
ggplot(Predict(h, sex, treat, age=50))
k <- contrast(h, list(treat='B', sex=levels(sex), age=50),
                 list(treat='A', sex=levels(sex), age=50))
class(k) <- 'data.frame'
ggplot(k, aes(y=exp(Contrast), x=sex)) + geom_point() + scale_y_log10(breaks=c(.5, .6, .7, .8, .9, 1, 1.1)) +
  geom_linerange(aes(ymin=exp(Lower), ymax=exp(Upper))) + 
  xlab('Sex') + ylab('B:A Hazard Ratio') + coord_flip()

## ----or-diff,w=5.5,h=4,cap='Absolute risk increase as a function of risk for control subject.  Numbers on curves are treatment:control odds ratios.',scap='Absolute risk increase as a function of risk'----
plot(0, 0, type="n", xlab="Risk for Subject Without Risk Factor",
     ylab="Increase in Risk",
     xlim=c(0,1), ylim=c(0,.6))   # Figure (*\ref{fig:ancova-or-diff}\ipacue*)
i <- 0
or <- c(1.1,1.25,1.5,1.75,2,3,4,5,10)
for(h in or) {
  i <- i + 1
  p <- seq(.0001, .9999, length=200)
  logit <- log(p/(1 - p))  # same as qlogis(p)
  logit <- logit + log(h)  # modify by odds ratio
  p2 <- 1/(1 + exp(-logit))# same as plogis(logit)
  d <- p2 - p
  lines(p, d, lty=i)
  maxd <- max(d)
  smax <- p[d==maxd]
  text(smax, maxd + .02, format(h), cex=.6)
}

## ----ordiff,h=4,w=5.5,cap='Absolute risk reduction by a new treatment as a function of background risk and an interacting factor',scap='Absolute risk reduction by background risk and interacting factor'----
require(Hmisc)
d <- expand.grid(X=0:1, brisk=seq(0.01, 0.99, length=150))
d <- upData(d,
            risk.standard = plogis(qlogis(brisk) + log(1.4) * X),
            risk.new      = plogis(qlogis(brisk) + log(1.4) * X +
                                     log(0.8) * (X == 0) +
                                     log(0.6) * (X == 1)),
            risk.diff     = risk.standard - risk.new,
            X = factor(X) )
ggplot(d, aes(x=risk.standard, y=risk.diff, color=X)) +
  geom_line() +
  xlab('Risk Under Standard Treatment') +
  ylab('Absolute Risk Reduction With New Treatment')   # (*\ipacue*)

## ----gusto-histdelt,w=5,h=4,cap='Distribution of absolute risk reduction with $t$--PA vs.\\ SK',scap='Absolute benefit vs.\\ baseline risk'----
delta <- with(gustomin, p.sk - p.tpa)
plot(density(delta), xlab='Mortality Difference',
     ylab='Probability Density', main='')    # Fig. (*\ref{fig:ancova-gusto-histdelt}*)
m <- mean(delta)
u <- par("usr")
arrows(m, u[3], m, 0, length=.1, lwd=2)
text(m, 2, 'Mean', srt=45, adj=0)
med <- median(delta)
arrows(med, u[3], med, 0, length=.1, lwd=2)
text(med, 2, 'Median', srt=45, adj=0)

## ----gusto-nomogram,results='asis',w=7,h=6,cap='Nomogram to predict SK - $t$--PA mortality difference, based on the difference between two binary logistic models.',scap='GUSTO-I nomogram'----
load('gusto.rda')
require(rms)
dd <- datadist(gusto); options(datadist='dd')
f <- lrm(day30 ~ tx + age * Killip + pmin(sysbp, 120) +
           lsp(pulse, 50) + pmi + miloc, data=gusto)
cat('{\\smaller ')
f
anova(f)
cat('}')   # (*\ipacue*)
cof <- coef(f)             # vector of regression coefficients
# For cof, X*beta without treatment coefficients estimates logit
# for SK+t-PA combination therapy (reference cell).  The coefficient for
# SK estimates the difference in logits from combo to SK.  The coefficient
# for tPA estimates the difference in tPA from combo.  The mortality
# difference of interest is mortality with SK minus mortality with tPA.
mort.sk   <- function(x) plogis(x + cof['tx=SK'])
mort.diff <- function(x)
  ifelse(x < 0, mort.sk(x) - plogis(x + cof['tx=tPA']), NA)
# only define when logit < 0 since U-shaped
n <- nomogram(f, fun=list(mort.sk, mort.diff),
              funlabel=c("30-Day Mortality\nFor SK Treatment",
                         "Mortality Reduction by t-PA"),
              fun.at=list(c(.001,.005,.01,.05,.1,.2,.5,.7,.9),
                c(.001,.005,.01,.02,.03,.04,.05)),
              pulse=seq(0,260,by=10), omit='tx', lp=FALSE)
plot(n, varname.label.sep=' ', xfrac=.27, lmgp=.2, cex.axis=.6)

## ----hr-vs-surv,w=5.5,h=4,cap='Relationship between baseline risk, relative treatment effect (hazard ratio --- numbers above curves) and absolute treatment effect.',scap='Baseline risk, hazard ratio, and absolute effect'----
plot(0, 0, type="n", xlab="Survival for Control Subject",
     ylab="Improvement in Survival",
     xlim=c(0,1), ylim=c(0,.7))     # Fig. (*\ref{fig:ancova-hr-vs-surv}\ipacue*)
i <- 0
hr <- seq(.1, .9, by=.1)
for(h in hr) {
  i <- i + 1
  p <- seq(.0001, .9999, length=200)
  p2 <- p ^ h
  d <- p2 - p
  lines(p, d, lty=i)
  maxd <- max(d)
  smax <- p[d==maxd]
  text(smax,maxd+.02, format(h), cex=.6)
}

## ----gusto-histcost,w=5,h=4,cap='Distribution of cost per life saved in GUSTO--I'----
cost.life <- 2400 / delta / 1e6
plot(density(cost.life), xlab='Cost Per Life Saved, $M', main='',
     ylab='Probability Density', xlim=c(0, 6))    # Fig. (*\ref{fig:ancova-gusto-histcost}\ipacue*)
m <- 2400 / mean(delta) / 1e6
u <- par("usr")
arrows(m, u[3], m, 0, length=.1, lwd=2)
text(m,.01,'Cost using\n  average\n    reduction',srt=45,adj=0)

## ----contrast,eval=FALSE-------------------------------------------------
## sites <- levels(site)
## contrast(fit, list(treat='b', site=sites),
##               list(treat='a', site=sites),
##          type='average', weights=table(site))

## ----contrastb,eval=FALSE------------------------------------------------
## contrast(fit, list(treat='b', site=sites),
##               list(treat='a', site=sites), type='average')

