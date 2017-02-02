## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('ancova', width=80)

## ----gustohistrisk,cap='Distribution of baseline risk in GUSTO-I.  Kernel density estimate of risk distribution for SK treatment.  Average risk is 0.07.  See also \\cite{ioa97imp}.',scap='Distribution of baseline risk in GUSTO-I'----
load('gustomin.rda')
with(gustomin,
     plot(density(p.sk), xlim=c(0, .4), xlab='Baseline Expected Risk',
          ylab='Probability Density', main='') )    # Fig. (*\ref{fig:ancova-gustohistrisk}*)

## ----or-diff,w=5.5,h=4,cap='Absolute risk reduction as a function of risk for control subject.  Numbers on curves are treatment:control odds ratios.',scap='Absolute risk reduction as a function of risk'----
plot(0, 0, type="n", xlab="Risk for Subject Without Risk Factor",
     ylab="Increase in Risk",
     xlim=c(0,1), ylim=c(0,.6))   # Figure (*\ref{fig:ancova-or-diff}*)
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
d <- expand.grid(X=0:1, brisk=seq(0.05, 0.3, length=150))
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
  ylab('Absolute Risk Reduction With New Treatment')

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
print(f, latex=TRUE)
latex(anova(f), file='', table.env=FALSE)
cat('}')
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
     xlim=c(0,1), ylim=c(0,.7))     # Fig. (*\ref{fig:ancova-hr-vs-surv}*)
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
     ylab='Probability Density', xlim=c(0, 6))    # Fig. (*\ref{fig:ancova-gusto-histcost}*)
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

