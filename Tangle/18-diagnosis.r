## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('dx', width=80)

## ----du,w=6,h=4.5,cap='Relative effect of total cholesterol for age 40 and 70; Data from Duke Cardiovascular Disease Databank, $n=2258$'----
require(rms)
getHdata(acath)
acath <- subset(acath, !is.na(choleste))
dd <- datadist(acath);  options(datadist='dd')
f <- lrm(sigdz ~ rcs(age,5)*sex, data=acath)
pre <- predict(f, type='fitted')
g <- lrm(sigdz ~ rcs(age,4)*sex + rcs(choleste,4) + rcs(age,4) %ia%
         rcs(choleste,4), data=acath)
ageg <- c(40, 70)
psig <- Predict(g, choleste, age=ageg)
s <- lrm(tvdlm ~ rcs(age,4)*sex + rcs(choleste,4) + rcs(age,4) %ia%
	rcs(choleste,4), data=acath)
psev <- Predict(s, choleste, age=ageg)
ggplot(rbind('Significant CAD'=psig, '3 Vessel or Left Main CAD'=psev),
	adj.subtitle=FALSE)

## ----du-acath,w=4.5,h=4.5, cap='Diagnostic Utility of Cholesterol for Diagnosing Significant CAD.  Curves are 0.1 and 0.9 quantiles from quantile regression using restricted cubic splines'----
post <- predict(g, type='fitted')
plot(pre, post, xlab='Pre-Test Probability (age + sex)',
     ylab='Post-Test Probability (age + sex + cholesterol)', pch=46)
abline(a=0, b=1, col=gray(.8))
lo <- Rq(post ~ rcs(pre, 7), tau=0.1)  # 0.1 quantile
hi <- Rq(post ~ rcs(pre, 7), tau=0.9)  # 0.9 quantile
at <- seq(0, 1, length=200)
lines(at, Predict(lo, pre=at)$yhat, col='red', lwd=1.5)
lines(at, Predict(hi, pre=at)$yhat, col='red', lwd=1.5)
abline(v=.5, col='red')

