require(rms)
d <- csv.get('dat.txt', header=TRUE, sep='')
d <- upData(d, baseline=baseline/100,
            m3=m3/100, m6=m6/100, m9=m9/100, m12=m12/100,
            tx=c(rep('hf', 9), rep('lf', 10), rep('placebo', 8)))

base <- subset(d, select=c(ptno, baseline))

z <- reshape(d, idvar='ptno', varying=c('baseline', 'm3', 'm6', 'm9', 'm12'), v.names='cpeptide',
             timevar='week', times=c(0, 3, 6, 9, 12), direction='long')

i <- with(z, order(ptno, week))
z <- z[i, ]
ggplot(z, aes(x=week, y=cpeptide, group=factor(ptno), color=tx)) + geom_line()
ggplot(z, aes(x=week, y=log10(cpeptide), group=factor(ptno), color=tx)) + geom_line()

z <- subset(z, week > 0)
z <- merge(base, z, by='ptno')



require(nlme)
dd <- datadist(z); options(datadist='dd')
f <- Gls(log10(cpeptide) ~ pol(week, 2) * tx + log10(baseline),
         correlation=corCAR1(form = ~ week | ptno), data=z)
f
anova(f)

wks <- seq(3, 12, length=100)
ggplot(Predict(f, week=wks, tx))
plotp(Predict(f, week=wks, tx))
