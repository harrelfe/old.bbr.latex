## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('change', width=80)

## ----diabetes,w=6,h=2.5,cap='$\\beta$-TG levels by diabetic status with a median line.  The left plot is on the original (non-transformed) scale and includes median lines.  The right plot displays the data on a log scale.',scap='$\\beta$-TG levels by diabetic status'----
d <- rbind(
  data.frame(status='normal',
             btg=c(4.1, 6.3, 7.8, 8.5, 8.9, 10.4, 11.5, 12.0, 13.8,
                   17.6, 24.3, 37.2)),
  data.frame(status='diabetic',
             btg=c(11.5, 12.1, 16.1, 17.8, 24.0, 28.8, 33.9, 40.7,
                   51.3, 56.2, 61.7, 69.2)))
require(ggplot2)
p1 <- ggplot(d, aes(x=status, y=btg)) +    # Fig. (*\ref{fig:change-diabetes}*)
  geom_dotplot(binaxis='y', stackdir='center', position='dodge') +
  geom_errorbar(stat = "hline", yintercept = "median", width=.25, size=1.3,
                aes(ymin=..y.., ymax=..y..)) +
  xlab('') + ylab(expression(paste(beta-TG, ' (ng/day/100 ml creatinine)'))) +
  coord_flip()
p2 <- ggplot(d, aes(x=status, y=btg)) +
  scale_y_log10(breaks=c(4,5,10,15,20,30,40,60,80)) +
  geom_dotplot(binaxis='y', stackdir='center', position='dodge') +
  xlab('') + ylab(expression(paste(beta-TG, ' (ng/day/100 ml creatinine)'))) +
  coord_flip()
arrGrob(p1, p2, ncol=2)

## ----diabetest-----------------------------------------------------------
t.test(btg ~ status, data=d)
t.test(log(btg) ~ status, data=d)

## ----diabetesw-----------------------------------------------------------
wilcox.test(btg ~ status, data=d)
wilcox.test(log(btg) ~ status, data=d)

## ----suppcr,w=5,h=3.75,cap='Estimated risk of hospital death as a function of day 3 serum creatinine and sex for 7772 critically ill ICU patients having day 1 serum creatinine $< 2$ and surviving to the start of day 3 in the ICU',scap='Hospital death as a function of creatinine'----
require(rms)
load('~/Analyses/SUPPORT/combined.sav')
combined <- subset(combined,
  select=c(id, death, d.time, hospdead, dzgroup, age, raceh, sex))
load('~/Analyses/SUPPORT/combphys.sav')
combphys <- subset(combphys, !is.na(crea1+crea3),
                   select=c(id,crea1,crea3,crea7,crea14,crea25,alb3,
                     meanbp3,pafi3,wblc3))
w <- merge(combined, combphys, by='id')
u <- 'mg/dl'
w <- upData(w, labels=c(crea1='Serum Creatinine, Day 1',
                 crea3='Serum Creatinine Day 3',
                 crea14='Serum Creatinine Day 14'),
            units=c(crea1=u, crea3=u, crea7=u, crea14=u, crea25=u))

w <- subset(w, crea1 < 2)
dd <- datadist(w); options(datadist='dd')

h <- lrm(hospdead ~ rcs(crea1, 5) + rcs(crea3, 5), data=w)
anova(h)   # (*\label{pg:change-anova}*)
h <- lrm(hospdead ~ sex * rcs(crea3, 5), data=w)
p <- Predict(h, crea3, sex, fun=plogis)
ggplot(p, ylab='Risk of Hospital Death')    # Fig. (*\ref{fig:change-suppcr}*)

## ----analysis,w=6,h=5,cap='Bland-Altman plots for three transformations'----
# Now add simulated some post data to the analysis of beta TG data
# Assume that the intervention effect (pre -> post effect) is
# multiplicative (x 1/4) and that there is a multiplicative error
# in the post measurements
set.seed(13)
d$pre  <- d$btg
d$post <- exp(log(d$pre) + log(.25) + rnorm(24, 0, .5))
# Make Bland-Altman (Tukey mean-difference) plots on the original and
# log scales
p1 <- ggplot(d, aes(x=(pre + post) / 2, y=post - pre, color=status)) +
  geom_point() + geom_smooth() + theme(legend.position='bottom')
# Use problematic asymmetric % change
p2 <- ggplot(d, aes(x=exp((log(pre) + log(post))/2), y=100*(post - pre)/pre,
                    color=status)) + geom_point() + geom_smooth() +
      xlab('Geometric Mean') + theme(legend.position='none') +
      ylim(-125, 0)
p3 <- ggplot(d, aes(x=exp((log(pre) + log(post))/2), y=log(post / pre),
                    color=status)) + geom_point() + geom_smooth() +
      xlab('Geometric Mean') + theme(legend.position='none') + ylim(-2.5, 0)
arrGrob(p1, p2, p3, ncol=2)   # Fig. (*\ref{fig:change-analysis}*)
with(d, {
     print(t.test(post - pre))
     print(t.test(100*(post - pre) / pre))       # improper
     print(t.test(log(post / pre)))
     print(wilcox.test(post - pre))
     print(wilcox.test(100*(post - pre) / pre))  # improper
     print(wilcox.test(log(post / pre)))
     } )

## ----baseball,w=5.5,h=4.5,bot=1,left=-3,rt=.5,cap='Initial batting averages as estimates of final batting averages for players, along with shrunken estimates that account for regression to the mean',scap='Baseball batting averages and regression to the mean'----
nam <- c('Roberto Clemente','Frank Robinson','Frank Howard','Jay Johnstone',
  	 'Ken Berry','Jim Spencer','Don Kessinger','Luis Alvarado',
		 'Ron Santo','Ron Swoboda','Del Unser','Billy Williams',
		 'George Scott','Rico Petrocelli','Ellie Rodriguez',
		 'Bert Campaneris','Thurman Munson','Max Alvis')
initial <- c(18,17,16,15,14,14,13,12,11,11,10,10,10,10,10,9,8,7)/45
season  <- c(345,297,275,220,272,270,265,210,270,230,265,258,306,265,225,
             283,320,200)/1000
initial.shrunk <- c(294,288,280,276,275,275,270,265,262,263,258,256,
                    257,256,257,252,245,240)/1000
plot(0,0,xlim=c(0,1),ylim=c(.15,.40),type='n',axes=F,xlab='',ylab='')
n  <- 18
x1 <- .5
x2 <- .75
x3 <- 1
points(rep(x1,n), initial)
points(rep(x2,n), initial.shrunk)
points(rep(x3,n), season)

for(i in 1:n) lines(c(x1,x2,x3),c(initial[i],initial.shrunk[i],season[i]),
                    col=i, lwd=2.5)
axis(2)
par(xpd=NA)
text(c(x1,x2+.01, x2+.25),rep(.12,3),c('First 45 ABs','Shrunken\nEstimates',
     'Rest of\nSeason'))
for(a in unique(initial)) {
  s <- initial==a
  w <- if(sum(s) < 4) paste(nam[s],collapse=', ') else {
	j <- (1:n)[s]
	paste(nam[j[1]],', ',nam[j[2]],', ',nam[j[3]],'\n',
		  nam[j[4]],', ',nam[j[5]],sep='')
  }
  text(x1-.02, a, w, adj=1, cex=.9)
}

