## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('descript')

## ----pdfcdf,w=6,h=3,cap='Example probability density (a) and cumulative probability distribution (b) for a positively skewed random variable (skewed to the right)',scap='Density and cumulative distribution functions'----
x <- seq(-3, 35, length=150)
par(mfrow=c(1,2)); xl <- expression(x)   # Fig. (*\ref{fig:descript-pdfcdf}*):
plot(x, dt(x, 4, 6), type='l', xlab=xl, ylab='Probability Density Function')
plot(x, pt(x, 4, 6), type='l', xlab=xl, ylab='Cumulative Distribution Function')

## ----normalhist,cap='Example of a continuous distribution that is symmetric: the Gaussian (normal) distribution with mean 0 and variance 1, along with a histogram from a sample of size 1000 from this distribution',scap='Symmetric continuous distribution'----
set.seed(1); x <- rnorm(1000)   # Fig. (*\ref{fig:descript-normalhist}*):
hist(x, nclass=40, prob=TRUE, col=gray(.9), xlab=xl, ylab='')
x <- seq(-4, 4, length=150)
lines(x, dnorm(x, 0, 1), col='blue', lwd=2)

## ----bimode,cap='Example of a bimodal distribution from sampling from a mixture of normal distributions with different means and variances and estimating the underlying density function.  Vertical red lines indicate true population means of the two component populations.  Such a distribution can occur naturally or by failing to condition on a binary characteristic such as sex.',scap='Bimodal distribution'----
set.seed(2)
x <- c(rnorm(500, mean=0, sd=1), rnorm(500, mean=6, sd=3))
hist(x, nclass=40, prob=TRUE, col=gray(.9), xlab=xl, ylab='')
lines(density(x), col='blue', lwd=2)
abline(v=c(0, 6), col='red')   # Fig. (*\ref{fig:descript-bimode}*)

## ----orda,cap='Distribution of number of days in the hospital in the year following diagnosis',scap='Count variable with clumping at zero'----
x <- 0:14
y <- c(.8, .04, .03, .02, rep(.01, 11))
plot(x, y, xlab=xl, ylab='', type='n')   # Fig. (*\ref{fig:descript-orda}*)
segments(x, 0, x, y)

## ----ordb,cap='Distribution of a functional status score that does not have points in the middle',scap='Ordinal variable with strange distribution'----
x <- 1:10
y <- c(.1, .13, .18, .19, 0, 0, .14, .12, .08, .06)
plot(x, y, xlab=xl, ylab='', type='n')   # Fig. (*\ref{fig:descript-ordb}*)
segments(x, 0, x, y)

## ----ordc,cap='Distribution of serum creatinine where the patient requiring dialysis is taken to have the worst renal function.  The variable is mostly continuous but is best analyzed as ordinal so that no assumption is made about how to score dialysis except for being worse than all non-dialysis patients. Data taken from NHANES where no patients were actually dialyzed.',scap='Continuous distribution with clumping at the end'----
require(Hmisc)
getHdata(nhgh)   # NHANES dataset    Fig. (*\ref{fig:descript-ordc}*):
scr <- pmin(nhgh$SCr, 5)   # truncate at 5 for illustration
scr[scr == 5 | runif(nrow(nhgh)) < .05] <- 5  # pretend 1/20 dialyzed
hist(scr, nclass=50, xlab='Serum Creatinine', ylab='Density', prob=TRUE)

## ----spaghetti,w=7,h=5,cap='Spaghetti plot showing all the raw data on the response variable for each subject, stratified by dose and study site (1--9).  Importantly, week 0 (baseline) measurements are included.',scap='Spaghetti plot'----
require(Hmisc)   # also loads ggplot2
getHdata(cdystonia)
ggplot(cdystonia, aes(x=week, y=twstrs, color=factor(id))) +
       geom_line() + xlab('Week') + ylab('TWSTRS-total score') +
       facet_grid(treat ~ site) +
       guides(color=FALSE) # Fig. (*\ref{fig:descript-spaghetti}*)

## ----pH,w=5,h=4,cap='Scatterplot of one measurement mode against another'----
getHdata(esopH)
contents(esopH)
xl <- label(esopH$conv)
yl <- label(esopH$orophar)
ggplot(esopH, aes(x=conv, y=orophar)) + geom_point(pch='.') +
  xlab(xl) + ylab(yl) +   # Fig. (*\ref{fig:descript-pH}*)
  geom_abline(intercept = 0, slope = 1)

## ----pHh,w=5.5,h=4.5,cap='Hexagonal binning replacing scatterplot for large $n$'----
require(hexbin)
ggplot(esopH, aes(x=conv, y=orophar)) +   # Fig. (*\ref{fig:descript-pHh}*)
  stat_binhex(aes(alpha=..count.., color=Hmisc::cut2(..count.., g=20)),
              bins=80, range=c(.1, .9)) +
  xlab(xl) + ylab(yl) +
  guides(alpha=FALSE, fill=FALSE, color=guide_legend(title='Frequency'))

## ----ecdf,w=6,h=4.5,cap='Empirical cumulative distributions of baseline variables  stratified by treatment in a randomized controlled trial. $m$ is the number of missing values.',scap='Empirical cumulative distribution functions'----
getHdata(pbc)
pbcr <- subset(pbc, drug != 'not randomized')
Ecdf(pbcr[,c('bili','albumin','protime','sgot')],  # Fig. (*\ref{fig:descript-ecdf}*)
    group=pbcr$drug, col=1:2,
    label.curves=list(keys='lines'))

## ----bwplot,w=5,h=3.25,cap='Box plots showing the distribution of serum creatinine  stratified by major diagnosis.  Dot: mean; vertical line: median; large box:interquartile range.  The 0.05 and 0.95 quantiles are also shown, which is not the way typical box plots are drawn but is perhaps more useful.  Asymmetry of distributions can be seen by both disagreement between $Q_{3}-Q_{2}$ and $Q_{2}-Q_{1}$ and by disagreement between $Q_{2}$ and $\\bar{x}$.',scap='Box plots with 0.05 and 0.95 quantiles'----
getHdata(support)   # Fig. (*\ref{fig:descript-bwplot}*)
bwplot(dzgroup ~ crea, data=support, panel=panel.bpplot,
       probs=c(.05,.25), xlim=c(0,8), xlab='Serum Creatinine')

## ----bpplt,w=5.5,h=3.75,cap='Schematic for extended box plot'------------
bpplt()   # Fig. (*\ref{fig:descript-bpplt}*)

## ----glybox,h=5,w=6,cap='Extended box plots for glycohemoglobin stratified by quartiles of age (vertical), two-tiles of waist circumference (horizontal), and sex (vertical)',scap='Extended box plots for glycohemoglobin'----
require(lattice)   # Fig. (*\ref{fig:descript-glybox}*):
getHdata(diabetes)
wst <- cut2(diabetes$waist, g=2)
levels(wst) <- paste('Waist', levels(wst))
bwplot(cut2(age,g=4) ~ glyhb | wst*gender, data=diabetes,
       panel=panel.bpplot, xlab='Glycosylated Hemoglobin', ylab='Age Quartile')

## ----cldif,cap='Means and nonparametric bootstrap 0.95 confidence limits for glycated hemoglobin for males and females, and confidence limits for males - females.  Lower and upper $x$-axis scales have same spacings but different centers.  Confidence intervals for differences are generally wider than those for the individual constituent variables.',scap='Showing group means and differences',mar=c(4,6,4,1),top=1----
attach(diabetes)
set.seed(1)
male   <- smean.cl.boot(glyhb[gender=='male'],   reps=TRUE)
female <- smean.cl.boot(glyhb[gender=='female'], reps=TRUE)
dif <- c(mean=male['Mean']-female['Mean'],
         quantile(attr(male, 'reps')-attr(female,'reps'), c(.025,.975)))
plot(0,0,xlab='Glycated Hemoglobin',ylab='',   # Fig. (*\ref{fig:descript-cldif}*)
     xlim=c(5,6.5),ylim=c(0,4), axes=F)
axis(1, at=seq(5, 6.5, by=0.25))
axis(2, at=c(1,2,3.5), labels=c('Female','Male','Difference'),
     las=1, adj=1, lwd=0)
points(c(male[1],female[1]), 2:1)
segments(female[2], 1, female[3], 1)
segments(male[2], 2,   male[3], 2)
offset <- mean(c(male[1],female[1])) - dif[1]
points(dif[1] + offset, 3.5)
segments(dif[2]+offset, 3.5, dif[3]+offset, 3.5)
at <- c(-.5,-.25,0,.25,.5,.75,1)
axis(3, at=at+offset, label=format(at))
segments(offset, 3, offset, 4.25, col=gray(.85))
abline(h=c(2 + 3.5)/2, col=gray(.85))

## ----rmsa,results='asis',w=8,h=5,cap='Partial effects in NHANES HbA1c model'----
require(rms)
getHdata(nhgh)   # NHANES data
dd <- datadist(nhgh); options(datadist='dd')
g        <- function(x) 0.09 - x ^ - (1 / 1.75)
ginverse <- function(y) (0.09 - y) ^ -1.75
f <- ols(g(gh) ~ rcs(age, 4) + re + sex + rcs(bmi, 4), data=nhgh)
# If not using LaTeX, just use print(f), anova(f)
cat('{\\small\n')
print(f, latex=TRUE)
latex(anova(f), dec.ss=3, dec.ms=3, file='', table.env=FALSE)
cat('}\n')
# Show partial effects of all variables in the model, on the original scale
ggplot(Predict(f, fun=ginverse),   # Fig. (*\ref{fig:descript-rmsa}*)
       ylab=expression(paste('Predicted Median ', HbA['1c'])))

## ----rmsb,w=7,h=3,top=2,cap='Partial effects chart on the transformed scale.  For age and BMI, effects are inter-quartile-range effects.  0.9, 0.95, and 0.99 confidence limits are shown.',scap='Partial effects chart for transformed glycohemoglobin'----
plot(summary(f))   # Fig. (*\ref{fig:descript-rmsb}*)

## ----rmsc,w=5,h=5,cap='Nomogram for predicting median HbA$_{1c}$.  To use the nomogram, use the top \\Co{Points} scale to convert each predictor value to a common scale.  Add the points and read this number on the \\Co{Total Points} scale, then read down to the median.',scap='Nomogram for predicting median HbA$_{1c}$'----
plot(nomogram(f, fun=ginverse, funlabel='Median HbA1c'))  # Fig. (*\ref{fig:descript-rmsc}*)

## ----dynamite,w=5.5,cap="Bar plot with error bars---``dynamite plot''"----
getHdata(FEV); set.seed(13)   
FEV <- subset(FEV, runif(nrow(FEV)) < 1/8)   # 1/8 sample
require(ggplot2)
s <- with(FEV, summarize(fev, llist(sex, smoke), smean.cl.normal))
ggplot(s, aes(x=smoke, y=fev, fill=sex)) +    # Fig. (*\ref{fig:descript-dynamite}*)
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  width=.1,
                  position=position_dodge(.9))

## ----tplot,w=6,h=2,cap='Jittered raw data and box plots.  Middle vertical lines indicate medians and diamonds indicate means. Horizontal lines indicate 0.1 to 0.9 quantiles when $n\\geq 10$.  The ink:information ratio for this plot is far better than a dynamite plot.',scap='Dot plot with superimposed box plots'----
require(ggplot2)   # Fig. (*\ref{fig:descript-tplot}*)
stats <- function(x) {
  z <- quantile(x, probs=c(.1, .25, .5, .75, .9))
  names(z) <- c('ymin', 'lower', 'middle', 'upper', 'ymax')
  if(length(x) < 10) z[c(1,5)] <- NA
  z
}
ggplot(FEV, aes(x=sex, y=fev)) +
  stat_summary(fun.data=stats, geom='boxplot', aes(width=.75), shape=5,
               position='dodge', col='lightblue') +
  geom_dotplot(binaxis='y', stackdir='center', position='dodge', alpha=.4) +
  stat_summary(fun.y=mean, geom='point', shape=5, size=4, color='blue') +
  facet_grid(~ smoke) +
  xlab('') + ylab(expression(FEV[1])) + coord_flip()

## ----vplot,w=6,h=2.5,cap='Jittered raw data and violin plots with median indicated by \\textbf{\\textcolor{blue}{+}}'----
ggplot(FEV, aes(x=sex, y=fev)) +
  geom_violin(width=.6, col='lightblue') +
  geom_dotplot(binaxis='y', stackdir='center', position='dodge', alpha=.4) +
  stat_summary(fun.y=median, geom='point', color='blue', shape='+', size=12) +
  facet_grid(~ smoke) +
  xlab('') + ylab(expression(FEV[1])) + coord_flip()

