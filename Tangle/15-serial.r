## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('serial', width=80)

## ----download------------------------------------------------------------
require(Hmisc)
require(data.table)   # elegant handling of aggregation
require(ggplot2)
d <- csv.get(paste('http://biostat.mc.vanderbilt.edu',
                   'dupontwd/wddtext/data/11.2.Long.Isoproterenol.csv',
                   sep='/'))
d <- upData(d, keep=c('id', 'dose', 'race', 'fbf'),
            race  =factor(race, 1:2, c('white', 'black')),
            labels=c(fbf='Forearm Blood Flow'),
            units=c(fbf='ml/min/dl'))
d <- data.table(d)
setkey(d, id, race)

## ----spag,w=7,h=4,cap="Spaghetti plots for isoproterenol data showing raw data stratified by race.  Next to each curve is the area under the curve divided by 400 to estimate the mean response function.  The area is computed analytically from a restricted cubic spline function fitted separately to each subject's dose-response curve.  Shadowing the raw data are faint lines depicting spline fits for each subject",scap='Spaghetti plots for isoproterenol data'----
# Fit subject-by-subject spline fits and either return the coefficients,
# the estimated area under the curve from [0,400], or evaluate each
# subject's fitted curve over a regular grid of 150 doses
# Area under curve is divided by 400 to get a mean function
require(rms)
g <- function(x, y, what=c('curve', 'coef', 'area')) {
  what <- match.arg(what)   # 'curve' is default
  knots <- c(20, 60, 150)
  f <- ols(y ~ rcs(x, knots))
  xs <- seq(0, 400, length=150)
  switch(what,
         coef = {k <- coef(f)
                 list(b0 = k[1], b1=k[2], b2=k[3])},
         curve= {x <- seq(0, 400, length=150)
                 list(dose=xs, fbf=predict(f, data.frame(x=xs)))},
         area = {antiDeriv = rcsplineFunction(knots, coef(f),
                   type='integral')
                 list(dose  = 400, fbf=y[x == 400],
                      area  = antiDeriv(400) / 400,
                      tarea = areat(x, y) / 400)} )
}
# Function to use trapezoidal rule to compute area under the curve
areat <- function(x, y) {
  i <- ! is.na(x + y)
  x <- x[i]; y <- y[i]
  i <- order(x)
  x <- x[i]; y <- y[i]
  if(! any(x == 400)) NA else
  sum(diff(x) * (y[-1] + y[-length(y)]))/2
}

w <- d[, j=g(dose, fbf), by = list(id, race)]   # uses data.table package
a <- d[, j=g(dose, fbf, what='area'), by = list(id, race)]

ggplot(d, aes(x=dose, y=fbf, color=factor(id))) +   # Fig. (*\ref{fig:serial-spag}\ipacue*)
       geom_line() + geom_line(data=w, alpha=0.25) +
       geom_text(aes(label = round(area,1)), data=a, size=2.5,
                 position=position_dodge(width=50)) +
       xlab('Dose') + ylab(label(d$fbf, units=TRUE, plot=TRUE)) +
       facet_grid(~ race) +
       guides(color=FALSE)

## ----auctwo,cap='AUC by curve fitting and by trapezoidal rule'-----------
ggplot(a, aes(x=tarea, y=area, color=race)) + geom_point() +
  geom_abline(col=gray(.8)) +
  xlab('Area by Trapezoidal Rule / 400') +
  ylab('Area by Spline Fit / 400')          # Fig. (*\ref{fig:serial-auctwo}\ipacue*)

## ----auc,w=5,h=1.5,cap='Mean blood flow computed from the areas under the spline curves, stratified by race, along with box plots',scap='Mean blood flow by race'----
ggplot(a, aes(x=race, y=area)) +    # Fig. (*\ref{fig:serial-auc}\ipacue*)
  geom_boxplot(alpha=.5, width=.25) + geom_point() + coord_flip() +
  ylab(expression(paste('Mean Forearm Blood Flow,  ', scriptstyle(ml/min/dl))))

## ----wil-----------------------------------------------------------------
wilcox.test(area ~ race, data=a, conf.int=TRUE)

## ----coefs---------------------------------------------------------------
h <- d[, j=g(dose, fbf, what='coef'), by = list(id, race)]
h

## ----lrm,results='asis'--------------------------------------------------
f <- lrm(race ~ b0 + b1 + b2, data=h, x=TRUE, y=TRUE)
print(f, latex=TRUE)

## ----boot,results='asis'-------------------------------------------------
set.seed(2)
v <- validate(f, B=1000)
latex(v, file='')

## ----glsa,cap='Residual plot for generalized least squares fit on untransformed  \\co{fbf}'----
require(nlme)
dd <- datadist(d); options(datadist='dd')
a <- Gls(fbf ~ race * rcs(dose, c(20,60,150)), data=d,
         correlation=corCAR1(form = ~ dose | id))
plot(fitted(a), resid(a))   # Fig. (*\ref{fig:serial-glsa}\ipacue*)

## ----glsb,results='asis'-------------------------------------------------
a <- Gls(log(fbf) ~ race * rcs(log(dose + 1), log(c(20,60,150)+1)), data=d,
         correlation=corCAR1(form = ~ dose | id))
latex(anova(a), file='', table.env=FALSE)

## ----glstwocorr----------------------------------------------------------
a <- Gls(log(fbf) ~ race * log(dose + 1), data=d,
         correlation=corCAR1(form = ~ dose | id))
d$time <- match(d$dose, c(0, 10, 20, 60, 150, 300, 400)) - 1
b <- Gls(log(fbf) ~ race * log(dose + 1), data=d,
         correlation=corCAR1(form = ~ time | id))
AIC(a);AIC(b)

## ----glsc,results='asis',w=6,h=4,cap='Checking assumptions of the GLS model that is linear after logging dose and blood flow.  The graph on the right is a QQ-plot to check normality of the residuals from the model, where linearity implies normality.',scap='Checking assumptions of GLS model'----
print(b, latex=TRUE)
latex(anova(b), file='', table.env=FALSE)
w <- data.frame(residual=resid(b), fitted=fitted(b))
p1 <- ggplot(w, aes(x=fitted, y=residual)) + geom_point()
p2 <- ggplot(w, aes(sample=residual)) + stat_qq()
gridExtra::grid.arrange(p1, p2, ncol=2)   # Figure (*\ref{fig:serial-glsc}\ipacue*)

## ----glsd,w=6,h=3.75,cap='Pointwise and simultaneous confidence bands for median dose-response curves by race',scap='Confidence bands for median dose-response curves',cache=TRUE----
dos <- seq(0, 400, length=150)
p <- Predict(b, dose=dos, race, fun=exp)
s <- Predict(b, dose=dos, race, fun=exp, conf.type='simultaneous')
ps <- rbind(Pointwise=p, Simultaneous=s)
ggplot(ps, ylab=expression(paste('Median Forearm Blood Flow,  ',
                          scriptstyle(ml/min/dl))))   # Fig. (*\ref{fig:serial-glsd}\ipacue*)

## ----glse,cap='White:black fold change for median response as a function of dose, with simultaneous confidence band',cache=TRUE----
k <- contrast(b, list(dose=dos, race='white'),
                 list(dose=dos, race='black'), conf.type='simultaneous')
k <- as.data.frame(k[c('dose', 'Contrast', 'Lower', 'Upper')])
ggplot(k, aes(x=dose, y=exp(Contrast))) + geom_line() +
  geom_ribbon(aes(ymin=exp(Lower), ymax=exp(Upper)), alpha=0.2, linetype=0,
              show_guide=FALSE) +
  geom_hline(yintercept=1, col='red', size=.2) +
  ylab('White:Black Ratio of Median FBF')    # Fig. (*\ref{fig:serial-glse}\ipacue*)

