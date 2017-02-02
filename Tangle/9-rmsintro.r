## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('rmsintro', width=80)

## ----eval=FALSE----------------------------------------------------------
## y ~ x1 + x2 + x3

## ----eval=FALSE----------------------------------------------------------
## f <- ols(y ~ age + sys.bp, data=mydata)

## ----eval=FALSE----------------------------------------------------------
## f$coefficients
## f$coef          # abbreviation
## coef(f)         # use the coef extractor function
## coef(f)[1]      # get intercept
## f$coef[2]       # get 2nd coefficient (1st slope)
## f$coef['age']   # get coefficient of age
## coef(f)['age']  # ditto

## ----eval=FALSE----------------------------------------------------------
## dd <- datadist(x1,x2,x3,...)   # generic form
## dd <- datadist(age, sys.bp, sex)
## dd <- datadist(mydataframe)    # for a whole data frame
## options(datadist='dd')         # let rms know where to find

## ----eval=FALSE----------------------------------------------------------
## dd <- datadist(dd, cholesterol, height)
## # Adds or replaces cholesterol, height summary stats in dd

## ----lead----------------------------------------------------------------
# For an Rmarkdown version of similar analyses see
# https://github.com/harrelfe/rscripts/raw/master/lead-ols.md
require(rms)    # also loads the Hmisc package
getHdata(lead)
# Subset variables just so contents() and describe() output is short
# Override units of measurement to make them legal R expressions
lead <- upData(lead,
               keep=c('ld72', 'ld73', 'age', 'maxfwt'),
               labels=c(age='Age'),
               units=c(age='years', ld72='mg/100*ml', ld73='mg/100*ml'))

contents(lead)
describe(lead)
dd <- datadist(lead); options(datadist='dd')
dd    # show what datadist computed
# Fit an ordinary linear regression model with 3 predictors assumed linear
f <- ols(maxfwt ~ age + ld72 + ld73, data=lead)
f         # same as print(f)
coef(f)   # retrieve coefficients
specs(f, long=TRUE)   # show how parameters are assigned to predictors,
                      # and predictor distribution summaries driving plots
g <- Function(f)  # create an R function that represents the fitted model
# Note that the default values for g's arguments are medians
g
# Estimate mean maxfwt at age 10, .1 quantiles of ld72, ld73 and .9 quantile of ld73
# keeping ld72 at .1 quantile
g(age=10, ld72=21, ld73=c(21, 47))  # more exposure in 1973 decreased y by 6
# Get the same estimates another way but also get std. errors
predict(f, data.frame(age=10, ld72=21, ld73=c(21, 47)), se.fit=TRUE)

## ----h=5,w=5,top=1-------------------------------------------------------
r <- resid(f)
par(mfrow=c(2,2))   # 2x2 matrix of plots
plot(fitted(f), r); abline(h=0)  # yhat vs. r
with(lead, plot(age,  r));    abline(h=0)
with(lead, plot(ld73, r));    abline(h=0)
qqnorm(r)           # linearity indicates normality
qqline(as.numeric(r))

## ----w=6,h=5-------------------------------------------------------------
ggplot(Predict(f))

## ------------------------------------------------------------------------
ggplot(Predict(f, age))   # plot age effect, using default range,
                          # 10th smallest to 10th largest age

## ------------------------------------------------------------------------
ggplot(Predict(f, age=3:15))  # plot age=3,4,...,15

## ------------------------------------------------------------------------
ggplot(Predict(f, age=seq(3,16,length=150)))   # plot age=3-16, 150 points

## ------------------------------------------------------------------------
ggplot(Predict(f, age, conf.type='individual'))

## ------------------------------------------------------------------------
p1 <- Predict(f, age, conf.int=0.99, conf.type='individual')
p2 <- Predict(f, age, conf.int=0.99, conf.type='mean')
p <- rbind(Individual=p1, Mean=p2)
ggplot(p)

## ------------------------------------------------------------------------
ggplot(Predict(f, ld73, age=3))

## ------------------------------------------------------------------------
ggplot(Predict(f, ld73, age=c(3,9)))  # add ,conf.int=FALSE to suppress conf. bands

## ----w=7,h=5-------------------------------------------------------------
bplot(Predict(f, ld72, ld73))

## ----w=6,h=4-------------------------------------------------------------
plot(nomogram(f))

## ------------------------------------------------------------------------
summary(f)         # inter-quartile-range effects
summary(f, age=5)  # adjust age to 5 when examining ld72,ld73
                   # (no effect since no interactions in model)
summary(f, ld73=c(20, 40))  # effect of changing ld73 from 20 to 40

## ------------------------------------------------------------------------
summary(f, age=5:6)    # starting age irrelevant since age is linear

## ----h=2.5,top=2---------------------------------------------------------
plot(summary(f))

## ------------------------------------------------------------------------
predict(f, data.frame(age=3, ld72=21, ld73=21))
# must specify all variables in the model

predict(f, data.frame(age=c(3, 10), ld72=21, ld73=c(21, 47)))
# predictions for (3,21,21) and (10,21,47)

newdat <- expand.grid(age=c(4, 8), ld72=c(21, 47), ld73=c(21, 47))
newdat
predict(f, newdat)     # 8 predictions

predict(f, newdat, conf.int=0.95)  # also get CLs for mean
predict(f, newdat, conf.int=0.95, conf.type='individual')  # CLs for indiv.

## ------------------------------------------------------------------------
# Model is b1 + b2*age + b3*ld72 + b4*ld73
b <- coef(f)
# For 3 year old with both lead exposures 21
b[1] + b[2]*3 + b[3]*21 + b[4]*21

## ------------------------------------------------------------------------
g <- Function(f)
g(age=c(3, 8), ld72=21, ld73=21)       # 2 predictions
g(age=3)              # 3 year old at median ld72, ld73

## ------------------------------------------------------------------------
an <- anova(f)
an                     # same as print(an)
print(an, 'names')     # print names of variables being tested
print(an, 'subscripts')# print subscripts in coef(f) (ignoring
                       # the intercept) being tested
print(an, 'dots')      # a dot in each position being tested
anova(f, ld72, ld73)   # combine effects into a 2 d.f. test

