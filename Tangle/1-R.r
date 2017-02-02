## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('r')

## ----getRsex,eval=FALSE--------------------------------------------------
## require(Hmisc)      # do this once per session (or library(Hmisc))
## options(url.method='libcurl')    # sometimes needed if using Windows
## getRs()             # list available scripts
## getRs(browse='browser')  # open scripts contents in your web browser
## scripts <- getRs()  # store directory of scripts in an object that can easily
##                     # be viewed on demand in RStudio (right upper pane)
## getRs('introda.r')  # download introda.r and open in script editor
## getRs(cats=TRUE)    # list available major and minor categories
## categories <- getRs(cats=TRUE)  # store results in a list for later viewing
## getRs(cats='reg')   # list all scripts in a major category containing 'reg'
## getRs('importREDCap.r', put='source')   # source() to define a function

## ----importred,eval=FALSE------------------------------------------------
## require(Hmisc)
## getRs('importREDCap.r', put='source')  # source() code to define function
## mydata <- importREDCap()  # by default operates on last downloaded export
## Save(mydata)   # Hmisc function to create mydata.rda in compressed format

## ----impcsv--------------------------------------------------------------
# What is between data <- .. and ') is exactly like an external .csv file
data <- textConnection('
Age in Years,sex,race,visit date,m/s
23,m,w,10/21/2014,1.1
14,f,b,10/22/2014,1.3
,f,w,10/15/2014,1.7
')
require(Hmisc)
d <- csv.get(data, lowernames=TRUE, datevars='visit.date',
             dateformat='%m/%d/%Y')
close(data)
# lowernames=TRUE: convert variable names to lower case
# Omit dateformat if dates are in YYYY-MM-DD format
contents(d)
d

## ----impcsv2-------------------------------------------------------------
require(Hmisc)
d <- csv.get('test.csv')
names(d)   # show names after modification by csv.get
contents(d)  # show labels created by csv.get
d <- upData(d,
            state=factor(state, 1:2, c('Alabama','Alaska')),
            rename=c(sys.bp='sbp'),
            labels=c(age = 'Age',
                     sbp = 'Systolic Blood Pressure'),
            drop='junk',   # for > 1: drop=c('junk1','junk2',...)
   units=c(sbp='mmHg'))
contents(d)
describe(d)
dim(d); nrow(d); ncol(d); length(d)  # length is no. of variables

## ----eval=FALSE----------------------------------------------------------
## d <- data.frame(age=c(10,20,30), sex=c('male','female','male'),
##                 sbp=c(120,125,NA))

## ----subset,eval=FALSE---------------------------------------------------
## young.males <- subset(d, sex == 'male' & age > 26)
## # If you want to exclude rows that are missing on sex or age:
## young.males <- subset(d, sex == 'male' & age > 26 & ! is.na(sex) &
##                         ! is.na(age))
## # f <- lrm(y ~ sex + age, data=subset(d, sex == 'male' & ...))
## # f <- lrm(y ~ sex + age, data=d, subset=sex == 'male' & age > 26 ...)

