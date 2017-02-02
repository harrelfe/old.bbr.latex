## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('obsvar', width=80)

## ----rex-----------------------------------------------------------------
d <- expand.grid(rep=1:2, observer=c('A','B','C'), subject=1:4)
d$y <- c(5,7, 8,5, 6,7,
         7,6, 8,6, 9,7,
         7,5, 4,6, 10,11,
         7,6, 5,6, 9,8)
d
# Function to compute mean absolute discrepancies
mad <- function(y, obs, subj) {
  nintra <- ninter <- sumintra <- suminter <- 0
  n <- length(y)
  for(i in 1 : (n - 1)) {
    for(j in (i + 1) : n) {
      if(subj[i] == subj[j]) {
        dif <- abs(y[i] - y[j])
        if(! is.na(dif)) {
          if(obs[i] == obs[j]) {
            nintra   <- nintra + 1
            sumintra <- sumintra + dif
          }
          else {
            ninter   <- ninter + 1
            suminter <- suminter + dif
          }
        }
      }
    }
  }
  c(nintra=nintra, intra=sumintra / nintra,
    ninter=ninter, inter=suminter / ninter)
}
  
# Compute statistics for first subject
with(subset(d, subject == 1), mad(y, observer, subject))
# Compute for all subjects
with(d, mad(y, observer, subject))

## ----robs----------------------------------------------------------------
source(paste('http://biostat.mc.vanderbilt.edu/wiki/pub/Main',
             'AnalysisOfObserverVariability/observerVariability_3_3.R',
             sep='/'))
require(Hmisc)
with(d, {
  intra <- intraVar(subject, observer, y)
  print(intra)
  summary(intra)
  set.seed(2)
  b=bootStrap(intra, by = 'subject', times=1000)
  # Get 0.95 CL for mean absolute intra-observer difference
  print(quantile(b, c(0.025, 0.975)))
  inter <- interVar(subject, observer, y)
  print(inter)
  summary(inter)
  b <- bootStrap(inter, by = 'subject', times=1000)
  # Get 0.95 CL for mean absolute inter-observer difference
  print(quantile(b, c(0.025, 0.975)))
})

