## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('multgroup', width=80)

## ----walking,h=2.25,w=5,cap='Age in months when infants first began walking by treatment group with mean lines',scap='Age of first walking'----
w <- rbind(
  data.frame(trt='Active',      months=c(9,9.5,9.75,10,13,9.5)),
  data.frame(trt='Passive',     months=c(11,10,10,11.75,10.5,15)),
  data.frame(trt='No Exercise', months=c(11.5,12,9,11.5,13.25,13)),
  data.frame(trt='8-Week Control', months=c(13.25,11.5,12,13.5,11.5,12.35)) )
aggregate(months ~ trt, w, function(x) c(Mean=mean(x), Variance=var(x)))
require(ggplot2)
require(data.table)
w <- data.table(w)
stats <- w[, j=list(months = mean(months), var=var(months)), by = trt]

ggplot(w, aes(x=trt, y=months)) +    # Fig. (*\ref{fig:multgroup-walking}*)
  geom_dotplot(binaxis='y', stackdir='center', position='dodge') +
  geom_errorbar(aes(ymin=..y.., ymax=..y..), width=.7, size=1.3,
                data=stats) +
  xlab('') + ylab('Months Until First Walking') + coord_flip()

## ----walk-anova----------------------------------------------------------
require(rms)
f <- ols(months ~ trt, data=w)
anova(f)

## ----fex-----------------------------------------------------------------
pf(58, 4, 1044)
1 - pf(58, 4, 1044)

## ----kw------------------------------------------------------------------
anova(ols(rank(months) ~ trt, data=w))
spearman2(months ~ trt, data=w)
kruskal.test(months ~ trt, data=w)

