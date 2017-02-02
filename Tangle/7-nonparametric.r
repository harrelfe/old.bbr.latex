## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('nonpar')

## ----wsr-----------------------------------------------------------------
sr <- c(1.5, -1.5, 4, 3)
z <- sum(sr) / sqrt(sum(sr ^ 2))
pval <- 2 * (1 - pnorm(abs(z)))
c(z=z, pval=pval)

## ----wsrsleep------------------------------------------------------------
drug1 <- c(.7, -1.6, -.2, -1.2, -.1, 3.4, 3.7, .8, 0, 2)
drug2 <- c(1.9, .8, 1.1, .1, -.1, 4.4, 5.5, 1.6, 4.6, 3.4)
wilcox.test(drug2, drug1, paired=TRUE)
wilcox.test(drug2 - drug1)
wilcox.test(drug2 - drug1, correct=FALSE)
sr <- c(3, 8, 4.5, 4.5, 0, 2, 7, 1, 9, 6)
z <- sum(sr) / sqrt(sum(sr ^ 2))
c(z=z, pval=2 * (1 - pnorm(abs(z))))
d <- data.frame(Drug=c(rep('Drug 1', 10), rep('Drug 2', 10),
                  rep('Difference', 10)),
                extra=c(drug1, drug2, drug2 - drug1))


## ----signtest------------------------------------------------------------
2 * (1 / 2) ^ 9    # 2 * to make it two-tailed

## ----wsrcompare----------------------------------------------------------
# Assume we are already starting with signed ranks as x
wsr <- function(x, ...) wilcox.test(x, ...)$p.value
sim <- function(x) {
  z <- sum(x) / sqrt(sum(x ^ 2))
  2 * (1 - pnorm(abs(z))) }
all <- function(x) round(c(
  continuity=wsr(x, correct=TRUE, exact=FALSE),
  nocontinuity=wsr(x, correct=FALSE, exact=FALSE),
  exact=wsr(x, exact=TRUE),
  simple=sim(x)), 4)
all(1:4)
all(c(-1, 2 : 4))
all(c(-2, c(1, 3, 4)))
all(c(-1, -2, 3 : 5))
all(c(-5, -1, 2, 3, 4, 6))

## ----calpro,w=6,h=3,cap='Fecal calprotectin by endoscopy severity rating. Red dotted line is the detection limit.  Ordinal disease categories should not have been combined.',scap='Fecal calprotectin by severity'----
#Fecal Calprotectin: 2500 is above detection limit
calpro <- c(2500, 244, 2500, 726, 86, 2500, 61, 392, 2500, 114, 1226,
            2500, 168, 910, 627, 2500, 781, 57, 483, 30, 925, 1027,
            2500, 2500, 38, 18)

# Endoscopy score: 1 = No/Mild, 2=Mod/Severe Disease
# Would have been far better to code dose as 4 ordinal levels
endo <- c(2, 1, 2, 2, 2, 1, 1, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2,
          2, 2, 2, 2, 2, 1, 1)
endo <- factor(endo, 1 : 2,
               c("No or Mild Activity", "Moderate or Severe Activity"))
require(ggplot2)   # Fig. (*\ref{fig:nonpar-calpro}*)
ggplot(data.frame(endo, calpro), aes(y=calpro, x=endo)) +
  geom_boxplot(color='lightblue', alpha=.85, width=.4) +
  geom_dotplot(binaxis='y', stackdir='center', position='dodge') +
    xlab('') + ylab('Fecal Calprotectin') + coord_flip() +
      geom_hline(aes(yintercept=2500, col='red'), linetype='dotted')
wilcox.test(calpro ~ endo)

## ----calpror,w=6,h=3,cap='Ranks of calprotectin'-------------------------
ggplot(data.frame(endo, calpro), aes(y=rank(calpro), x=endo)) + #Fig (*\ref{fig:nonpar-calpror}*)
  geom_dotplot(binaxis='y', stackdir='center', position='dodge') +
    xlab('') + ylab('Rank of Fecal Calprotectin') + coord_flip()

## ----cendo---------------------------------------------------------------
require(Hmisc)
# Convert endo to a binary variable
somers2(calpro, endo=='Moderate or Severe Activity')

## ----ccode,eval=FALSE----------------------------------------------------
## mean.rank <- mean(rank(x)[y == 1])
## c.index <- (mean.rank - (n1 + 1)/2) / (n - n1)

## ----hlest---------------------------------------------------------------
female <- c(120, 118, 121, 119)
male   <- c(124, 120, 133)
differences <- outer(male, female, '-')
differences
median(differences)
wilcox.test(male, female, conf.int=TRUE)

## ----checkhl,cap='Wilcoxon $P$-value vs.\\ hypothesized male-female difference.  Horizontal line is $P=0.05$.  Vertical lines from left to right are the lower 0.95 confidence limit from \\Co{wilcox.test}, the median difference, the Hodges-Lehman estimator as computed by \\Co{wilcox.test}, and the upper 0.95 confidence limit from \\Co{wilcox.test}.',scap='Wilcoxon $P$-value vs.\\ hypothesized difference.'----
dif  <- seq(-3, 15, by=.1)
n    <- length(dif)
pval <- numeric(n)
for(i in 1 : n) pval[i] <- wilcox.test(male - dif[i], female)$p.value
ggplot(data.frame(dif, pval), aes(x=dif, y=pval)) +
  geom_step() +
  geom_hline(yintercept=.05, col='red', linetype='dotted') +
  geom_vline(xintercept=c(4.5, 4.791, -1, 15), col='red', linetype='dotted') +
  xlab('Difference') + ylab('P-value')

## ----diffmedboot---------------------------------------------------------
diffs <- numeric(1000)
set.seed(13)
for(i in 1 : 1000) diffs[i] <-
  median(sample(male, replace=TRUE)) - median(sample(female, replace=TRUE))
ggplot(data.frame(diffs), aes(x=diffs)) + xlab('Differences in Medians') +
  geom_histogram(bin_widith=.01, color='blue', fill='white')
quantile(diffs, c(0.025, 0.975))

