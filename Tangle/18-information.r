## ----echo=FALSE----------------------------------------------------------
require(Hmisc)
knitrSet('info', width=80)

## ----pneuwho,w=5,h=3.75,cap='Estimated risk of pneumonia with respect to two predictors in WHO ARI study from \\citet{har98dev}.  Tick marks show data density of respiratory rate stratified by cough.  Any cutpoint for the rate \\textbf{must} depend on cough to be consistent with optimum decision making, which must be risk-based.',scap='Risk of pneumonia with two predictors'----
require(rms)
getHdata(ari)
r <- ari[ari$age >= 42, Cs(age, rr, pneu, coh, s2)]
abn.xray <- r$s2==0
r$coh <- factor(r$coh, 0:1, c('no cough','cough'))
f <- lrm(abn.xray ~ rcs(rr,4)*coh, data=r)
anova(f)
dd <- datadist(r); options(datadist='dd')
p <- Predict(f, rr, coh, fun=plogis, conf.int=FALSE)
ggplot(p, rdata=r,     # Fig. (*\ref{fig:info-pneuwho}*)
       ylab='Probability of Pneumonia',
       xlab='Adjusted Respiratory Rate/min.',
       ylim=c(0,.7), legend.label='')

## ----thresholds,echo=FALSE,w=8,h=3,mfrow=c(1,2),cap='Two kinds of thresholds.  The pattern on the left represents a discontinuity in the first derivative (slope) of the function relating a marker to outcome.  On the right there is a lowest-order discontinuity.',scap='Two kinds of thresholds'----
plot(0:1, 0:1, type='n', axes=FALSE, xlab='Marker', ylab='Outcome')
axis(1, labels=FALSE)
axis(2, labels=FALSE)
lines(c(.1, .5), c(.25, .25))
lines(c(.5, .9), c(.25, .9))

plot(0:1, 0:1, type='n', axes=FALSE, xlab='Marker', ylab='Outcome')
axis(1, labels=FALSE)
axis(2, labels=FALSE)
lines(c(.1, .5), c(.25, .25))
lines(c(.5, .9), c(.75, .75))
lines(c(.5, .5), c(.25, .75), col=gray(.9))   # Fig. (*\ref{fig:info-thresholds}*)

## ----coxsim--------------------------------------------------------------
set.seed(1)
n <- 1000
age <- rnorm(n, mean=50, sd=12)
describe(age)
cens <- 15 * runif(n)
h  <- 0.02 * exp(0.04 * (age - 50))
dt <- -log(runif(n))/h
e  <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
S  <- Surv(dt, e)
coef(cph(S ~ age))   # close to true value of 0.04 used in simulation
exp(coef(cph(S ~ age >= 50)))   # >=50 : < 50 hazard ratio estimate
exp(coef(cph(S ~ age >= 50, subset=age < 60)))
exp(coef(cph(S ~ age >= 50, subset=age < 55)))
exp(coef(cph(S ~ age >= 50, subset=age > 40)))
exp(coef(cph(S ~ age >= 50, subset=age > 40 & age < 60)))

## ----psa,w=4.5,h=3.5,cap='Relationship between post-op PSA level and 2-year recurrence risk.  Horizontal lines represent the only prognoses provided by the new staging system. Data are courtesy of M Kattan from JNCI 98:715; 2006.  Modification of AJCC staging by Roach \\emph{et al.} 2006.',scap='Continuous PSA vs.\ risk'----
load('~/doc/Talks/infoAllergy/kattan.rda')
attach(kattan)
t   <- t.stg
gs  <- bx.glsn
psa <- preop.psa
t12 <- t.stg %in% Cs(T1C,T2A,T2B,T2C)

s <- score.binary(t12 & gs<=6 & psa<10,
                  t12 & gs<=6 & psa >=10 & psa < 20,
                  t12 & gs==7 & psa < 20,
                  (t12 & gs<=6 & psa>=20) |
                  (t12 & gs>=8 & psa<20),
                  t12 & gs>=7 & psa>=20,
                  t.stg=='T3')
levels(s) <- c('none','I', 'IIA', 'IIB', 'IIIA', 'IIIB', 'IIIC')
u <- is.na(psa + gs) | is.na(t.stg)
s[s=='none'] <- NA
s <- s[drop=TRUE]
s3 <- s
levels(s3) <- c('I','II','II','III','III','III')
table(s3)
units(time.event) <- 'month'
dd <- datadist(data.frame(psa, gs)); options(datadist='dd')
S <- Surv(time.event, event=='YES')
label(psa) <- 'PSA'; label(gs) <- 'Gleason Score'
f <- cph(S ~ rcs(sqrt(psa), 4), surv=TRUE, x=TRUE, y=TRUE)
p <- Predict(f, psa, time=24, fun=function(x) 1 - x)
h <- cph(S ~ s3, surv=TRUE)
z <- 1 - survest(h, times=24)$surv
ggplot(p, rdata=data.frame(psa), ylab='2-year Disease Recurrence Risk') +
  geom_hline(yintercept=unique(z), col='red', size=0.2)   # Fig. (*\ref{fig:info-psa}*)

## ----spectrum,w=5.75,h=6,cap='Prognostic spectrum from various models with model $\\chi^2$ - d.f., and generalized $c$-index.  The mostly vertical segmented line connects different prognostic estimates for the same man.',scap='Prognostic spectrum from various models'----
f <- cph(S ~ rcs(sqrt(psa),4) + pol(gs,2), surv=TRUE)    
g <- function(form, lab) {
  f <- cph(form, surv=TRUE, subset=!u)
  cat(lab,'\n'); print(coef(f))
  s <- f$stats
  cat('N:', s['Obs'],'\tL.R.:', round(s['Model L.R.'],1),
      '\td.f.:',s['d.f.'],'\n\n')
  prob24 <- 1 - survest(f, times=24)$surv
  prn(sum(!is.na(prob24)))
  p2 <<- c(p2, prob24[2])  # save est. prognosis for one subject
  p1936 <<- c(p1936, prob24[1936])
  C <- rcorr.cens(1-prob24, S[!u,])['C Index']
  data.frame(model=lab, chisq=s['Model L.R.'], d.f.=s['d.f.'],
             C=C, prognosis=prob24)
}
p2 <- p1936 <- NULL
w <-          g(S ~ t.stg, 'Old Stage')
w <- rbind(w, g(S ~ s3, 'New Stage'))
w <- rbind(w, g(S ~ s, 'New Stage, 6 Levels'))
w <- rbind(w, g(S ~ pol(gs,2),        'Gleason'))
w <- rbind(w, g(S ~ rcs(sqrt(psa),4), 'PSA'))
w <- rbind(w, g(S ~ rcs(sqrt(psa),4) + pol(gs,2), 'PSA+Gleason'))
w <- rbind(w, g(S ~ rcs(sqrt(psa),4) + pol(gs,2) + t.stg,
                'PSA+Gleason+Old Stage'))

w$z <- paste(w$model, '\n',
             'X2-d.f.=',round(w$chisq-w$d.f.),
             '  C=', sprintf("%.2f", w$C), sep='')
w$z <- with(w, factor(z, unique(z)))
require(lattice)
stripplot(z ~ prognosis, data=w, lwd=1.5,    # Fig. (*\ref{fig:info-spectrum}*)
          panel=function(x, y, ...) {
            llines(p2, 1:7, col=gray(.6))
            ## llines(p1936, 1:7, col=gray(.8), lwd=2)
            ## panel.stripplot(x, y, ..., jitter.data=TRUE, cex=.5)
            for(iy in unique(unclass(y))) {
              s <- unclass(y)==iy
              histSpike(x[s], y=rep(iy,sum(s)), add=TRUE, grid=TRUE)
            }
            panel.abline(v=0, col=gray(.7))
          },
          xlab='Predicted 2-year\nDisease Recurrence Probability')

## ----cutpointExists,w=5,h=4----------------------------------------------
require(rms)
set.seed(4)
n <- 900
X <- rnorm(n, 100, 20)
dd <- datadist(X); options(datadist='dd')

p <- ifelse(X < 100, .2, .8)
y <- ifelse(runif(n) <= p, 1, 0)

f <- lrm(y ~ rcs(X, c(90,95,100,105,110)))
hs <- function(yval, side)
  histSpikeg(yhat ~ X, data=subset(data.frame(X, y), y == yval),
             side = side, ylim = c(0, 1),
             frac = function(f) .03 * f / max(f))
ggplot(Predict(f, fun=plogis), ylab='Probability of Response') +
  hs(0, 1) + hs(1, 3) + geom_vline(xintercept=100, col=gray(.7))

## ----devlogistread,echo=FALSE--------------------------------------------
d <- read.table(textConnection(
"0.29017806553  -26.5670446561
0.251361744241  -11.2913339428
0.314187070796  -8.22259324454
0.416584193946  8.4673672441
0.236463771182  18.8611965218
0.307040274389  18.4124441893
0.369895819955  22.6370621076
0.346460975923  26.2542778789
0.291689016145  31.2268163515
0.228924127615  30.4698300936
0.401595563849  35.1522660482
0.378221157841  41.0812362598
0.307644654634  41.5299885923
0.229286755762  44.3403567354
0.292202739353  50.8767290941
0.253114446954  55.7495448261
0.331623440888  58.7185627838
0.394267453368  54.8520401611
0.402471915205  68.6727054326
0.449704231417  75.3088005319
0.39481139559  75.6578301238
0.35572310319  80.5306458558
0.300709391314  76.2561665672
0.237974721796  76.6550575294
0.167307561552  73.6361782015
0.395143804725  88.3724795455
0.340311406922  91.0332635778
0.293199966759  89.0206773592
0.246179183633  90.475722801
0.207121110247  96.5044157532
0.301283552547  98.2178337501
0.277939365552  105.302681182
0.23882085414  109.019619694
0.278271774687  118.017330604
0.380215612653  117.36913279
0.434927134407  110.084839877
0.39614103213  126.516427811
0.396443222253  138.075200012
0.302310998965  137.517659235
0.223832224044  135.704518498
0.216232142452  145.001397629
0.294710917374  146.814538367
0.295103764533  161.840942229
0.280175572461  190.837595473
0.272726585932  205.913860705
0.729864694372  -8.55349142913
0.730317979557  8.78466687316
0.651688109574  1.1921400349
0.644088027983  10.4890191664
0.581323139453  9.7320329085
0.652262270808  23.1538072178
0.573843933911  23.6524209206
0.636971450588  38.2799338204
0.621559754319  48.7825515423
0.676694342245  57.6805397116
0.621801506418  58.0295693035
0.56690867059  58.3785988955
0.724410162654  82.8106703332
0.740335582131  91.9579653539
0.61477558606  89.2881156179
0.630670786525  97.2795334185
0.615319528281  110.093905581
0.560608006527  117.378198494
0.639117000461  120.347216451
0.678416825946  123.56554126
0.624037713327  143.564483595
0.624370122463  156.279133017
0.618945809756  248.799171999
0.697303708628  245.988803856
1.38326017814  383.82376272
1.27238662204  342.910241979
1.60225736022  360.46597717
1.60101838071  313.075011143
1.17650169604  275.311822433
1.23136431286  273.806915621
1.29421985842  278.031533539
1.27823400092  266.572484078
1.24662491406  257.524911798
1.30929914556  254.814266396
1.19155076416  250.938678069
1.24644359999  250.589648477
1.28547145437  243.405078305
1.22255547077  236.868705946
1.29310175497  235.264076394
1.29282978386  224.861181412
1.36337606805  223.25655186
1.19073485083  219.729993125
1.24562768666  219.380963533
1.29252759374  213.302409211
1.30012767533  204.005530079
1.28414181782  192.546480618
1.23696993964  188.222139959
1.28386984671  182.143585637
1.2445095832  176.613506388
1.2521096648  167.316627256
1.29916066693  167.017459034
1.23591227421  147.766437254
1.28290283832  145.155514592
1.21989619769  135.151510573
1.29041226288  132.3910038
1.29798212545  121.938247448
1.25853120491  112.940536539
1.30555198803  111.485491097
1.19561522132  106.40416418
1.25029652406  97.9639940469
1.29749862126  103.444211926
1.29725686916  94.1971941647
1.30491738878  87.2120694735
1.40686122674  86.5638716599
1.45367047678  77.0176856769
1.38306375456  76.3105607893
1.31927141961  36.2537490462
1.31124827185  29.3683470956
1.23286015397  31.0228380185
1.40355224489  259.995316053
1.40303852169  240.34540331
1.39341376627  172.198508692
1.46417158355  178.68501968
1.67017459034  258.300029464
1.62306315018  256.287443245
1.61491912637  244.778532414
1.63042147967  237.743546352
1.66996305726  250.208888922
1.71719537347  256.844984022
1.71689318335  245.28621182
1.7244026079  232.521701028
1.56726374398  221.960156232
1.61434496513  222.816865231
1.62194504673  213.519986099
1.61389167995  205.478706928
1.59796626047  196.331411908
1.59775472739  188.240271367
1.67644503539  198.144552645
1.6767170065  208.547447627
1.7315191853  204.730786374
1.62891052906  179.949685345
1.69161497956  178.394917162
1.70714755188  172.515808321
1.76197994969  169.855024289
1.62085716228  171.908406174
1.57359462706  164.116433854
1.63620842053  159.094034011
1.62013190599  144.16735289
1.63560404028  135.976489608
1.69031556204  128.692196696
1.62749023548  125.623455997
1.57262761867  127.128362809
1.58022770026  117.831483678
1.62730892141  118.688192676
1.73664130788  100.65197519
1.62652322709  88.6353849524
1.62628147499  79.3883671912
1.64123988607  51.547591167
2.38189276783  481.520318508
2.38122794956  456.091019665
2.26224058867  404.824465312
2.32497525818  404.42557435
2.37971699895  398.297158657
2.33254512076  393.972817998
2.28534302356  388.492600119
2.2384733355  395.727031662
2.17582932302  399.593554285
2.62260231024  388.660315638
2.61355171606  342.475088202
2.62893319332  330.81659326
2.62078916951  319.307682428
2.72243081735  307.100712413
2.27580892519  323.813337161
2.2834694448  316.82821247
2.27538585901  307.631056079
2.28295572159  297.178299727
2.22064411825  313.759471772
2.56492932529  282.668640975
2.62769421381  283.425627233
2.61958040901  273.072593622
2.29007229899  269.387385073
2.61858318161  234.928645357
2.71962044921  199.604130939
2.28792674911  187.320102442
2.29465047935  144.502783927
2.29350215688  100.579449561
2.62988509221  67.2267256945
2.72664636956  168.345584625
2.73421623214  157.892828273
3.27140450414  605.344232324
3.71527646619  483.446780542
3.28198115845  409.901259377"), col.names=c('x','y'))

## ----devlogist-----------------------------------------------------------
# Points from published graph were defined in code not printed
g <- trunc(d$x)
g <- factor(g, 0:3, c('very well', 'fair to good', 'poor', 'very poor'))
remiss <- 1 * (g == 'very well')
CDAI <- d$y
label(CDAI) <- "Crohn's Disease Activity Index"
label(remiss) <- 'Remission'
dd <- datadist(CDAI,remiss); options(datadist='dd')
f <- lrm(remiss ~ rcs(CDAI,4))
ggplot(Predict(f, fun=plogis), ylab='Probability of Remission')

