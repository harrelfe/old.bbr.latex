f <- c(
  'rrms-1'='R rms: Model Formula Language',
  'rrms-2'='R rms: Operating on the Fit Object',
  'rrms-3'='R rms: datadist Function',
  'rrms-4'='R rms: lead Dataset',
  'rrms-5'='R rms: lead Residuals',
  'rrms-6'='R rms: lead Partial Effect Plots',
  'rrms-7'='R rms: lead Nomogram',
  'rrms-8'='R rms: lead Effect Estimation',
  'rrms-9'='R rms: lead Predicted Values',
  'rrms-10'='R rms: lead ANOVA',
  'reg-1'='Simple and Multiple Regression Models: Background',
  'reg-2'='Simple Linear Regression: Concepts',
  'reg-3'='Simple Linear Regression: Formal Model Statement and Parameter Estimation',
  'reg-4'='Simple Linear Regression: Inference About Parameters',
  'reg-5'='Simple Linear Regression: Interval Estimation',
  'reg-6'='Simple Linear Regression: Goodness of Fit',
  'reg-7'='Proper Transformations and Percentiling',
  'reg-8'='Multiple Linear Regression: Model and Parameter Estimation',
  'reg-9'='MLR: Interpretation of Parameters',
  'reg-9a'='MLR: Example: Estimation of Body Surface Area',
  'reg-10'='MLR: Degrees of Freedom and Testing Total Association',
  'reg-11'='MLR: Testing Partial Effects',
  'reg-12'='MLR: Assessing Goodness of Fit',
  'reg-13'='MLR: Binary Predictor',
  'reg-14'='MLR: Analysis of Covariance',
  'reg-15'='MLR: Correlation Coefficient Revisted',
  'reg-16'='MLR: Using Regression for ANOVA: lead Dataset',
  'reg-17'='MLR: More on ANOVA with Regression and ANCOVA',
  'reg-18'='MLR: Two-way ANOVA',
  'reg-19'='Heterogeneity of Treatment Effect',
  'reg-20'='Interaction Between Categorical and Continuous Predictor',
  'reg-21'='Internal vs. External Model Validation',
  'anc-1'='ANCOVA: Linear Models',
  'anc-2'='ANCOVA: Nonlinear Models',
  'anc-3'='ANCOVA: GUSTO-I Example',
  'anc-4'='ANCOVA: Nonlinear Models, General',
  'anc-5'='ANCOVA: Why are Adjusted Estimates Right?',
  'anc-6'='ANCOVA: Differential and Absolute Treatment Effects',
  'anc-7'='ANCOVA: Specifying Interactions',
  'anc-8'='ANCOVA: Strategy for Analyzing Differential Treatment Effect',
  'anc-9'='ANCOVA: Absolute vs. Relative Treatment Effects Revisited',
  'anc-10'='ANCOVA: Estimating Absolute Treatment Effects without Interaction',
  'anc-11'='ANCOVA: Estimating Absolute Treatment Effects with Interaction',
  'anc-12'='ANCOVA: Absolute Treatment Effects for GUSTO-I',
  'anc-13'='ANCOVA: Absolute Treatment Effect on Survival Probability',
  'anc-14'='ANCOVA: Cost-Effectiveness Ratios',
  'serial-1'='Serial Data: Introduction and Analysis Options',
  'serial-2'='Serial Data: Summary Measures',
  'serial-3'='Serial Data: Case Study - Data and Summary Measures',
  'serial-4'='Serial Data: Case Study - Nonparametric Tests',
  'serial-5'='Serial Data: Case Study - Generalized Least Squares, General Dose-Response',
  'serial-6'='Serial Data: GLS Linear in Logs',
  'serial-7'='Serial Data: GLS Final Estimates'

)

all <- character()
for(w in names(f)) {
  fn <- paste(w, 'm3u', sep='.')
  z <- paste('#EXTINF:-1,Frank Harrell - ', f[w], '\n',
             'http://hbiostat.org/audio/bbr/', w, '.mp3\n',
             sep='')
  cat(z, file=fn)
  all <- c(all, z)
}
cat(all, file='all.m3u', sep='')
