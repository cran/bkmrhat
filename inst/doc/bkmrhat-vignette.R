## ----metals data, echo=TRUE, results='markup', message=FALSE------------------
library("bkmr")
library("bkmrhat")
library("coda")
Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet") # for future package

set.seed(111)
dat <- bkmr::SimData(n = 50, M = 5, ind=1:3, Zgen="realistic")
y <- dat$y
Z <- dat$Z
X <- cbind(dat$X, rnorm(50))
head(cbind(y,Z,X))

## ----1 vs 1+ chains, cache=FALSE, results='markup'----------------------------

# enable parallel processing (up to 4 simultaneous processes here)
future::plan(strategy = future::multiprocess, workers=4, .skip=TRUE)

# single run of 4000 observations from bkmr package
set.seed(111)
system.time(kmfit <- suppressMessages(kmbayes(y = y, Z = Z, X = X, iter = 4000, verbose = FALSE, varsel = FALSE)))


# 4 runs of 1000 observations from bkmrhat package
set.seed(111)
system.time(kmfit5 <- suppressMessages(kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 1000, verbose = FALSE, varsel = FALSE)))


## ----diagnostics 1, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# Using rstan functions (set burnin/warmup to zero for comparability with coda numbers given later
#  posterior summaries should be performed after excluding warmup/burnin)
singlediag = kmbayes_diagnose(kmfit, warmup=0, digits_summary=2)


# Using rstan functions (multiple chains enable R-hat)
multidiag = kmbayes_diagnose(kmfit5, warmup=0, digits_summary=2)

# using coda functions, not using any burnin (for demonstration only)
kmfitcoda = as.mcmc(kmfit, iterstart = 1)
kmfit5coda = as.mcmc.list(kmfit5, iterstart = 1)

# single chain trace plot
traceplot(kmfitcoda)

## ----diagnostics 2, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# multiple chain trace plot
traceplot(kmfit5coda)

## ----diagnostics 2b, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# multiple cross-correlation plot (combines all samples)
crosscorr(kmfit5coda)
crosscorr.plot(kmfit5coda)

## ----diagnostics 2c, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# multiple chain trace plot
#autocorr(kmfit5coda) # lots of output
autocorr.plot(kmfit5coda)

## ----diagnostics 3, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# Gelman's r-hat using coda estimator (will differ from rstan implementation)
gelman.diag(kmfit5coda)
# effective sample size
effectiveSize(kmfitcoda)
effectiveSize(kmfit5coda)

## ----post summaries 1, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
# posterior kernel marginal densities using `mcmc` and `mcmc` objects
densplot(kmfitcoda)

## ----post summaries 2, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
# posterior kernel marginal densities using `mcmc` and `mcmc` objects
densplot(kmfit5coda)

## ----post summaries 3, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# posterior summaries using `mcmc` and `mcmc` objects
summary(kmfitcoda)
summary(kmfit5coda)

# highest posterior density intervals using `mcmc` and `mcmc` objects
HPDinterval(kmfitcoda)
HPDinterval(kmfit5coda)

# combine multiple chains into a single chain
fitkmccomb = kmbayes_combine(kmfit5)


# For example:
summary(fitkmccomb)


mean.difference <- suppressWarnings(OverallRiskSummaries(fit = fitkmccomb, y = y, Z = Z, X = X, 
                                      qs = seq(0.25, 0.75, by = 0.05), 
                                      q.fixed = 0.5, method = "exact"))
mean.difference

with(mean.difference, {
  plot(quantile, est, pch=19, ylim=c(min(est - 1.96*sd), max(est + 1.96*sd)), 
       axes=FALSE, ylab= "Mean difference", xlab = "Joint quantile")
  segments(x0=quantile, x1=quantile, y0 = est - 1.96*sd, y1 = est + 1.96*sd)
  abline(h=0)
  axis(1)
  axis(2)
  box(bty='l')
})

## ----varsel, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

set.seed(111)
system.time(kmfitbma.list <- suppressWarnings(kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 1000, verbose = FALSE, varsel = TRUE)))

bmadiag = kmbayes_diagnose(kmfitbma.list)

# posterior exclusion probability of each chain
lapply(kmfitbma.list, function(x) t(ExtractPIPs(x)))


kmfitbma.comb = kmbayes_combine(kmfitbma.list)
summary(kmfitbma.comb)
ExtractPIPs(kmfitbma.comb) # posterior inclusion probabilities

mean.difference2 <- suppressWarnings(OverallRiskSummaries(fit = kmfitbma.comb, y = y, Z = Z, X = X,                                       qs = seq(0.25, 0.75, by = 0.05), 
                                      q.fixed = 0.5, method = "exact"))
mean.difference2

with(mean.difference2, {
  plot(quantile, est, pch=19, ylim=c(min(est - 1.96*sd), max(est + 1.96*sd)), 
       axes=FALSE, ylab= "Mean difference", xlab = "Joint quantile")
  segments(x0=quantile, x1=quantile, y0 = est - 1.96*sd, y1 = est + 1.96*sd)
  abline(h=0)
  axis(1)
  axis(2)
  box(bty='l')
})


## ----post diagnostics, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

set.seed(111)
system.time(kmfitbma.list <- suppressWarnings(kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 1000, verbose = FALSE, varsel = TRUE)))

meandifference_par = OverallRiskSummaries_parallel(kmfitbma.list, y = y, Z = Z, X = X ,qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, method = "exact")

head(meandifference_par)
nchains = length(unique(meandifference_par$chain))

with(meandifference_par, {
  plot.new()
  plot.window(ylim=c(min(est - 1.96*sd), max(est + 1.96*sd)), 
              xlim=c(min(quantile), max(quantile)),
       ylab= "Mean difference", xlab = "Joint quantile")
  for(cch in seq_len(nchains)){
    width = diff(quantile)[1]
    jit = runif(1, -width/5, width/5)
   points(jit+quantile[chain==cch], est[chain==cch], pch=19, col=cch) 
   segments(x0=jit+quantile[chain==cch], x1=jit+quantile[chain==cch], y0 = est[chain==cch] - 1.96*sd[chain==cch], y1 = est[chain==cch] + 1.96*sd[chain==cch], col=cch)
  }
  abline(h=0)
  axis(1)
  axis(2)
  box(bty='l')
  legend("bottom", col=1:nchains, pch=19, lty=1, legend=paste("chain", 1:nchains), bty="n")
})

regfuns_par = PredictorResponseUnivar_parallel(kmfitbma.list, y = y, Z = Z, X = X ,qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, method = "exact")

head(regfuns_par)
nchains = length(unique(meandifference_par$chain))

# single variable
with(regfuns_par[regfuns_par$variable=="z1",], {
  plot.new()
  plot.window(ylim=c(min(est - 1.96*se), max(est + 1.96*se)), 
              xlim=c(min(z), max(z)),
       ylab= "Predicted Y", xlab = "Z")
  pc = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  pc2 = c("#0000001A", "#E69F001A", "#56B4E91A", "#009E731A", "#F0E4421A", "#0072B21A", "#D55E001A", "#CC79A71A", "#9999991A")
  for(cch in seq_len(nchains)){
   ribbonX = c(z[chain==cch], rev(z[chain==cch]))
   ribbonY = c(est[chain==cch] + 1.96*se[chain==cch], rev(est[chain==cch] - 1.96*se[chain==cch]))
   polygon(x=ribbonX, y = ribbonY, col=pc2[cch], border=NA)
   lines(z[chain==cch], est[chain==cch], pch=19, col=pc[cch]) 
  }
  axis(1)
  axis(2)
  box(bty='l')
  legend("bottom", col=1:nchains, pch=19, lty=1, legend=paste("chain", 1:nchains), bty="n")
})


## ----continue, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# install dev version of bkmr to allow true continued fits.
#install.packages("devtools")
#devtools::install_github("jenfb/bkmr")

set.seed(111)
# run 100 initial iterations for a model with only 2 exposures
Z2 = Z[,1:2]
kmfitbma.start <- suppressWarnings(kmbayes(y = y, Z = Z2, X = X, iter = 500, verbose = FALSE, varsel = FALSE))
kmbayes_diag(kmfitbma.start)

# run 2000 additional iterations
moreiterations = kmbayes_continue(kmfitbma.start, iter=2000)
kmbayes_diag(moreiterations)
TracePlot(moreiterations, par="beta")
TracePlot(moreiterations, par="r")



