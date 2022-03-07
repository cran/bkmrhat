
library(future)
library(bkmr)
library(bkmrhat)

set.seed(111)
dat <- bkmr::SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X
set.seed(111)
Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
future::plan(strategy = future::multiprocess)
fitkm.list <- kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 50,
                               verbose = FALSE, varsel = TRUE)

fit1 = fitkm.list[[1]]
fit2 = kmbayes_continue(fit1, iter=100)

fitkm.list2 <- kmbayes_parallel_continue(fitkm.list, iter=100)
d1 = kmbayes_diag(fitkm.list)
d2 = kmbayes_diag(fitkm.list2)

closeAllConnections()
