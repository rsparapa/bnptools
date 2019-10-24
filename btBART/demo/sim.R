## original file is/was ...
## ~btuyishimire/allfiles/RESEARCH/NEWBART/TEST/A1CDATA/a1cdata_full.R
library(btBART)
library(BART3)
##library(BayesTree)
library(lmeVarComp)
##source("~btuyishimire/allfiles/RESEARCH/Rfunctions/lpml.R")

set.seed(20L)
x <- matrix(runif(200L), 100L, 2L)
x1 <- x[, 1L]; x2 <- x[, 2L]

y <- 2*x1 + x2 + 10*x1*x2 + rnorm(100L)

ndpost <- 1000; burn <- 100

post1 <- wbart( x.train  = x,
                y.train  = y,   
                ndpost   = ndpost, 
                nskip    = burn)

sigs1 <- post1$sigma[(burn+1):(burn+ndpost)]
mus1 <- post1$yhat.train
out1 <- lpml(y,mus1,sigs1)

post2 <- twbarts( x.train1 = x1,
                  x.train2 = x2,
                  y.train  = y,   
                  ndpost   = ndpost, 
                  nskip    = burn)

sigs2 <- post2$sigma[(burn+1):(burn+ndpost)]
mus2 <- post2$yhat.train
out2 <- lpml(y,mus2,sigs2)

add_test1 <<- test.additivity(x, y)
pbf_ex1 <<- exp(out1$LPML - out2$LPML)
add_test1
pbf_ex1

set.seed(20L)
x <- matrix(runif(200L), 100L, 2L)
x1 <- x[, 1L]; x2 <- x[, 2L]

y <- 2*x1 + x2  + rnorm(100L)

post1 <- wbart( x.train  = x,
                y.train  = y,   
                ndpost   = ndpost, 
                nskip    = burn)

sigs1 <- post1$sigma[(burn+1):(burn+ndpost)]
mus1 <- post1$yhat.train
out1 <- lpml(y,mus1,sigs1)

post2 <- twbarts( x.train1 = x1,
                  x.train2 = x2,
                  y.train  = y,   
                  ndpost   = ndpost, 
                  nskip    = burn)

sigs2 <- post2$sigma[(burn+1):(burn+ndpost)]
mus2 <- post2$yhat.train
out2 <- lpml(y,mus2,sigs2)

add_test2 <<- test.additivity(x, y)
pbf_ex2 <<- exp(out1$LPML - out2$LPML)
add_test2
pbf_ex2

set.seed(20L)
x <- matrix(runif(200L), 100L, 2L)
x1 <- x[, 1L]; x2 <- x[, 2L]

y <- 2*x1 + x2 + 3*x1*x2 + rnorm(100L)

post1 <- wbart( x.train  = x,
                y.train  = y,   
                ndpost   = ndpost, 
                nskip    = burn)

sigs1 <- post1$sigma[(burn+1):(burn+ndpost)]
mus1 <- post1$yhat.train
out1 <- lpml(y,mus1,sigs1)

post2 <- twbarts( x.train1 = x1,
                  x.train2 = x2,
                  y.train  = y,   
                  ndpost   = ndpost, 
                  nskip    = burn)

sigs2 <- post2$sigma[(burn+1):(burn+ndpost)]
mus2 <- post2$yhat.train
out2 <- lpml(y,mus2,sigs2)

add_test3 <- test.additivity(x, y)
pbf_ex3 <- exp(out1$LPML - out2$LPML)
add_test3
pbf_ex3
