## code to prepare `mixcure` dataset goes here

set.seed(82322)
n <- 2000
x1 <- runif(n); x2 <- runif(n); x3 <- runif(n); x4 <- runif(n); x5 <- runif(n); x6 <- runif(n); x7 <- runif(n); x8 <- runif(n); x9 <- runif(n); x10 <- runif(n)
bp <- sin(pi*x1*x2) - 4*(x6-0.5)^2 + 1*(x7>0.3) + 0.5*x8 - 0.8  #p mean
p <- pnorm(bp)  #non-cure probability
status <- rbinom(n, 1, prob = p)
mean <- -cos(4*pi*x2*x5) + 3*(2*x3-0.6)^2 - 2*(x1-0.2)*pmax(log(x4), -1) - 0.7
ltime <- rnorm(n, mean)
time <- exp(ltime)
ctime <- rexp(n, rate = 0.01)  #censoring_time~Exp(0.01)
mixcure <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6, x7 = x7, x8 = x8, x9 = x9, x10 = x10, r_time = time, p = p, ltime=mean, bp = bp)
mixcure$Ncure <- status
cuttime <- quantile(time, probs=0.95) 
mixcure$c_time <- pmin(ctime, cuttime)  #admin censoring
mixcure$obstime <- ifelse(mixcure$Ncure == 0, mixcure$c_time, pmin(mixcure$c_time, mixcure$r_time))
mixcure$event <- 1-(mixcure$obstime == mixcure$c_time)
usethis::use_data(mixcure, overwrite = TRUE)
