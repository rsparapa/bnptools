\name{normwrapper}
\alias{normwrapper}
\docType{package}
\title{
A wrapper for the DPdensity function using the LIO prior parameter settings
}

\description{ A wrapper for the DPdensity function from the R package
DPpackage. This wrapper automatically uses the Low Information Omnibus
prior in generating a posterior density sample for a Dirichlet mixture
of Gaussians distributions.  }

\usage{
normwrapper(y, y50=NULL, y95=NULL, ngrid=1000, grid=NULL,
mcmc=list(nburn = 1000, nsave = 1000, nskip = 10, ndisplay = 100))
}

\arguments{
\item{y}{A vector or matrix giving the data from which the density estimate is to be computed.}
\item{y50}{A vector of specified medians of coordinates the data generating distribution.}
\item{y95}{A vector of specified 95th percentiles of coordinates of the data generating distribution.}
\item{mcmc}{A list giving the MCMC parameters. The list must include the following integers:
nburn giving the number of burn-in scans (the default value is 1000), nskip giving the thinning interval (the default value is 10), nsave giving the total number of scans (the default value is 1000), and ndisplay giving the number of saved scans to be displayed on screen (the function reports on the
screen when every ndisplay iterations have been carried out. The default value is 100).}
\item{ngrid}{Number of grid points where the density estimate is evaluated. This is only used if the dimension of y is lower or equal to 2. The default value is 1000.}
\item{grid}{Matrix of dimension ngrid*nvar of grid points where the density estimate is evaluated. This is only used if the dimension of y is lower or equal to 2. The default value is NULL, under which the grid is chosen according to the range of the data.}
}
\details{
This generic function fits a Dirichlet process mixture of Gaussians model for density estimation (Escobar
and West, 1995). Using estimates for the medians \bold{m} and 95th percentiles \bold{c} of coordinates of the data's distribution, the original data is rescaled as \eqn{\bold{z}_i = 2} \eqn{Diag(\bold{c})^{-1} ( \bold{y_i-m})}. Then, we fit the Gaussian DPM model 
\deqn{
\begin{array}{rl}
\mathbf{z}_i | \boldsymbol{\mu}_i, \mathbf{T}_i &\sim No(\boldsymbol{\mu}_i,\mathbf{T}_i), \\
(\boldsymbol{\mu}_i, \mathbf{T}_i) | G &\sim G, \\
G|G_0 &\sim DP(\nu, G_0), \\
G_0|\lambda,\boldsymbol{\Psi} &= NoWi(\mathbf{m}_\mu, \lambda, k_T, \boldsymbol{\Psi}), \\
\lambda &\sim Ga(a_\lambda, b_\lambda), \\
\boldsymbol{\Psi} &\sim Wi(k_\psi, \mathbf{W}_\psi), \\
\nu &\sim Ga(a_\nu,b_\nu), \\
\end{array}}
where \eqn{No(\bold{m},\bold{U})} denotes a normal distribution with mean \bold{m} and precision matrix \bold{U}, \eqn{Ga(a,b)} denotes a Gamma distribution with shape parameter a and rate parameter b, \eqn{Wi(k,\bold{W})} denotes a Wishart distribution with degrees of freedom k and rate matrix \bold{W} (expectation \eqn{k \bold{W}^{-1}}), and \eqn{NoWi(\bold{m}, \lambda, k, \bold{\Psi})} denotes a Normal-Wishart distribution with location \bold{m}, precision factor \eqn{\lambda}, degrees of freedom k, and rate matrix \bold{W}.
The LIO prior specifies the hyperparameters as \eqn{\bold{m}_\mu = \bold{0}, k_T = p+2, a_\lambda = 3/2, b_\lambda = v^2/2, k_\psi = p}, and \eqn{\bold{W}_\psi = p \bold{I}}, where \eqn{v^2 = 100 p (n-1)/ [n} \eqn{\chi_{p,0.99}^{2} ]} and \eqn{p = dim} \eqn{\bold{y}_i}. Details of this choice of prior are provided in our forthcoming paper (Shi et al., 2017).
}
\value{
An object of class "DPdensity", whose elements include:

\item{y}{Data set used for estimation}
\item{mcmc}{List of MCMC parameters}
\item{save.state$randsave}{Matrix containing MCMC samples of cluster parameters}
\item{grid1}{First coordinates of points at which the density is estimated}
\item{grid2}{Second coordinates of points at which the density is estimated}
\item{dens}{Density estimates at points in grid}
\item{fun1}{Marginal density estimates at points in grid1}
\item{fun2}{Marginal density estimates at points in grid2}

Since Jara's DPdensity function generates the MCMC samples, the DPdensity object contains many other elements that we have left unmodified. Interpretation of these elements might differ from that in the original function due to the rescaling process performed in the wrapper.
}
\source{
Escobar,M.D. and West,M. (1995) Bayesian density estimation and inference using mixtures, \emph{Journal of the American Statistical Association} \bold{90} Num 430, 577--588
	
Jara,A. and Hanson,T.E. and Quintana,F.A. and Muller,P. and Rosner,G.L.(2011) DPpackage: Bayesian semi-and nonparametric modeling in R, \emph{Journal of Statistical Software} \bold{40} Num 5, 1 

Shi,Y. and Martens,M. and Banerjee,A. and Laud,P. (2017) Low Information Omnibus Priors for Dirichlet Process Mixture Models(Manuscript)
}

\examples{
\dontrun{
library(DPWeibull)
# Scalar data from gamma(2,1)

n <-  200
y <- rgamma(n,2,1)
# Specify percentiles
fit <- normwrapper(y=y,y50=1,y95=4)

plot(fit$dens~fit$grid1,xlim=c(0,8),type="l")
curve(dgamma(x,2),xlim=c(0,8),lty=2,add=TRUE)
rug(y)

############################################################################

# Bivariate t / normal mixture
library(mvtnorm)

df1 <- Inf
mu1 <- c(2,0)
T1 <- 3*solve(matrix(c(1,1,1,4),nrow=2))
df2 <- 5
mu2 <- c(0,0)
T2 <- (df2-2)/df2*matrix(c(1,0,0,1),c(2,2))

n <- 400
ratio <- 0.5
n1 <- rbinom(1,n,ratio)
n2 <- n-n1
Y1 <- rmvt(n1,df=df1,sigma=solve(T1),delta=mu1)
Y2 <- rmvt(n2,df=df2,sigma=solve(T2),delta=mu2)
Y <- rbind(Y1,Y2)

# MCMC settings
nburn = 1000
nsave = 1000
nskip = 0
ndisplay = 1000
mcmc = list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

ngrid <- 1000
mesh1 <- 7/100;
mesh2 <- 6/100
grid <- cbind(seq(-3,4,mesh1),seq(-3,3,mesh2))

# Use sample percentiles
fit <- normwrapper(y=Y,mcmc=mcmc,ngrid=ngrid,grid=grid)

image(x=fit$grid1,y=fit$grid2,z=fit$dens,col=terrain.colors(12))
contour(x=fit$grid1,y=fit$grid2,z=fit$dens,add=TRUE)

y50 <- c(0,0)
y95 <- c(3,3)

# Specify percentiles
fit2 <- normwrapper(y=Y,y50=y50,y95=y95,mcmc=mcmc,ngrid=ngrid,grid=grid)

image(x=fit2$grid1,y=fit2$grid2,z=fit2$dens,col=terrain.colors(12))
contour(x=fit2$grid1,y=fit2$grid2,z=fit2$dens,add=TRUE)

############################################################################

# Air quality data, a real data example
set.seed(13)
data("airquality")
Y <- cbind(airquality$Ozone,airquality$Solar.R)
Y <- Y[!is.na(rowSums(Y)),]
n <- nrow(Y)
p <- 2

# MCMC settings
nburn <- 1000
nsave <- 1000
nskip <- 0
ndisplay <- 1000
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
ngrid <- 10000
grid <- NULL

# Use sample percentiles
fit <- normwrapper(y=Y,mcmc=mcmc,ngrid=ngrid,grid=grid)

# Scatter plot
plot(Y[,2]~Y[,1],cex = 0.5,xlab = "Ozone (parts per billion)",
     ylab = "Solar radiation (Langleys)")

# Contour plot of bivariate density
image(x=fit$grid1,y=fit$grid2,z=fit$dens,col=terrain.colors(12),
      xlab = "Ozone (parts per billion)",
      ylab = "Solar radiation (Langleys)",main = "Density estimate")
contour(x=fit$grid1,y=fit$grid2,z=fit$dens,add=TRUE)

# Marginal density plots
plot(fit$fun1~fit$grid1,type="l",xlab = "Ozone (parts per billion)",
     ylab = "Density",main = "Marginal density estimate")
plot(fit$fun2~fit$grid2,type="l",xlab = "Solar radiation (Langleys)",
     ylab = "Density",main = "Marginal density estimate")
}
}
