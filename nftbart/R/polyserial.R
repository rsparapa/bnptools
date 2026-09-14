# last modified 2021-12-07 by J. Fox

polyserial <- function(x, y, ML=FALSE, control=list(), std.err=FALSE, 
                       maxcor=.9999, bins=4, start, thresholds=FALSE){
	f <- function(pars){
		rho <- pars[1]
		if (abs(rho) > maxcor) rho <- sign(rho)*maxcor
		cts <- if (length(pars) == 1) c(-Inf, cuts, Inf)
			else c(-Inf, pars[-1], Inf)
        if (any(diff(cts) < 0)) return(Inf)
		tau <- (matrix(cts, n, s+1, byrow=TRUE) - matrix(rho*z, n, s+1))/
			sqrt(1 - rho^2)
		- sum(log(dnorm(z)*(pnorm(tau[cbind(indices, y+1)]) - pnorm(tau[cbind(indices, y)]))))
	}
	if (!is.numeric(x)) stop("x must be numeric")
	valid <- complete.cases(x, y)
	x <- x[valid]
	y <- y[valid]
	z <- scale(x)
	tab <- table(y)
	n <- sum(tab)
	s <- length(tab)
	if (s < 2) {
		warning("y has fewer than 2 levels")
		return(NA)
	}
	if (sum(0 != tab) < 2){
		warning("y has fewer than 2 levels with data")
		return(NA)
	}
	nzeros <- sum(zeros <- 0 == tab)
	if (nzeros > 0){
	  warning("the following ", if (nzeros == 1) " level" else " levels", " of y",
	          if (nzeros == 1) " has" else " have", " no cases: ", names(tab)[zeros])
	}
	indices <- 1:n
	cuts <- qnorm(cumsum(tab)/n)[-s]
	y <- as.numeric(as.factor(y))
	rho <- sqrt((n - 1)/n)*sd(y)*cor(x, y)/sum(dnorm(cuts))
	if (!missing(start) && (ML || std.err)) {
	  if (is.list(start)){
	    rho <- start$rho
	    cuts <- start$thresholds
	  } else {
	    rho <- start
	  }
	  if (!is.numeric(rho) || length(rho) != 1)
	    stop("start value for rho must be a number")
	  if (!is.numeric(cuts) || length(cuts) != s - 1) 
	    stop("start values for thresholds must be ", s - 1, "numbers")
	}
  if (abs(rho) > maxcor) {
    warning("initial correlation inadmissible, ", rho, ", set to ", sign(rho)*maxcor)
    rho <- sign(rho)*maxcor
  }
	if (ML) {
		result <- optim(c(rho, cuts), f, control=control, hessian=std.err)
		if (result$par[1] > 1){
			result$par[1] <- maxcor
			warning(paste("inadmissible correlation set to", maxcor))
		}
		else if (result$par[1] < -1){
			result$par[1] <- -maxcor
			warning(paste("inadmissible correlation set to -", maxcor, sep=""))
		}
		if (std.err){
			chisq <- chisq(y, z, result$par[1], result$par[-1], bins=bins)
			df <- s*bins - s  - 1
			result <- list(type="polyserial",
				rho=result$par[1],
				cuts=result$par[-1],
				var=solve(result$hessian),
				n=n,
				chisq=chisq,
				df=df,
				ML=TRUE)
			##class(result) <- "polycor"
			return(result)
		} else if (thresholds){
		  result <- list(type="polyserial",
		                 rho=result$par[1],
		                 cuts=result$par[-1],
		                 var=NA,
		                 n=n,
		                 chisq=NA,
		                 df=NA,
		                 ML=TRUE)
		  ##class(result) <- "polycor"
		  return(result)
		}
		else return(as.vector(result$par[1]))  
	}
	else if (std.err){
		result <- optim(rho, f, control=control, hessian=TRUE, method="BFGS")
		if (result$par > 1){
			result$par <- maxcor
			warning(paste("inadmissible correlation set to", maxcor))
		}
		else if (result$par < -1){
			result$par <- -maxcor
			warning(paste("inadmissible correlation set to -", maxcor, sep=""))
		}
		chisq <- chisq(y, z, rho, cuts, bins=bins)
		df <- s*bins - s  - 1
		result <- list(type="polyserial",
			rho=result$par,
			cuts=cuts,
			var=1/result$hessian,
			n=n,
			chisq=chisq,
			df=df,
			ML=FALSE)
		##class(result) <- "polycor"
		return(result)
	}
  else if (thresholds){
    result <- list(type="polyserial",
                   rho=rho,
                   cuts=cuts,
                   var=NA,
                   n=n,
                   chisq=NA,
                   df=NA,
                   ML=FALSE)
    ##class(result) <- "polycor"
    return(result)
  }
  else return(rho)
}
