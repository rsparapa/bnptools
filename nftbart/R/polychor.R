# last modified 2022-01-10 by J. Fox

polychor <- function (x, y, ML=FALSE, control=list(), std.err=FALSE, 
                      maxcor=0.9999, start, thresholds=FALSE){
        f <- function(pars) {
            if (length(pars) == 1){
                rho <- pars
                if (abs(rho) > maxcor) rho <- sign(rho)*maxcor
                row.cuts <- rc
                col.cuts <- cc
            }
            else {
                rho <- pars[1]
                if (abs(rho) > maxcor) rho <- sign(rho)*maxcor
                row.cuts <- pars[2:r]
                col.cuts <- pars[(r+1):(r+c-1)]
                if (any(diff(row.cuts) < 0) || any(diff(col.cuts) < 0)) return(Inf)
            }
            P <- binBvn(rho, row.cuts, col.cuts)
            - sum(tab * log(P))
        }
        tab <- if (missing(y)) x else table(x, y)
        zerorows <- apply(tab, 1, function(x) all(x == 0))
        zerocols <- apply(tab, 2, function(x) all(x == 0))
        zr <- sum(zerorows)
        if (0 < zr) warning(paste(zr, " row", suffix <- if(zr == 1) "" else "s",
                                  " with zero marginal", suffix," removed", sep=""))
        zc <- sum(zerocols)
        if (0 < zc) warning(paste(zc, " column", suffix <- if(zc == 1) "" else "s",
                                  " with zero marginal", suffix, " removed", sep=""))
        tab <- tab[!zerorows, ,drop=FALSE]  
        tab <- tab[, !zerocols, drop=FALSE] 
        r <- nrow(tab)
        c <- ncol(tab)
        if (r < 2) {
            warning("the table has fewer than 2 rows")
            return(NA)
        }
        if (c < 2) {
            warning("the table has fewer than 2 columns")
            return(NA)
        }
        n <- sum(tab)
        rc <- qnorm(cumsum(rowSums(tab))/n)[-r]
        cc <- qnorm(cumsum(colSums(tab))/n)[-c]
        
        if (!missing(start) && (ML || std.err)) {
            if (is.list(start)){
                rho <- start$rho
                rc <- start$row.thresholds
                cc <- start$col.thresholds
            } else {
                rho <- start
            }
            if (!is.numeric(rho) || length(rho) != 1)
                stop("start value for rho must be a number")
            if (!is.numeric(rc) || length(rc) != r - 1) 
                stop("start values for row thresholds must be ", r - 1, "numbers")
            if (!is.numeric(cc) || length(cc) != c - 1) 
                stop("start values for column thresholds must be ", c - 1, "numbers")
        }
        
        if (ML) {
            result <- optim(
                c(if (missing(start)) optimise(f, interval=c(-1, 1))$minimum else rho, rc, cc), 
                f, control=control, hessian=std.err
            )
            if (result$par[1] > 1){
                result$par[1] <- maxcor
                warning(paste("inadmissible correlation set to", maxcor))
            }
            else if (result$par[1] < -1){
                result$par[1] <- -maxcor
                warning(paste("inadmissible correlation set to -", maxcor, sep=""))
            }
            if (std.err) {
                chisq <- 2*(result$value + sum(tab * log((tab + 1e-6)/n)))
                df <- length(tab) - r - c
                result <- list(type="polychoric",
                               rho=result$par[1],
                               row.cuts=result$par[2:r],
                               col.cuts=result$par[(r+1):(r+c-1)],
                               var=solve(result$hessian),
                               n=n,
                               chisq=chisq,
                               df=df,
                               ML=TRUE)
                ##class(result) <- "polycor"
                return(result)
            } else if (thresholds){
                    result <- list(type="polychoric",
                                   rho=result$par[1],
                                   row.cuts=result$par[2:r],
                                   col.cuts=result$par[(r+1):(r+c-1)],
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
            result <- optim(if (missing(start)) 0 else rho, 
                            f, control=control, hessian=TRUE, method="BFGS")
            if (result$par > 1){
                result$par <- maxcor
                warning(paste("inadmissible correlation set to", maxcor))
            }
            else if (result$par < -1){
                result$par <- -maxcor
                warning(paste("inadmissible correlation set to -", maxcor, sep=""))
            }
            chisq <- 2*(result$value + sum(tab *log((tab + 1e-6)/n)))
            df <- length(tab) - r - c 
            result <- list(type="polychoric",
                           rho=result$par,
                           row.cuts=rc,
                           col.cuts=cc,
                           var=1/result$hessian,
                           n=n,
                           chisq=chisq,
                           df=df,
                           ML=FALSE)
            ##class(result) <- "polycor"
            return(result)
        } else {
                rho <- optimise(f, interval=c(-maxcor, maxcor))$minimum
                if (thresholds){
                        result <- list(type="polychoric",
                                       rho=rho,
                                       row.cuts=rc,
                                       col.cuts=cc,
                                       var=NA,
                                       n=n,
                                       chisq=NA,
                                       df=NA,
                                       ML=FALSE)
                        ##class(result) <- "polycor"
                        return(result) 
                } else {
                        return(rho)
                }
        }
}
