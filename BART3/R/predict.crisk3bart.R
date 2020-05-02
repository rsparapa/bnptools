
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2018 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

predict.crisk3bart <- function(object, newdata,
                               newdata2=newdata, newdata3=newdata,
                               mc.cores=1, openmp=(mc.cores.openmp()>0),
                               mult.impute=mc.cores, seed=99, ...) {

    p <- length(object$treedraws$cutpoints)
    if(p!=ncol(newdata))
        stop(paste0('The number of columns in newdata must be equal to ', p))

    p2 <- length(object$treedraws2$cutpoints)
    if(p2!=ncol(newdata2))
        stop(paste0('The number of columns in newdata2 must be equal to ',p2))

    p3 <- length(object$treedraws3$cutpoints)
    if(p3!=ncol(newdata3))
        stop(paste0('The number of columns in newdata3 must be equal to ',p3))

    n=nrow(newdata)
    n2=nrow(newdata2)
    n3=nrow(newdata3)

    if(n!=n2)
        stop(paste0('The number of rows in newdata2 must be equal to ', n))

    if(n!=n3)
        stop(paste0('The number of rows in newdata3 must be equal to ', n))

    if(.Platform$OS.type == "unix") mc.cores.detected <- detectCores()
    else mc.cores.detected <- NA

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected)
        mc.cores <- mc.cores.detected

    if(length(object$binaryOffset)==0) object$binaryOffset=object$offset
    if(length(object$binaryOffset2)==0) object$binaryOffset2=object$offset2
    if(length(object$binaryOffset3)==0) object$binaryOffset3=object$offset3

    miss=apply(is.na(newdata), 2, sum)
    miss2=apply(is.na(newdata2), 2, sum)
    miss3=apply(is.na(newdata3), 2, sum)
    names(miss)=names(object$treedraws$cutpoints)
    names(miss2)=names(object$treedraws2$cutpoints)
    names(miss3)=names(object$treedraws3$cutpoints)
    miss.=((sum(miss)+sum(miss2)+sum(miss3))>0)

    if(!miss.) mult.impute=1
    else {
        set.seed(seed)
        newdata.=newdata
        newdata2.=newdata2
        newdata3.=newdata3
    }

    pred<-list()
    tx.test<-list()
    tx.test2<-list()
    tx.test3<-list()

    ## hot deck missing imputation
    for(k in 1:mult.impute) {
        if(miss.) {
            if(k==1) {
                warning(paste0('missing elements of x imputed with hot decking ',
                               mult.impute, ' times'))

                check=(p==p2 && p==p3)
                check2=(p==p2)
                check3=(p==p3)
            } else {
                newdata=newdata.
                newdata2=newdata2.
                newdata3=newdata3.
            }

            for(i in 1:n)
                for(j in 1:p) {
                    if(k==1) {
                        if(check) check=((is.na(newdata[i, j]) && is.na(newdata2[i, j]) && is.na(newdata3[i, j])) ||
                                         (!is.na(newdata[i, j]) && !is.na(newdata2[i, j]) && !is.na(newdata3[i, j]) &&
                                          newdata[i, j]==newdata2[i, j] && newdata[i, j]==newdata3[i, j]))
                        if(check2) check2=((is.na(newdata[i, j]) && is.na(newdata2[i, j])) ||
                                           (!is.na(newdata[i, j]) && !is.na(newdata2[i, j]) &&
                                            newdata[i, j]==newdata2[i, j]))
                        if(check3) check3=((is.na(newdata[i, j]) && is.na(newdata3[i, j])) ||
                                           (!is.na(newdata[i, j]) && !is.na(newdata3[i, j]) &&
                                            newdata[i, j]==newdata3[i, j]))
                    }
                    if(is.na(newdata[i, j]) && miss[j]==n) {
                        h=length(object$treedraws$cutpoints[[j]])
                        if(h==1) newdata[i, j]=object$treedraws$cutpoints[[j]][1]
                        else newdata[i, j]=object$treedraws$cutpoints[[j]][sample.int(h, 1)]
                    }
                    else while(is.na(newdata[i, j])) {
                        h=sample.int(n, 1)
                        newdata[i, j]=newdata[h, j]
                    }
                }

            if(check) {
                newdata2=newdata
                newdata3=newdata
            }
            else {
                if(check2) newdata2=newdata
                else if(check3) newdata3=newdata

                for(i in 1:n) {
                    if(!check2)
                        for(j in 1:p2) {
                            if(is.na(newdata2[i, j]) && miss2[j]==n) {
                                h=length(object$treedraws2$cutpoints[[j]])
                                if(h==1) newdata2[i, j]=object$treedraws2$cutpoints[[j]][1]
                                else newdata2[i, j]=object$treedraws2$cutpoints[[j]][sample.int(h, 1)]
                            }
                            else while(is.na(newdata2[i, j])) {
                                h=sample.int(n, 1)
                                newdata2[i, j]=newdata2[h, j]
                            }
                        }
                    if(!check3)
                        for(j in 1:p3) {
                            if(is.na(newdata3[i, j]) && miss3[j]==n) {
                                h=length(object$treedraws3$cutpoints[[j]])
                                if(h==1) newdata3[i, j]=object$treedraws3$cutpoints[[j]][1]
                                else newdata3[i, j]=object$treedraws3$cutpoints[[j]][sample.int(h, 1)]
                            }
                            else while(is.na(newdata3[i, j])) {
                                h=sample.int(n, 1)
                                newdata3[i, j]=newdata3[h, j]
                            }
                        }
                }
            }
        }

        pred[[k]]=mc.crisk3.pwbart(newdata, newdata2, newdata3,
                                   object$treedraws, object$treedraws2,
                                   object$treedraws3,
                                   object$binaryOffset, object$binaryOffset2,
                                   object$binaryOffset3,
                                   mc.cores=mc.cores, type=object$type, ...)

        if(mult.impute>1) {
            if(k==1) {
                pred[[1]]$yhat.test =pred[[1]]$yhat.test/mult.impute
                pred[[1]]$yhat.test2=pred[[1]]$yhat.test2/mult.impute
                pred[[1]]$yhat.test3=pred[[1]]$yhat.test3/mult.impute
                pred[[1]]$prob.test =pred[[1]]$prob.test/mult.impute
                pred[[1]]$prob.test2=pred[[1]]$prob.test2/mult.impute
                pred[[1]]$prob.test3=pred[[1]]$prob.test3/mult.impute
                pred[[1]]$surv.test =pred[[1]]$surv.test/mult.impute
                pred[[1]]$cif.test  =pred[[1]]$cif.test/mult.impute
                pred[[1]]$cif.test2 =pred[[1]]$cif.test2/mult.impute
                pred[[1]]$cif.test3 =pred[[1]]$cif.test3/mult.impute
            } else {
                pred[[1]]$yhat.test =pred[[1]]$yhat.test +pred[[k]]$yhat.test /mult.impute
                pred[[k]]$yhat.test =NULL
                pred[[1]]$yhat.test2=pred[[1]]$yhat.test2+pred[[k]]$yhat.test2/mult.impute
                pred[[k]]$yhat.test2=NULL
                pred[[1]]$yhat.test3=pred[[1]]$yhat.test3+pred[[k]]$yhat.test3/mult.impute
                pred[[k]]$yhat.test3=NULL
                pred[[1]]$prob.test =pred[[1]]$prob.test +pred[[k]]$prob.test /mult.impute
                pred[[k]]$prob.test =NULL
                pred[[1]]$prob.test2=pred[[1]]$prob.test2+pred[[k]]$prob.test2/mult.impute
                pred[[k]]$prob.test2=NULL
                pred[[1]]$prob.test3=pred[[1]]$prob.test3+pred[[k]]$prob.test3/mult.impute
                pred[[k]]$prob.test3=NULL
                pred[[1]]$surv.test =pred[[1]]$surv.test +pred[[k]]$surv.test /mult.impute
                pred[[k]]$surv.test =NULL
                pred[[1]]$cif.test  =pred[[1]]$cif.test  +pred[[k]]$cif.test  /mult.impute
                pred[[k]]$cif.test  =NULL
                pred[[1]]$cif.test2 =pred[[1]]$cif.test2 +pred[[k]]$cif.test2 /mult.impute
                pred[[k]]$cif.test2 =NULL
                pred[[1]]$cif.test3 =pred[[1]]$cif.test3 +pred[[k]]$cif.test3 /mult.impute
                pred[[k]]$cif.test3 =NULL
            }

            tx.test[[k]]=pred[[k]]$tx.test
            pred[[k]]$tx.test=NULL
            tx.test2[[k]]=pred[[k]]$tx.test2
            pred[[k]]$tx.test2=NULL
            tx.test3[[k]]=pred[[k]]$tx.test3
            pred[[k]]$tx.test3=NULL
        }
    }

    pred[[1]]$mult.impute=mult.impute
    pred[[1]]$miss=miss
    pred[[1]]$miss2=miss2
    pred[[1]]$miss3=miss3

    if(mult.impute>1) {
        pred[[1]]$cif.test.mean  = apply(pred[[1]]$cif.test, 2, mean)
        pred[[1]]$cif.test2.mean = apply(pred[[1]]$cif.test2, 2, mean)
        pred[[1]]$cif.test3.mean = apply(pred[[1]]$cif.test3, 2, mean)
        pred[[1]]$surv.test.mean = apply(pred[[1]]$surv.test, 2, mean)
        pred[[1]]$tx.test=tx.test
        pred[[1]]$tx.test2=tx.test2
        pred[[1]]$tx.test3=tx.test3
    }

    return(pred[[1]])
}

