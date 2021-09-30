hbart = function(
x.train,
y.train,
x.test=matrix(0.0,0,0),
ntree=200,
ntreeh=40,
ndpost=1000, nskip=100,
k=5, ##2,
power=2.0, base=.95,
tc=1,
sigmav=rep(1,length(y.train)),
fmean=mean(y.train),
overallsd = sd(y.train),
overallnu=10,
chv = cor(x.train,method="spearman"),
pbd=.7,
pb=.5,
stepwpert=.1,
probchv=.1,
minnumbot=5,
printevery=100,
numcut=100,
xicuts=NULL,
nadapt=1000,
adaptevery=100,
summarystats=TRUE
)
{
#require(Rcpp)

#--------------------------------------------------
nd = ndpost
burn = nskip
m = ntree
mh = ntreeh
#--------------------------------------------------
#data
n = length(y.train)
p = ncol(x.train)
np = nrow(x.test)
x = t(x.train)
xp = t(x.test)
y.train = y.train-fmean
#--------------------------------------------------
#cutpoints
if(is.null(xicuts)) xicuts=xicuts(x.train, numcut=numcut)
## if(!is.null(xicuts)) # use xicuts
## {
##    xi=xicuts
## }
## else # default to equal numcut per dimension
## {
##    xi=vector("list",p)
##    minx=apply(x,1,min)
##    maxx=apply(x,1,max)
##    for(i in 1:p)
##    {
##       xinc=(maxx[i]-minx[i])/(numcut+1)
##       xi[[i]]=(1:numcut)*xinc+minx[i]
##    }
## }
#--------------------------------------------------
rgy = range(y.train)
#tau =  (rgy[2]-rgy[1])/(sqrt(m)*k) this is not consistent with BART
tau =  (rgy[2]-rgy[1])/(2*sqrt(m)*k)

#--------------------------------------------------
overalllambda = overallsd^2
#--------------------------------------------------
powerh=power
baseh=base
if(length(power)>1) {
   powerh=power[2]
   power=power[1]
}
if(length(base)>1) {
   baseh=base[2]
   base=base[1]
}
#--------------------------------------------------
pbdh=pbd
pbh=pb
if(length(pbd)>1) {
   pbdh=pbdh[2]
   pbd=pbd[1]
}
if(length(pb)>1) {
   pbh=pb[2]
   pb=pb[1]
}
#--------------------------------------------------
stepwperth=stepwpert
if(length(stepwpert)>1) {
   stepwperth=stepwpert[2]
   stepwpert=stepwpert[1]
}
#--------------------------------------------------
probchvh=probchv
if(length(probchv)>1) {
   probchvh=probchv[2]
   probchv=probchv[1]
}
#--------------------------------------------------
minnumboth=minnumbot
if(length(minnumbot)>1) {
   minnumboth=minnumbot[2]
   minnumbot=minnumbot[1]
}

#--------------------------------------------------
#call
res=.Call("cpsambrt",
   x,
   y.train,
   xp,
   m,
   mh,
   nd,
   burn,
   nadapt,
   adaptevery,
   tau,
   overalllambda,
   overallnu,
   base,
   power,
   baseh,
   powerh,
   tc,
   sigmav,
   chv,
   pbd,
   pb,
   pbdh,
   pbh,
   stepwpert,
   stepwperth,
   probchv,
   probchvh,
   minnumbot,
   minnumboth,
   printevery,
#   numcut,
##   xi,
   xicuts,
   summarystats,
   PACKAGE="hbart"
)

res$x.train=x.train
res$y.train=y.train+fmean
if(summarystats) {
    ## names(res$mu.varcount)=names(xicuts)
    ## res$mu.varprob=res$mu.varcount/sum(res$mu.varcount)
    ## names(res$sd.varcount)=names(xicuts)
    ## res$sd.varprob=res$sd.varcount/sum(res$sd.varcount)
    dimnames(res$mu.varcount)[[2]]=names(xicuts)
    res$mu.varcount[2:ndpost, ]=res$mu.varcount[2:ndpost, ]-res$mu.varcount[1:(ndpost-1), ]
    res$mu.varcount.mean=apply(res$mu.varcount, 2, mean)
    res$mu.varprob=res$mu.varcount.mean/sum(res$mu.varcount.mean)
    dimnames(res$sd.varcount)[[2]]=names(xicuts)
    res$sd.varcount[2:ndpost, ]=res$sd.varcount[2:ndpost, ]-res$sd.varcount[1:(ndpost-1), ]
    res$sd.varcount.mean=apply(res$sd.varcount, 2, mean)
    res$sd.varprob=res$sd.varcount.mean/sum(res$sd.varcount.mean)
}
res$ntree=ntree
res$ntreeh=ntreeh
res$ndpost=ndpost
## class(xi)="BARTcutinfo"
## res$xicuts=xi
res$xicuts=xicuts
## res$treedraws=list()
## res$treedraws$cutpoints=xicuts
## res$treedraws$f.trees=res$f.trees
## res$f.trees=NULL
## res$treedraws$s.trees=res$s.trees
## res$s.trees=NULL

attr(res, 'class') <- 'hbart'

    res$pred=predict(res, res$x.train, soffset=0)
    res$soffset=0.5*log(mean(res$pred$s.test.mean^2)/
                   mean(res$pred$s.train.mean^2)) ## for stability
    ##if(!is.finite(res$b0)) res$b0=NA
    
return(res)
}

