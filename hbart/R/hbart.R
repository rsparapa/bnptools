hbart = function(
x.train,
y.train,
x.test=matrix(0.0,0,0),
ntree=200,
ntreeh=40,
ndpost=1000, nskip=100,
k=2,
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
summarystats=FALSE
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
if(!is.null(xicuts)) # use xicuts
{
   xi=xicuts
}
else # default to equal numcut per dimension
{
   xi=vector("list",p)
   minx=apply(x,1,min)
   maxx=apply(x,1,max)
   for(i in 1:p)
   {
      xinc=(maxx[i]-minx[i])/(numcut+1)
      xi[[i]]=(1:numcut)*xinc+minx[i]
   }
}
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
   xi,
   summarystats,
   PACKAGE="hbart"
)

res$x.train=x.train
res$y.train=y.train+fmean
res$ntree=ntree
res$ntreeh=ntreeh
res$ndpost=ndpost
class(xi)="BARTcutinfo"
res$xicuts=xi

attr(res, 'class') <- 'hbart'

return(res)
}

