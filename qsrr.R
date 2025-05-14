library(Qtools)
library(modi)
library(numDeriv)
library(compiler)
library(np)

expit=function(x) exp(x)/(1+exp(x))

qsrr=function(y,x,tau1=0.8,tau2=0.2,n.sample=NULL,n.sample2=NULL,se=FALSE,weights=NULL) {

n=length(y)
if(is.null(weights)) {weights=rep(1,n)}
if(is.null(n.sample)) {n.sample=n}
if(is.null(n.sample2)) {n.sample2=n}
if(n<n.sample) {n.sample=n}
if(n<n.sample2) {n.sample2=n}
p=ncol(x)

# step 1

# estimate conditional CDF and conditional QF

if(!is.data.frame(x)) {x=as.data.frame(x)}
yo=sort(unique(y))

ys=yo[seq(1,n,length=n.sample+2)][-c(1,n.sample+2)]
cfs.ns=sapply(ys,function(th) coef(glm(ifelse(y<=th,1,0)~.-1,data=x,weights=weights,family=binomial)))
cfs=t(apply(cfs.ns,1,function(x) approx(ys,x,yo)$y))

x=as.matrix(x)

pseudoy=sapply(1:n, function(i) {
ex=sort(as.vector(expit(x[i,]%*%cfs)))
w1=n
w2=1
if(any(ex>=tau1)) {w1=which(ex>=tau1)}
if(any(ex<=tau2)) {w2=which(ex<=tau2)}
sum(yo[w1]*diff(c(tau1,ex[w1])))/sum(yo[w2]*diff(c(0,ex[w2])))})

if(any(pseudoy<=0)) {pseudoy[pseudoy<=0]=min(pseudoy[pseudoy>0])}
betahat=coef(lm(I(log(pseudoy))~.-1,data=as.data.frame(x),weights=weights))

sds=rep(NA,p)
if(se) {

sX=solve(t(sweep(x,1,w,"*"))%*%x)
se.1=sX%*%t(x)%*%diag(w*(pseudoy-exp(as.numeric(x%*%betahat)))^2)%*%x%*%sX

# se.y=sX%*%t(x)%*%diag(log(y)-as.numeric(x%*%betahat))%*%x%*%sX

pb=function(vc,ys) {
cfs=matrix(as.vector(vc),nrow=p)
cfs=t(apply(cfs,1,function(x) approx(ys,x,yo)$y))
pseudoy=sapply(1:n, function(i) {
ex=sort(as.vector(expit(xmat[i,]%*%cfs)))
w1=n
w2=1
if(any(ex>=tau1)) {w1=which(ex>=tau1)}
if(any(ex<=tau2)) {w2=which(ex<=tau2)}
sum(yo[w1]*diff(c(tau1,ex[w1])))/sum(yo[w2]*diff(c(0,ex[w2])))})
coef(lm(I(log(pseudoy))~.-1,data=x,weights=weights))} 

pbj=function(vc,h,ys) {
cfs.ns2[,h]=vc 
cfs=t(apply(cfs.ns2,1,function(x) approx(ys,x,yo)$y))
pseudoy=sapply(1:n, function(i) {
ex=sort(as.vector(expit(xmat[i,]%*%cfs)))
w1=n
w2=1
if(any(ex>=tau1)) {w1=which(ex>=tau1)}
if(any(ex<=tau2)) {w2=which(ex<=tau2)}
sum(yo[w1]*diff(c(tau1,ex[w1])))/sum(yo[w2]*diff(c(0,ex[w2])))})
coef(lm(I(log(pseudoy))~.-1,data=x,weights=weights))}

ys2=yo[seq(1,n,length=n.sample2+2)][-c(1,n.sample2+2)]
x=as.data.frame(x)

# zero jacobian for flat ecdf, can be pruned 

de=diff(ecdf(y)(ys2))
if(any(de==0)) {
wNA=which(de==0)+1
ys2=ys2[-wNA]}

cfs.ns2=sapply(ys2,function(th) coef(glm(ifelse(y<=th,1,0)~.-1,data=x,family=binomial,weights=weights)))

xmat=as.matrix(x)
ex=expit(xmat%*%cfs.ns2)

# zero jacobian for non-crossing conditional quantiles, can be pruned
de=apply(ex,2,min)>tau2 & apply(ex,2,max)<tau1
if(any(de))  {
wNA=which(de)
ys2=ys2[-wNA]
ex=ex[,-wNA]
cfs.ns2=cfs.ns2[,-wNA]}

vars=matrix(0,prod(dim(cfs.ns2)),prod(dim(cfs.ns2)))
for(k in 1:length(ys2)) {
vars[1:p+(k-1)*p,1:p+(k-1)*p]=vcov(glm(ifelse(y<=ys2[k],1,0)~.-1,data=x,family=binomial,weights=weights))}

gr=jacobian(pb,as.vector(cfs.ns2),method="simple",ys=ys2)
se.2=gr%*%vars%*%t(gr)

sds=sqrt(diag(se.2+se.1))}

return(list(betahat=betahat,se=sds))}

qsrr=cmpfun(qsrr) 
