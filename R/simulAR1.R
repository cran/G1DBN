simulAR1<-function(p,n,edgeProp,minA,maxA,minB,maxB,minSig,maxSig){
#library library(MASS) 

# generated data 
data= matrix(0,p,n)
B=matrix(0,p,1)


# error
sigmaEps=runif(p,minSig,maxSig)

# B
B<-runif(p,minB,maxB)
neg<-sample(1:p,p/2,replace=FALSE)
B[neg]<-runif(length(neg),-maxB,-minB)



# edges
unp2<-seq(1,p^2,1)
aretes<-sample(unp2, round(edgeProp*p^2), replace=F)
aretesPlus<-sample(aretes,length(aretes)/2,replace=F)
aretesMoins<-aretes[which((aretes %in% aretesPlus)==F)]

A<-array(0,p^2)
nbPos<-round(length(aretes)/2)
A[which((unp2 %in% aretesPlus))]<-round(runif(nbPos,minA,maxA),1)
A[which((unp2 %in% aretesMoins))]<-round(runif((length(aretes)-nbPos),-maxA,-minA),1)
A<-matrix(A,p,p)


# AR(1) : random initialization for 1st time point 
data[,1]=B+rnorm(p,0,sigmaEps*10)

#  data generation
for (i in 2:n){
data[,i]=B+A %*% data[,(i-1)] + rnorm(p,0,sigmaEps)
}

list(data=t(data),A=A,B=B,sig=sigmaEps)
}
