GfromG1<-function(S1,data,method='ls',alpha1,alpha2=1,predictor=NULL,target=NULL){

data<-t(data)
r=length(S1[,1]) # nb of target genes 
d=length(S1[1,]) # nb of predictor genes
n=length(data[1,]) # nb of time points 

tar=1:r
pred=1:d

if(length(predictor)>0){pred=predictor}
if(length(target)>0){tar=target}


S2=matrix(1,r,d)
nbpar<-apply(S1<alpha1,1,sum) # nb de parents de chaque gene cible dans G(1)

if (max(nbpar)>(n-1)){
print(paste("Warning: Threshold alpha1 is two high, some nodes have more than ",n-1," parents in the inferred DAG G(1)"))}else{

# for genes having parents in G1
for (i in which(nbpar>=1)){
selec1<-which(S1[i,]<alpha1)


y=as.matrix(data[tar[i],2:n])
lag1=matrix(t(data[pred[selec1],1:(n-1)]),n-1,length(selec1))

temp<-matrix(0,n-1,length(selec1)+1)
temp[,1]<-data[tar[i],2:n] # without time poin

                        
temp[,2:(length(selec1)+1)]<-t(as.matrix(data[pred[selec1],1:(n-1)]))   # without time n

# vector of ros wothout NA
nomiss<-which((seq(1,n-1,by=1) %in% c(which(is.na(temp))-(trunc((which(is.na(temp))-1)/(n-1)) )*(n-1)))==F)
 
y=temp[nomiss,1]                               
lag1=temp[nomiss,2:(length(selec1)+1)]   # predictors matrix

################### Least Square ########################
if(method=='ls'){ls<-lm(y~lag1)
if(length(summary(ls)$coeff[,"Pr(>|t|)"])!=(length(selec1)+1)){print(paste("pb regression gene",i))
}else{S2[i,selec1]<-summary(ls)$coeff[2:(length(selec1)+1),"Pr(>|t|)"]}
}
################## Tukey bisquare ######################
if(method=='tukey'){
tuk<-rlm(y~lag1,method='MM')
if(length(summary(tuk)$coef[,"t value"])!=(length(selec1)+1)){print(paste("pb regression gene",i))
}else{S2[i,selec1]<-pt(abs(summary(tuk)$coef[2:(length(selec1)+1),"t value"]),n-length(selec1)-1,lower.tail = F, log.p = FALSE)*2}
}
###################### Huber ###########################
if(method=='huber'){
hub<-rlm(y~lag1)
if(length(summary(hub)$coef[,"t value"])!=(length(selec1)+1)){print(paste("pb regression gene",i))
}else{S2[i,selec1]<-pt(abs(summary(hub)$coef[2:(length(selec1)+1),"t value"]),n-length(selec1)-1,lower.tail = F, log.p = FALSE)*2}
}


}
}

nbparG<-apply(S2<alpha2,1,sum)
########## sortie :  ##############
list(S2=S2,maxInG1=max(nbpar),edgesNbG1=sum(nbpar),edgesNbG=sum(nbparG))
}
