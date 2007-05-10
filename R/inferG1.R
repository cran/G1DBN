inferG1<-function(data,ls=T,tukey=F,huber=F,predictor=NULL,target=NULL){

#library library(MASS) 
data<-t(data)
p=length(data[,1])  # nb de genes
n=length(data[1,]) # nb de repetitions pour chaque gene

tar=1:p
pred=1:p

if(length(predictor)>0){pred=predictor}
if(length(target)>0){tar=target}

r=length(tar)
d=length(pred)

pmax=NULL
pmaxH=NULL
pmaxT=NULL
# score matrices :
if(ls==T)pmax<-matrix(0,r,d)
if(huber==TRUE)pmaxH<-matrix(0,r,d)
if(tukey==TRUE)pmaxT<-matrix(0,r,d)


for (i in 1:r){
print(i)

temp<-matrix(0,n-1,3)

temp[,1]<-data[tar[i],2:n]                   # without time point 1

for (j in  c(1:(d-1))){

for (k in c(1:d)[-c(1:j)]){        #  for all k > j

                         
temp[,2:3]<-t(data[c(pred[j],pred[k]),1:(n-1)])   # without time n
# vector of ros wothout NA
nomiss<-which((seq(1,n-1,by=1) %in% c(which(is.na(temp))-(trunc((which(is.na(temp))-1)/(n-1)) )*(n-1)))==F)
 
y=temp[nomiss,1]                               
lag1=temp[nomiss,2:3]   # regressor matrice 



#################### Least Square ######################
if(ls==T){
lm.3<-lm(y~lag1)
phat<-abs(summary(lm.3)$coef[,"Pr(>|t|)"])  

# aij(k) : aij given k 
if(!is.na(phat[2]) && phat[2]>pmax[i,j]){
pmax[i,j]<-phat[2]}

# aik(j) : aik given j 
if(!is.na(phat[3]) && phat[3]>pmax[i,k]){
pmax[i,k]<-phat[3]}
}


################## Tukey bisquare ######################
if(tukey==TRUE){
bisq.3<-rlm(y~lag1,method='MM')
prob<-pt(abs(summary(bisq.3)$coef[,"t value"]),n-4,lower.tail = F, log.p = FALSE)*2

# coefficient aij(k) : aij given k 
if(  prob[2]>pmaxT[i,j]){
pmaxT[i,j]<-prob[2]}

# coefficient aik(j) : aik given j (
if(!is.na(prob[3]) && prob[3]>pmaxT[i,k]){
pmaxT[i,k]<-prob[3]}
}

########################## Huber ######################
if(huber==TRUE){
hub.3<-rlm(y~lag1)

prob<-pt(abs(summary(hub.3)$coef[,"t value"]),n-4,lower.tail = F, log.p = FALSE)*2

# coefficient aij(k) : aij given k 
if(!is.na(prob[2]) && prob[2]>pmaxH[i,j]){
pmaxH[i,j]<-prob[2]}

# coefficient aik(j) : aik given j 
if(!is.na(prob[3]) && prob[3]>pmaxH[i,k]){
pmaxH[i,k]<-prob[3]}
}

} # end k
} # end j
} # end i

list(ls=pmax,tukey=pmaxT,huber=pmaxH)
}
