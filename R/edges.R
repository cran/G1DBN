edges<-function(score,predNames=NULL,targetNames=NULL,validMat=NULL,roc=FALSE,threshold=1,nb=NULL,prec=3){
r=length(score[,1]) # nb of target genes 
d=length(score[1,]) # nb of predictor genes

target=1:r
pred=1:d
rocx=NULL
rocy=NULL
edgesNames=NULL
if(length(predNames)>0)pred=predNames
if(length(targetNames)>0)target=targetNames

scoreOrd<-order(score,decreasing=F)

scoreTri<-sort(score,decreasing=F)
nbar=min(nb,sum(score<threshold),length(scoreTri))

#nbar=length(scoreTri)

#m<-array(c(t(validMat)),r*d)

m<-c(validMat)
y<-array(0,nbar)

edgesRef<-matrix(0,nbar,3+(length(validMat)>0)*1)
if((length(predNames)+length(targetNames))>0){edgesNames<-edgesRef}

lag<-ceiling(scoreOrd[1:nbar]/r)  # posCol -> j (ceiling : renvoie la plus petite valeur supérieure)
y<-scoreOrd[1:nbar]-(lag-1)*r  # posLine -> i

for (i in 1:nbar){
edgesRef[i,]<-cbind(lag[i],y[i],round(score[y[i],lag[i]],prec),m[scoreOrd[i]]==1)
if((length(predNames)+length(targetNames))>0){edgesNames[i,]<-cbind(pred[lag[i]],target[y[i]],round(score[y[i],lag[i]],prec),(m[scoreOrd[i]]==1)*1)}
}


if(roc==T){
if(length(validMat)==0){print("Please indicate a validation matrix")}else{
edgesPos=which((abs(validMat)>0))

testEdge<-scoreOrd%in% edgesPos

rocx<-array(0,r*d+1) # for x coord roc curve
rocy<-array(0,r*d+1) # for y coord

for (i in 1:(r*d)){
if(testEdge[i]==T){
rocy[i+1]<-rocy[i]+1
rocx[i+1]<-rocx[i]}
else{
rocx[i+1]<-rocx[i]+1
rocy[i+1]<-rocy[i]}
}
rocx<-as.matrix(rocx)
rocy<-as.matrix(rocy)
}
}

list(list=edgesRef,nameslist=edgesNames,rocx=rocx,rocy=rocy)
}

