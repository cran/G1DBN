\name{edges}
\alias{edges}
\title{Edges listing and evaluation}

\description{This function allows to order the edges by decreasing likelihood. ROC curves can also be computed if a validation matrix is specified. 
}


\usage{
out<-edges(score,predNames=NULL,targetNames=NULL,validMat=NULL,roc=FALSE,
threshold=1, nb=NULL, prec=3)
}

\arguments{\item{score}{matrix with r columns (=target genes) and d rows (=predictor genes) containing the scores resulting from an estimation procedure (either inferG1 or GfromG1).} 
	\item{predNames}{An optional array (d) giving a list of names for the predictor genes, default=NULL.} 
	\item{targetNames}{An optional array (r) giving a list of names for
the target genes, default=NULL.} 
	\item{validMat}{An optional matrix specifying the validated edges (1 if an edge is validated, 0 otherwise).}
	\item{roc}{Boolean, TRUE wether ROC curves should be computed, default=NULL. This option require to specify a validation matrix validMat.} 
	\item{threshold}{An optional real setting the maximal value for edge selection, default=1.} 
	\item{nb}{An optional integer setting the maximal number of selected edges, default=NULL.} 
	\item{prec}{An optional integer setting the number of decimal places for score display, default=3.} 
}

\value{A list with out\$ref and out\$names matrices containing a list of edges ordered by decreasing likelihood (First column: parent, second column: child, third column: the corresponding score). If a validation matrix validMat is specified, a 4th column indicates wether the edge is validated (1 for a validated edge, 0 otherwise).

out\$ref lists edges according to the index value of each variable and out\$names uses the names from arrays predNames and targetNames. 

out\$rocx and out\$rocy contain the coordinates for plotting ROC curves if a validation matrix is specified.} 

\keyword{utilities}

\seealso{inferG1, GfromG1.}

\author{Lebre Sophie (\url{http://stat.genopole.cnrs.fr/~slebre}).}

\examples{#generate AR(1) time series 
AR<-simulAR1(p=10,n=50,edgeProp=0.02,minA=0.5,maxA=1.5,minB=0,maxB=1,minSig=0.1,maxSig=0.8)
AR

# compute score S1
pmaxG1<-inferG1(AR$data, ls=TRUE, tukey=FALSE, huber=FALSE, predictor=NULL, target=NULL)
pmaxG1$ls
# compute score S2 from S1
Gls<-GfromG1(pmaxG1$ls, AR$data, method='ls',alpha1=0.2,alpha2=1,predictor=NULL,
target=NULL)
Gls$score

### Results (The validation matrix is obtained from the simulation matrix AR$A). 
resG1<-edges(score=pmaxG1$ls, targetNames=NULL, predNames=NULL,
validMat = (abs(AR$A)>0)*1, roc=TRUE, threshold=1, nb=NULL, prec=4)
resG<-edges(score=Gls$S2, targetNames=NULL, predNames=NULL,
validMat=(abs(AR$A)>0)*1, roc=TRUE, threshold=1, nb=NULL, prec=4)

# Edges list
resG1$list
resG$list

# ROC curve
plot(resG1$rocx/max(resG1$rocx),resG1$rocy/max(resG1$rocy), type="l",
 xlab="False Positive", ylab="True Positive", main="ROC curve")
lines(resG$rocx/max(resG$rocx), resG$rocy/max(resG$rocy), col=2,lty=2)
leg=c("Step 1", "Step 2")
legend(0.8,0.2, leg, lty=c(1,2), col=c(1,2))
}

