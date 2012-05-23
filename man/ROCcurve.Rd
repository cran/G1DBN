\name{ROCcurve}
\alias{ROCcurve}
\title{ROC curves computation}

\description{Given a score matrix and a validation matrix, this function
allows to compute the corresponding ROC curve by returning a list with
the respective \eqn{x} and \eqn{y} coordinates of the ROC curve. The
score matrix has been computed from a network inference algorithm (e.g. DBNScoreStep1 or DBNScoreStep2, Shrinkage, Lasso, ...).  }

\usage{ROCcurve(score,validMat,dec=FALSE)
}

\arguments{\item{score}{matrix with \eqn{r} rows (=target genes) and
\eqn{d} columns (=predictor genes) containing the scores resulting from
an estimation procedure (e.g. DBNScoreStep1 or DBNScoreStep2, Shrinkage,
Lasso, ...).}
	\item{validMat}{An optional matrix specifying the validated edges (1 if an edge is validated, 0 otherwise).}
	\item{dec}{boolean, FALSE if the smallest score points out the most significant edge, default=FALSE.} 
}

\value{A list with out$x and out$y contain respctively the \eqn{x} and the \eqn{y} coordinates of the ROC curve.} 

\keyword{utilities}

\author{Lebre Sophie
   (\url{http://icube-bfo.unistra.fr/en/index.php/Sophie_Lebre}),
 
  
Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet}).}

\seealso{DBNScoreStep1, DBNScoreStep2, PRcurve.}

\examples{
## generate the validation matrix
## number of genes
p <- 20
## the network - adjacency Matrix
MyNet <- SimulNetworkAdjMatrix(p,0.05,c(-1,0.5,0.5,1))

## generate the time series 
## initializing the B vector
B <- runif(p,-1,1)
## initializing the variance of the noise
sigmaEps <- runif(p,0.05,0.5)
## initializing the process Xt
X0 <- B + rnorm(p,0,sigmaEps*10)
## number of time points
n <- 20

## the AR(1) times series process
Xn <- SimulGeneExpressionAR1(MyNet$A,B,X0,sigmaEps,n)


## compute score S1 
out<-DBNScoreStep1(Xn)
roc1<-ROCcurve(score=out$S1ls,validMat=abs(MyNet$AdjMatrix)>0,dec=FALSE)

## compute score S2 from S1 
## depending on the generated data, the threshold alpha1 has to be chosen differently.
alpha1=0.8
S2<-DBNScoreStep2(S1=out$S1ls,data=Xn,alpha1=alpha1)
roc2_0.8<-ROCcurve(score=S2,validMat=abs(MyNet$AdjMatrix)>0,dec=FALSE)

alpha1=0.6
S2<-DBNScoreStep2(S1=out$S1ls,data=Xn,alpha1=alpha1)
roc2_0.6<-ROCcurve(score=S2,validMat=abs(MyNet$AdjMatrix)>0,dec=FALSE)

alpha1=0.4
S2<-DBNScoreStep2(S1=out$S1ls,data=Xn,alpha1=alpha1)
roc2_0.4<-ROCcurve(score=S2,validMat=abs(MyNet$AdjMatrix)>0,dec=FALSE)

alpha1=0.2
S2<-DBNScoreStep2(S1=out$S1ls,data=Xn,alpha1=alpha1)
roc2_0.2<-ROCcurve(score=S2,validMat=abs(MyNet$AdjMatrix)>0,dec=FALSE)


TP=sum(abs(MyNet$AdjMatrix)>0)
FN=p^2-TP

plot(roc1$x/FN,roc1$y/TP,type="l", main="ROC curve after both Step1 and Step2",
     ylab="True Positive Rate", xlab="False Negative Rate",lwd=2,lty=2)
lines(roc2_0.8$x/FN, roc2_0.8$y/TP, col=3,lwd=2)
lines(roc2_0.6$x/FN, roc2_0.6$y/TP, col=4,lwd=2)
lines(roc2_0.4$x/FN, roc2_0.4$y/TP, col=5,lwd=2)
lines(roc2_0.2$x/FN, roc2_0.2$y/TP, col=6,lwd=2)
lines(0:1,c(0,0),lty=3)
lines(0:1,c(1,1),lty=3)
lines(c(0,0),0:1,lty=3)
lines(c(1,1),0:1,lty=3)
leg=c("Step 1", "Step 2 (alpha1=0.8)", "Step 2 (alpha1=0.6)",
      "Step 2 (alpha1=0.4)", "Step 2 (alpha1=0.2)")
legend(0.568,0.265, leg, lty=c(2,1,1,1,1), col=c(1,3,4,5,6),lwd=array(2,5))
}

