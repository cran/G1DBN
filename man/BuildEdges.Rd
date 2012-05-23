\name{BuildEdges}
\alias{BuildEdges}
\title{Edges listing and evaluation}

\description{Given a score matrix, this function builds the list of  the edges of the associated network. The edges are ordered according to their scores.
The score matrix has been computed from a network inference algorithm (e.g. DBNScoreStep1 or DBNScoreStep2, Shrinkage, Lasso, ...). An optional threshold can be specified, as well as a maximal number of edges.
}


\usage{BuildEdges(score,threshold=1,nb=NULL,
                  targetNames=NULL,predNames=NULL,prec=3,dec=FALSE)
}

\arguments{\item{score}{matrix with \eqn{r} rows (=target genes) and
\eqn{d} columns (=predictor genes) containing the scores resulting from
an estimation procedure (e.g. DBNScoreStep1 or DBNScoreStep2, Shrinkage,
Lasso, ...).}
	\item{threshold}{An optional real setting the maximal value for edge selection, default=1.} 
	\item{nb}{An optional integer setting the maximal number of selected edges, default=NULL.} 
	\item{targetNames}{An optional array (\eqn{r}) giving a list of names for
the target genes, default=NULL.}
 	\item{predNames}{An optional array (\eqn{d}) giving a list of names for the predictor genes, default=NULL.} 
	\item{prec}{An optional integer setting the number of decimal places for score display, default=3.} 
	\item{dec}{boolean, FALSE if the smallest score points out the most significant edge, default=FALSE.} 
}

\value{A matrix containing a list of edges ordered according to the
score (First column: predictor, second column: target, third column:
corresponding score). Predictors and targets are referred to through the
names given by targetNames or predNames when specified.  }

\author{Lebre Sophie
  (\url{http://icube-bfo.unistra.fr/en/index.php/Sophie_Lebre}),
  
Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet}).}

\keyword{utilities}

\seealso{DBNScoreStep1, DBNScoreStep2, BuildNetwork}


\examples{
library(G1DBN)
## ======================================
## SIMULATING THE NETWORK

## number of genes
p <- 10
## the network - adjacency Matrix
MyNet <- SimulNetworkAdjMatrix(p,0.05,c(-1.5,-0.5,0.5,1.5))
MyNet

## ======================================
## SIMULATING THE TIME SERIES EXPERIMENTS

## number of time points
n <- 20
## initializing the B vector
B <- runif(p,-1,1)
## initializing the variance of the noise
sigmaEps <- runif(p,0.1,0.5)
## initializing the process Xt
X0 <- B + rnorm(p,0,sigmaEps*10)
## the times series process
Xn <- SimulGeneExpressionAR1(MyNet$A,B,X0,sigmaEps,n)


## ======================================
## NETWORK INFERENCE WITH G1DBN

## STEP 1 - The first step score matrix
S1 <- DBNScoreStep1(Xn, method='ls')

## Building the edges of the network inferred after Step1
alpha1 <- 0.5
G1 <- BuildEdges(S1$S1ls,threshold=alpha1,dec=FALSE)
G1
## STEP 2- The second step score matrix
S2 <- DBNScoreStep2(S1$S1ls, Xn, method='ls', alpha1)

## Building the edges of the network inferred after Step2
alpha2 <- 0.05
G2 <- BuildEdges(S2,threshold=alpha2,dec=FALSE)
G2

## Building the edges of the simulation Graph
Gsimul <- BuildEdges(MyNet$AdjMatrix,threshold=0,dec=TRUE)
Gsimul
}