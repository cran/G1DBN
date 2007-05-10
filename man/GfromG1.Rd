\name{GfromG1}
\alias{GfromG1}
\title{Full order dependence DAG G inference}

\description{This function infers the scores of each potential edge of a Dynamic Bayesian Network (defined by full order dependence DAG G) from a score matrix S1 (obtained with function inferG1)
which describes the score of the edges for the 1st order dependence DAG G(1). 
This is the second step of the inference procedure described in the references: 1st step inferG1 allows
to reduce the number of potential edges, GfromG1 performs the last step selection. The smallest score points out the most significant edge.
}

\usage{out<-GfromG1(S1,data,method='ls',alpha1,alpha2=1,predictor=NULL,target=NULL)
}

\arguments{\item{S1}{a matrix with r rows (=target genes) and d columns (=predictor genes) containing score S1 (maximal p-value) obtained with function inferG1.} 
	\item{data}{the matrix with n rows (=time points) and p columns (=genes) containing the corresponding gene expression time series.} 
	\item{method}{currently M estimation with either LS, Tukey bisquare or Huber estimator, c('ls','tukey','huber'), default='ls'.} 
	\item{alpha1}{Step 1 threshold for edge selection in G(1).}
	\item{alpha2}{Step 2 threshold for edge selection in G, default=1.}
	\item{predictor}{To be specified if the possible predictor genes should be restricted to a subset of d<p genes: an array included in [1,p] defining the position of the d predictor genes in the data matrix, default=NULL.} 
	\item{target}{To be specified if the possible target genes should be reduced to a subset of r<p genes: an array included in [1,p] defining the position of the r target genes in the data matrix, default=NULL.} 
}

\value{A list with out\$S2 a matrix (r rows, d columns) containing the scores obtained with the chosen M estimator,
out\$maxInG1 the maximal number of parent in DAG G(1) for the chosen threshold alpha1,
out\$NbG1 the number of edges in 1st order dependence DAG G(1),
out\$NbG the number of edges in full order dependence DAG G.
}

\author{Lebre Sophie (\url{http://stat.genopole.cnrs.fr/~slebre}).}


\references{Lebre, S. 2007. Inferring Dynamic Bayesian Networks with low order dependencies. Preprint available at \url{http://hal.archives-ouvertes.fr/hal-00142109}.}


\keyword{graphs}
\keyword{ts}
\keyword{models}
\keyword{regression}
\keyword{robust}

\seealso{inferG1, edges.}

\examples{# load G1DBN Library
library(G1DBN)

data(arth800line)
data<-arth800line
id<-c(60, 141, 260, 333, 365, 424, 441, 512, 521, 578, 789, 799)

# compute score S1 
pmaxG1<-inferG1(data,ls=TRUE,tukey=FALSE,huber=FALSE,predictor=id,target=id)
round(pmaxG1$ls,2)

# compute score S2 from S1 
GwithLS<-GfromG1(pmaxG1$ls, data, method='ls', alpha1=0.1, alpha2=0.01, predictor=NULL,
target=NULL)
GwithLS

resG<-edges(score=GwithLS$S2,targetNames=id,predNames=id,validMat=NULL,roc=FALSE,
threshold=0.001,nb=NULL,prec=6)
resG$nameslist

# As the number of genes is reduced to 10 here, this results slightly differ 
# from the results obtained in the Preprint cited in References.

}



