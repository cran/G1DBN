\name{inferG1}
\alias{inferG1}
\title{First order dependence graph G(1) inference}

\description{This function infers the score S1 of each potential edge of a Dynamic Bayesian Network (defined by full order dependence DAG G) by considering 1st order dependencies. The smallest score points out the most significant edge for the 1st order dependence DAG G(1). inferG1 is the 1st step of the estimation procedure described in the references.}

\usage{out<-inferG1(data,ls=T,huber=F,tukey=F,kendall=F,predictor=NULL,target=NULL)
}

\arguments{
	\item{data}{a matrix with n rows (=time points) and p columns (=genes) containing the gene expression time series.} 
	\item{ls}{boolean, TRUE to obtain a score matrix computed by using Least Square estimator, default=TRUE.} 
	\item{huber}{boolean, TRUE to obtain a score matrix computed by using Huber estimator, default=FALSE.}
	\item{tukey}{boolean, TRUE to obtain a score matrix computed by using Tukey bisquare (or biweight) estimator, default=FALSE.} 
	\item{predictor}{To be specified if the possible predictor genes should be reduced to a subset of d<p genes: an array included in [1,p] defining the position of the d predictor genes in the data matrix,
	default=NULL.} 
	\item{target}{To be specified if the possible target genes should be reduced to a subset of r<p genes: an array included in [1,p] defining the
	position of the r target genes in the data matrix, default=NULL.} 
}

\value{A list with out\$ls a matrix with min(r,p) columns (=target genes) and min(d,p) rows (=predictor genes) containing the scores S1 obtained with least square estimator (out\$ls[i,j] is the score for the
edge j -> i), 
out\$huber a matrix containing scores S1 obtained with Huber estimator, 
out\$tukey a matrix containing scores S1 obtained with Tukey bisquare (or biweight) estimator. 
}

\references{Lebre, S. 2007. Inferring Dynamic Bayesian Networks with low order
 dependencies. Preprint available at \url{http://hal.archives-ouvertes.fr/hal-00142109}.}

\author{Lebre Sophie (\url{http://stat.genopole.cnrs.fr/~slebre}).}

\seealso{GfromG1, edges.}

\keyword{graphs}
\keyword{ts}
\keyword{models}
\keyword{regression}
\keyword{robust}

\note{For a large number of target genes, it may be interesting to parallel run the procedure inferG1 for each target gene.}

\examples{# load G1DBN Library
library(G1DBN)

data(arth800line)
data<-arth800line
id<-c(60, 141, 260, 333, 365, 424, 441, 512, 521, 578, 789, 799)

# compute score S1 
pmaxG1<-inferG1(data,ls=TRUE,tukey=FALSE,huber=FALSE,predictor=id,target=id)
round(pmaxG1$ls,2)

resG1<-edges(score=pmaxG1$ls,targetNames=id,predNames=id,validMat=NULL,roc=FALSE,
threshold=0.001,nb=NULL,prec=6)
resG1$nameslist

# compute score S2 from S1 
GwithLS<-GfromG1(pmaxG1$ls,data,method='ls',alpha1=0.1,alpha2=0.01,predictor=NULL,
target=NULL)
GwithLS


resG<-edges(score=GwithLS$S2,targetNames=id,predNames=id,validMat=NULL,roc=FALSE,
threshold=0.001,nb=NULL,prec=6)
resG$nameslist
}
