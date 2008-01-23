\name{DBNScoreStep1}
\alias{DBNScoreStep1}
\title{First order dependence graph G(1) inference}

\description{Given a time series dataset for \eqn{p} genes, this
function infers a 1st order dependence score matrix S1 (\eqn{p \times p})
which contains the score of each edge of a Dynamic
Bayesian Network (DAG G(1)) describing first order dependencies between successives variables.  The smallest score points out the most significant edge for the 1st order dependence DAG G(1). The sets of both predictor and target genes can be reduced to different subsets of the \eqn{p} genes. DBNScoreStep1 is the first step of the estimation procedure described in the references. See function DBNScoreStep2 to perform the second step selection and infer a score matrix describing full order dependencies.}

\usage{out <- DBNScoreStep1(data,method='ls',predPosition=NULL,targetPosition=NULL)
}

\arguments{\item{data}{a matrix with \eqn{n} rows (=time points) and
\eqn{p} columns (=genes) containing the gene expression time series.}
	\item{method}{currently M estimation with either LS, Tukey
	  bisquare or Huber estimator, \code{c('ls','tukey','huber')},
	  default='ls'.}
	\item{predPosition}{To be specified to reduce the set of
possible predictor genes to a subset of \eqn{d<p} genes: an array
included in \eqn{[1,p]} defining the position of the d predictor genes
in the data matrix (\eqn{n \times p}), default=NULL.}
	\item{targetPosition}{To be specified to reduce the set of
possible target genes to a subset of \eqn{r<p} genes: an array included
in \eqn{[1,p]} defining the position of the \eqn{r} target genes in the
data matrix (\eqn{n \times p}), default=NULL.}  }

\value{A list with out\$S1ls a matrix with \eqn{\mathrm{min}(r,p)} rows
(=target genes) and \eqn{\mathrm{min}(d,p)} columns (=predictor genes)
containing the scores S1 obtained with least square estimator, out\$S1huber a matrix
containing scores S1 obtained with Huber estimator, out\$S1tukey a
matrix containing scores S1 obtained with Tukey bisquare (or biweight)
estimator.(out\$S1ls[i,j] is the score for the edge \eqn{j \leftarrow i} pointing
out from predictor \eqn{j} toward target \eqn{i}.)}

\references{L\`ebre, S. 2007. Inferring Dynamic Bayesian Networks with low order
 dependencies. Preprint available at
 \url{http://hal.archives-ouvertes.fr/hal-00142109}.}

\author{L\`ebre Sophie (\url{http://www3.imperial.ac.uk/theoreticalgenomics/people/slebre/}),

Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet}).}

\seealso{DBNScoreStep2, BuildEdges, PRcurve.}

\keyword{graphs}
\keyword{ts}
\keyword{models}
\keyword{regression}
\keyword{robust}

\note{For a large number of target genes, it is of interest to parallel
run the procedure DBNScoreStep1 for each target gene by running \eqn{p}
the following jobs for \eqn{i=1\dots p},

\code{outi <- DBNScoreStep1(data, target=i)}.}

\examples{
## load G1DBN Library
library(G1DBN)

data(arth800line)
data<-as.matrix(arth800line)
idx<-c(60, 141, 260, 333, 365, 424, 441, 512, 521, 578, 789, 799)
names<-c("carbohydrate/sugar transporter","ATGPX2","putative integral
membrane prot" ,
"AT3G05900", "At3g27350", "At1g16720","ATISA3/ISA3","AT4G32190",
"catalase", "plasma membrane intrinsic prot", "At4g16146", "DPE2")

## compute score S1 
out<-DBNScoreStep1(data,method='ls', targetPosition=idx,predPosition=idx)
round(out$S1ls,2)


alpha1=0.5
edgesG1idx<-BuildEdges(score=out$S1ls,threshold=alpha1,
                       targetNames=idx,predNames=idx,prec=6)
edgesG1names<-BuildEdges(score=out$S1ls,threshold=alpha1,
                         targetNames=names,predNames=names,prec=6)
edgesG1idx[1:15,]
edgesG1names[1:15,]


## compute score S2 from S1 
S2<-DBNScoreStep2(out$S1ls,data,method='ls',alpha1=alpha1,
                  predPosition=idx,targetPosition=idx)
S2

alpha2=0.05
edgesGidx<-BuildEdges(score=S2,threshold=alpha2,
                      targetNames=idx,predNames=idx,prec=6)
edgesGidx

## As the number of genes is reduced to 10 here, this results slightly differ
## from the results obtained in the Preprint cited in References.


## encoding as the adjacency matrix graph
Step1InferredNet <- BuildNetwork(edgesG1idx,idx)

## encoding as the adjancecy matrix graph
Step2InferredNet <- BuildNetwork(edgesGidx,idx)

## The Inferred Nets
## -----------------

#after Step 2
split.screen(c(1,2))
screen(1)
attach(Step2InferredNet)
pos<-gplot(t(AdjMatrix), vertex.cex=1.5, diag=TRUE,
      displaylabel=TRUE, usecurv=TRUE, label=names,
      boxed.label=FALSE, main="Inferred network - Step 2")
detach(Step2InferredNet)

#after Step 1
screen(2)
attach(Step1InferredNet)
gplot(t(AdjMatrix), vertex.cex=1.5, diag=TRUE, coord=pos,
      displaylabel=TRUE, usecurv=TRUE, label=names,
      boxed.label=FALSE, main="Inferred network - Step 1")
detach(Step1InferredNet)
close.screen(all = TRUE)

}