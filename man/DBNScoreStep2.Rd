\name{DBNScoreStep2}
\alias{DBNScoreStep2}
\title{Full order dependence DAG G score matrix inference from a 1st
  order dependence score matrix S1}

\description{Given a time series dataset for \eqn{p} genes, a 1st
order dependence score matrix S1 (obtained with function DBNScoreStep1) and a threshold alpha1 for edge selection in matrix S1, this function infers the score of each edge of a Dynamic Bayesian Network (DAG G) describing full order dependencies between successive variables. This is the second step of the inference procedure described in the references. 1st step DBNScoreStep1 allows to reduce the number of potential edges, DBNScoreStep2 performs the last step selection. The
smallest score points out the most significant edge.  }

\usage{out<-DBNScoreStep2(S1,data,method='ls',alpha1,predPosition=NULL,
                                                     targetPosition=NULL)}

\arguments{\item{S1}{a matrix with \eqn{r} rows (=target genes) and
\eqn{d} columns (=predictor genes) containing score S1 (maximal p-value)
obtained with function DBNScoreStep1.}
        \item{data}{a matrix with \eqn{n} rows (=time points) and
\eqn{p} columns (=genes) containing the gene expression time series.}
	\item{method}{one of 'ls' (default), 'huber','tukey'. This specifies the regression method.}
	\item{alpha1}{Threshold for edge selection in the 1st order dependence score matrix S1. Edges having a score greater than alpha1 are pruned and quoted 'NA' is the resulting score matrix S2.}
	\item{predPosition}{To be specified if the number \eqn{d} of predictor
genes in score matrix S1 is lower than the number \eqn{p} of genes in the data: an
array included in \eqn{[1,p]} defining the position of the \eqn{d} predictor genes
in the data matrix (\eqn{n \times p}), default=NULL.}
	\item{targetPosition}{To be specified if the number \eqn{r} of
target genes in score  matrix S1 is lower than the number \eqn{p} of genes in the
data: an array included in \eqn{[1,p]} defining the position of the \eqn{r} target
genes in the data matrix (\eqn{n \times p}), default=NULL.}  }

\value{A matrix (\eqn{r} rows, \eqn{d} columns) containing the scores S2 obtained after the second step inference with the chosen M estimator. The score of the edges pruned after the first step inference is 'NA'.}

\author{L\`ebre Sophie (\url{http://www3.imperial.ac.uk/theoreticalgenomics/people/slebre}), 

Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet/}).}

\references{L\`ebre, S. 2007. Inferring Dynamic Bayesian Networks with low order dependencies. Preprint available at \url{http://hal.archives-ouvertes.fr/hal-00142109}.}


\keyword{graphs}
\keyword{ts}
\keyword{models}
\keyword{regression}
\keyword{robust}

\seealso{DBNScoreStep1, BuildEdges.}
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

split.screen(c(1,2))

#after Step 2
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

