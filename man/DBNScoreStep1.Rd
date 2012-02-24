\name{DBNScoreStep1}
\alias{DBNScoreStep1}
\title{First order dependence graph G(1) inference}

\description{Given a time series dataset for \eqn{p} genes, this
function infers a 1st order dependence score matrix S1 (\eqn{p \times p})
which contains the score of each edge of a Dynamic
Bayesian Network (DAG G(1)) describing first order dependencies between successives variables.  The smallest score points out the most significant edge for the 1st order dependence DAG G(1). The sets of both predictor and target genes can be reduced to different subsets of the \eqn{p} genes. DBNScoreStep1 is the first step of the estimation procedure described in the references. See function DBNScoreStep2 to perform the second step selection and infer a score matrix describing full order dependencies.}

\usage{DBNScoreStep1(data,method='ls',predPosition=NULL,targetPosition=NULL)
}

\arguments{\item{data}{a matrix with \eqn{n} rows (=time points) and
\eqn{p} columns (=genes) containing the gene expression time series.}
	\item{method}{currently M estimation with either LS, Tukey
	  bisquare or Huber estimator:
	  
	  \code{c('ls','tukey','huber')}, default='ls'.}
	\item{predPosition}{To be specified to reduce the set of
possible predictor genes to a subset of \eqn{d<p} genes: an array
included in \eqn{[1,p]} defining the position of the d predictor genes
in the data matrix (\eqn{n \times p}), default=NULL.}
	\item{targetPosition}{To be specified to reduce the set of
possible target genes to a subset of \eqn{r<p} genes: an array included
in \eqn{[1,p]} defining the position of the \eqn{r} target genes in the
data matrix (\eqn{n \times p}), default=NULL.}  }

\value{A list with out$S1ls a matrix with \eqn{\mathrm{min}(r,p)} rows
(=target genes) and \eqn{\mathrm{min}(d,p)} columns (=predictor genes)
containing the scores S1 obtained with least square estimator, out$S1huber a matrix
containing scores S1 obtained with Huber estimator, out$S1tukey a
matrix containing scores S1 obtained with Tukey bisquare (or biweight)
estimator.(out$S1ls[i,j] is the score for the edge \eqn{j \leftarrow i} pointing
out from predictor \eqn{j} toward target \eqn{i}.)}

\references{Lebre, S. 2009. Inferring dynamic bayesian network with low order independencies, Statistical Applications in Genetics and Molecular Biology, 2009: Vol. 8: Iss. 1, Article 9.}

\author{Lebre Sophie
  (\url{http://newlsiit.u-strasbg.fr/bfo_en/index.php/Sophie_Lebre}),
  
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
id<-c(60, 141, 260, 333, 365, 424, 441, 512, 521, 578, 789, 799)
names<-c("carbohydrate/sugar transporter","ATGPX2","putative integral
membrane prot" ,
"AT3G05900", "At3g27350", "At1g16720","ATISA3/ISA3","AT4G32190",
"catalase", "plasma membrane intrinsic prot", "At4g16146", "DPE2")

## compute score S1 
out<-DBNScoreStep1(data,method='ls', targetPosition=id,predPosition=id)
round(out$S1ls,2)


## Threshold for the selection of the edges after Step 1
alpha1=0.5
## Build the edges with id as label
edgesG1id<-BuildEdges(score=out$S1ls,threshold=alpha1,
                       targetNames=id,predNames=id,prec=6)
## Build the edges with names as label
edgesG1names<-BuildEdges(score=out$S1ls,threshold=alpha1,
                         targetNames=names,predNames=names,prec=6)
edgesG1id[1:15,]
edgesG1names[1:15,]


## compute score S2 from S1 
S2<-DBNScoreStep2(out$S1ls,data,method='ls',alpha1=alpha1,
                  predPosition=id,targetPosition=id)
S2

## Threshold for the selection of the edges after Step 2
alpha2=0.05
## Build the edges with id as label
edgesG2id<-BuildEdges(score=S2,threshold=alpha2,
                      targetNames=id,predNames=id,prec=6)
## Build the edges with names as label
edgesG2names<-BuildEdges(score=S2,threshold=alpha2,
                         targetNames=names,predNames=names,prec=6)
edgesG2id
edgesG2names


## As the number of genes is reduced to 10 here, this results slightly differ
## from the results obtained in the paper (Lebre, 2009) cited in References.


## ======================================
## PLOTTING THE RESULTS...
## ______________________________________
\dontrun{
## The Inferred Nets
## -----------------

## Nodes coordinates are calculated according to the global structure of the graph
all_parents=c(edgesG1id[,1], edgesG2id[,1])
all_targets=c(edgesG1id[,2], edgesG2id[,2])
posEdgesG1=1:dim(edgesG1id)[1]
posEdgesG2=(dim(edgesG1id)[1]+1):length(all_targets)

## Global network with all the edges
netAll = graph.edgelist(cbind(as.character(all_parents),as.character(all_targets )))

## Nodes coordinates
nodeCoord=layout.fruchterman.reingold(netAll)


split.screen(c(1,2))

# after Step 1
screen(1)
# set the edges list
netG1 = graph.edgelist(cbind(as.character(edgesG1id[,1]),as.character(edgesG1id[,2])))
# set the object for plotting the network with global coordinates of all nodes
G1toPlot=delete.edges(netAll, E(netAll)[posEdgesG2] )
# plot the network
plot(G1toPlot, layout=nodeCoord, vertex.label = 
get.vertex.attribute(G1toPlot , name="name"), edge.arrow.size = 0.2,
main="G1DBN Inferred network:\n Step 1")

# after Step 2
screen(2)
# set the edges list

netG2 = graph.edgelist(cbind(as.character(edgesG2id[,1]),as.character(edgesG2id[,2])))
# set the object for plotting the network with global coordinates of all nodes
G2toPlot=delete.edges(netAll, E(netAll)[posEdgesG1] )
# plot the network
plot(G2toPlot, layout=nodeCoord, vertex.label = 
get.vertex.attribute(G2toPlot , name="name"),edge.arrow.size = 0.2,
main="G1DBN Inferred network:\n Step 2")

close.screen(all = TRUE)
}
}


