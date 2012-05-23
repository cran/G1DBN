\name{DBNScoreStep2}
\alias{DBNScoreStep2}
\title{Full order dependence DAG G score matrix inference from a 1st
  order dependence score matrix S1}

\description{Given a time series dataset for \eqn{p} genes, a 1st
order dependence score matrix S1 (obtained with function DBNScoreStep1) and a threshold alpha1 for edge selection in matrix S1, this function infers the score of each edge of a Dynamic Bayesian Network (DAG G) describing full order dependencies between successive variables. This is the second step of the inference procedure described in the references. 1st step DBNScoreStep1 allows to reduce the number of potential edges, DBNScoreStep2 performs the last step selection. The
smallest score points out the most significant edge.  }

\usage{DBNScoreStep2(S1,data,method='ls',alpha1,predPosition=NULL,
                                          targetPosition=NULL)
}

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

\references{Lebre, S. 2009. Inferring dynamic bayesian network with low order independencies, Statistical Applications in Genetics and Molecular Biology, 2009: Vol. 8: Iss. 1, Article 9.}

\author{Lebre Sophie
 (\url{http://icube-bfo.unistra.fr/en/index.php/Sophie_Lebre}),
 
  
Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet}).}



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
get.vertex.attribute(G1toPlot , name="name"),edge.arrow.size = 0.2,
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

