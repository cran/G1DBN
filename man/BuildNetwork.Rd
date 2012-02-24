\name{BuildNetwork}
\alias{BuildNetwork}
\title{Network object creation}

\description{Given a list of scored edges and a REQUIRED vector of labels
(e.g, \code{1:p}), this function builds an object "Network". This is useful for exporting a network described by a list of edges to a network described by an adjacency matrix.}

\usage{BuildNetwork(Edges,Labels,nonedges.val=NA)
}

\arguments{\item{Edges}{a \eqn{n \times 3} matrix (each line contains a
couple of vertices plus the associated score)}
  \item{Labels}{a vector of labels for the vertices}
  \item{nonedges.val}{optional. Value attributed to the not existing edges in the generated score matrix out$Score, default=NA}
}

\value{a list that contains out$Vertices$Num the number of vertices, out$Vertices$Labels a vector of labels of the vertices, out$Vertices$Connected a vector of the connected vertices, out$Edges$Prop the proportion of edges, out$Edges$Num the number of edges, the graph of the network (out$AdjMatrix an adjacency matrix and out$Edges$List  a list of edges) and out$Score a score matrix.}

\seealso{BuildEdges}

\author{Lebre Sophie (\url{http://newlsiit.u-strasbg.fr/bfo_en/index.php/Sophie_Lebre}),

Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet}).}


\keyword{datagen}

\examples{
library(G1DBN)

## ========================
## SIMULATING THE NETWORK
## ________________________

## number of genes
p <- 10
## number of time points
n <- 20
## proportion of genes
geneProp <- 0.05
## the network - adjacency Matrix
MyNet <- SimulNetworkAdjMatrix(p,geneProp,c(-1.5,-0.5,0.5,1.5))

cat("\n==========================================\n")
cat("SIMULATION\n")

## ======================================
## SIMULATING THE TIME SERIES EXPERIMENTS
## ______________________________________
##
## Autoregressive model
##
cat("Time series experiments with")

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
## ______________________________________
##
cat("\n==========================================\n")
cat("NETWORK INFERENCE\n")
cat("Using a Dynamic Bayesian Network model\n")

## STEP 1
## ------
cat("STEP 1...\n")
S1 <- DBNScoreStep1(Xn, method='ls')

## STEP 2
## ------
cat("STEP 2...\n")
alpha1=0.5
S2 <- DBNScoreStep2(S1$S1ls, data=Xn, method='ls', alpha1=alpha1)

## ======================================
## POST TREATMENTS

## building the inferred Graph
G1 <- BuildEdges(S1$S1ls,threshold=alpha1,dec=FALSE)

## encoding as the adjancecy matrix graph
Step1InferredNet <- BuildNetwork(G1,1:p)
Step1InferredNet

#Step 2
alpha2=0.05
G2 <- BuildEdges(S2,threshold=alpha2,dec=FALSE)
Step2InferredNet <- BuildNetwork(G2,1:p)


## ======================================
## PLOTTING THE RESULTS...
## ______________________________________
\dontrun{
cat("\n==========================================\n")
cat("SUMMARY\n")
cat("Plotting the results...\n")
split.screen(c(1,3))


## The Original graph and data
## ---------------------------
# set the edges list of the simulated network
G0 <- BuildEdges(MyNet$AdjMatrix,threshold=0.9,dec=TRUE)

## Nodes coordinates are calculated according to the global structure of the network
all_parents=c(G0[,1],G1[,1], G2[,1])
all_targets=c(G0[,2],G1[,2], G2[,2])
posEdgesG0=1:dim(G0)[1]
posEdgesG1=(dim(G0)[1]+1):(dim(G0)[1]+dim(G1)[1])
posEdgesG2=(dim(G0)[1]+dim(G1)[1]+1):length(all_parents)

## Global network with all the edges
netAll =
graph.edgelist(cbind(as.character(all_parents),as.character(all_targets)))

## Nodes coordinates
nodeCoord=layout.fruchterman.reingold(netAll)


#after Step 1
screen(1)
# set the edges list
netG1 = graph.edgelist(cbind(as.character(G1[,1]),as.character(G1[,2])))
# set the object for plotting the network with global coordinates of all nodes
G1toPlot=delete.edges(netAll, E(netAll)[c(posEdgesG0,posEdgesG2)-1] )
# plot the network
plot(G1toPlot, layout=nodeCoord, vertex.label = 
get.vertex.attribute(G1toPlot, name="name"), edge.arrow.size = 0.2,
main="G1DBN Inferred network:\n Step 1")

#after Step 2
screen(2)
# set the edges list
netG2 = graph.edgelist(cbind(as.character(G2[,1]),as.character(G2[,2])))
# set the object for plotting the network with global coordinates of all nodes
G2toPlot=delete.edges(netAll, E(netAll)[c(posEdgesG0,posEdgesG1)-1 ] )
# plot the network
plot(G2toPlot, layout=nodeCoord, vertex.label = 
get.vertex.attribute(G2toPlot, name="name"),edge.arrow.size = 0.2,
main="G1DBN Inferred network:\n Step 2")

screen(3)
net0 = graph.edgelist(cbind(as.character(G0[,1]),as.character(G0[,2])))
# set the object for plotting the network with global coordinates of all nodes
G0toPlot=delete.edges(netAll, E(netAll)[c(posEdgesG1,posEdgesG2)-1] )
plot(G0toPlot, layout=nodeCoord, vertex.label =
get.vertex.attribute(G0toPlot, name="name"), edge.arrow.size = 0.2,
main="Simulated network:")

close.screen(all = TRUE)
}
cat("")
cat("\nDONE !\n")
}
