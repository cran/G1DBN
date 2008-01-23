\name{SimulNetworkAdjMatrix}
\alias{SimulNetworkAdjMatrix}
\title{Network object generation}

\description{This function builds a object "Network" by simulating a
matrix of valued adjacencies from a number of vertices, a proportion
of edges and the range of the uniform distribution that is used to build
the adjacency matrix. An optional vector of labels may be given.}

\usage{out<-SimulNetworkAdjMatrix(Num,EdgesProp,Range,Labels=1:Num)}

\arguments{\item{Num}{number of genes}
  \item{EdgesProp}{edges proportion in the network}
  \item{Range}{vector with 4 elements specifying range values for the adjacency matrix generation (minimum negative value, maximum negative value, minimum positive value, maximum positive value)}
  \item{Labels}{an optional vector of labels for the edges}
}

\value{a list that contains out\$Vertices\$Num the number of vertices,  out\$Vertices\$Labels a vector of labels of the vertices, out\$Vertices\$Regulated a vector of the regulated vertices, out\$Edges\$Prop the proportion of edges, out\$Edges\$Num the number of edges, out\$AdjMatrix an adjacency matrix  (binary) and out\$A a valued adjacency matrix.}


\author{L\`ebre Sophie (\url{http://www3.imperial.ac.uk/theoreticalgenomics/people/slebre/}),

Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet}).}


\seealso{SimulGeneExpressionAR1, BuildEdges}


\keyword{datagen}

\examples{
library(G1DBN)
## number of genes
p <- 10
## the network - adjacency Matrix
MyNet <- SimulNetworkAdjMatrix(p,0.05,c(-1,0,0,1))
MyNet

## initializing the B vector
B <- runif(p,0,0.5)
## initializing the variance of the noise
sigmaEps <- runif(p,0.1,0.8)
## initializing the process Xt
X0 <- B + rnorm(p,0,sigmaEps*10)
## number of time points
n <- 20

## the AR(1) times series process
Xn <- SimulGeneExpressionAR1(MyNet$AdjMatrix,B,X0,sigmaEps,n)
}
