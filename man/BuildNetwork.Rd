\name{BuildNetwork}
\alias{BuildNetwork}
\title{Network object creation}

\description{Given a list of scored edges and a REQUIRED vector of labels
(e.g, \code{1:p}), this function builds an object "Network". This is useful for exporting a network described by a list of edges to a network described by an adjacency matrix.}

\usage{out <- BuildNetwork(Edges,Labels,nonedges.val=NA)
}

\arguments{\item{Edges}{a \eqn{n \times 3} matrix (each line contains a
couple of vertices plus the associated score)}
  \item{Labels}{a vector of labels for the vertices}
  \item{nonedges.val}{optional. Value attributed to the not existing edges in the generated score matrix out\$Score, default=NA}
}

\value{a list that contains out\$Vertices\$Num the number of vertices, out\$Vertices\$Labels a vector of labels of the vertices, out\$Vertices\$Connected a vector of the connected vertices, out\$Edges\$Prop the proportion of edges, out\$Edges\$Num the number of edges, the graph of the network (out\$AdjMatrix an adjacency matrix and out\$Edges\$List  a list of edges) and out\$Score a score matrix.}

\seealso{BuildEdges}

\author{L\`ebre Sophie (\url{http://www3.imperial.ac.uk/theoreticalgenomics/people/slebre}), 

Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet/}).}

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
cat("SIMULATION\n\n")

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
## NETWORK INFERENCE WITH DBN
## ______________________________________
##
cat("\n==========================================\n")
cat("NETWORK INFERENCE\n\n")
cat("Using a Dynamic Bayesian Network model\n\n")

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
G <- BuildEdges(S2,threshold=alpha2,dec=FALSE)
Step2InferredNet <- BuildNetwork(G,1:p)


## ======================================
## PLOTTING THE RESULTS...
## ______________________________________

cat("\n==========================================\n")
cat("SUMMARY\n\n")
cat("Plotting the results...\n")


## The Original graph and data
## ---------------------------
attach(MyNet)
pos <- gplot(abs(AdjMatrix), vertex.cex=1.5, diag=TRUE,
      displaylabel=TRUE, usecurv=TRUE,
      boxed.label=FALSE, main="Simulated network")
detach(MyNet)

## The Inferred Nets
## -----------------
#after Step 2

split.screen(c(1,2))
screen(1)
attach(Step2InferredNet)
gplot(abs(AdjMatrix), vertex.cex=1.5, diag=TRUE, coor=pos,
      displaylabel=TRUE, usecurv=TRUE,
      boxed.label=FALSE, main="Inferred network - Step 2")
detach(Step2InferredNet)

#after Step 1
screen(2)
attach(Step1InferredNet)
gplot(abs(AdjMatrix), vertex.cex=1.5, diag=TRUE, coor=pos,
      displaylabel=TRUE, usecurv=TRUE, 
      boxed.label=FALSE, main="Inferred network - Step 1")
detach(Step1InferredNet)

close.screen(all = TRUE)

cat("")
cat("\nDONE !\n")
}
