\name{SimulGeneExpressionAR1}
\alias{SimulGeneExpressionAR1}
\title{First order multivariate Auto-Regressive time series generation}

\description{This function generates multivariate time series according to the following first order Auto-Regressive process,
\deqn{X(t)= A X(t-1) + B + \varepsilon(t),}
where \eqn{\varepsilon(t)} follows a zero-centered multivariate gaussian distribution whose variance matrix \eqn{S} is diagonal.
}

\usage{out<-SimulGeneExpressionAR1(A,B,X0,SigmaEps,n)
}

\arguments{\item{A}{a matrix (\eqn{p \times p})}
	\item{B}{a column vector (\eqn{p \times 1})}
        \item{X0}{a column vector (\eqn{p \times 1}) containing the values of the process at time 0}
        \item{SigmaEps}{a column vector (\eqn{p \times 1}) containing the values of the diagonal of covariance matrix S}
	\item{n}{the desired length of the time serie.}
}

\value{A matrix, with \eqn{n} rows (=length) and \eqn{p} columns
  (=dimension), containing the generated time series,
}

\seealso{SimulNetworkAdjMatrix}

\author{L\`ebre Sophie (\url{http://www3.imperial.ac.uk/theoreticalgenomics/people/slebre}), 

Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet/}).}


\keyword{datagen}
\keyword{ts}

\examples{
library(G1DBN)
## number of genes
p <- 20
## the network - adjacency Matrix
MyNet <- SimulNetworkAdjMatrix(p,0.05,c(-1,0,0,1))

## initializing the B vector
B <- runif(p,0,0.5)
## initializing the variance of the noise
sigmaEps <- runif(p,0.1,0.8)
## initializing the process Xt
X0 <- B + rnorm(p,0,sigmaEps*10)
## number of time points
n <- 20

## the AR(1) times series process
Xn <- SimulGeneExpressionAR1(MyNet$A,B,X0,sigmaEps,n)
}
