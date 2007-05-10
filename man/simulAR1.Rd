\name{simulAR1}
\alias{simulAR1}
\title{First order multivariate Auto-Regressive process generation}
\description{This function generates multivariate time series according to the following first order Auto-Regressive process,

X(t)= A X(t-1) + B + e(t),

where matrix A has size p x p and arrays X(t), B and e(t) have length p.
e(t) follows a zero-centered multivariate gaussian distribution whose variance matrix S is diagonal.
First, matrix A, array B and diagonal of S are randomly generated.
Each diagonal term S[i,i] is uniformly generated from U([minSig,maxSig]).
The elements of matrix A and array B are uniformly generated from U([-maxA,-minA],[minA,maxA])
and U([minB,maxB]) respectively.
Second, the time series data are generated according the so defined AR(1) model.
}

\usage{out<-simulAR1(p,n,edgeProp,minA,maxA,minB,maxB,minSig,maxSig)
}

\arguments{\item{p}{the desired dimension of the multivariate time series.}
	\item{n}{the desired length of the time serie.}
	\item{edgeProp}{the desired proportion of non zero coefficient in the AR transition matrix.}
	\item{minA}{the minimum value for matrix A elements generation.}
	\item{maxA}{the maximum value for matrix A elements generation.}
	\item{minB}{the minimum value for matrix B elements generation.}
	\item{maxB}{the maximum value for matrix B elements generation.}
	\item{minSig}{the minimum value for the diagonal of covariance matrix S generation.}
	\item{maxSig}{the maximum value for the diagonal of covariance matrix S generation.}
}

\value{A list with out\$data a matrix, with n rows (=length) and p columns (=dimension), containing the generated time series,
out\$A the AR generated matrix A (p x p),
out\$B the AR generated vector B (p),
out\$sig the generated diagonal (p) of covariance matrix S.
}

\author{Lebre Sophie (\url{http://stat.genopole.cnrs.fr/~slebre}).}

\keyword{datagen}
\keyword{ts}

\examples{#generate AR(1) time series 
AR<-simulAR1(p=10,n=50,edgeProp=0.02,minA=0.5,maxA=1.5,minB=0,maxB=1,minSig=0.1,maxSig=0.8)}
