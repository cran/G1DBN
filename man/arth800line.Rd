\name{arth800line}
\alias{arth800line}
\title{Arabidopsis Thaliana temporal gene expression data}

\description{This data set describes the temporal expression of 800 genes of
   A. thaliana during the diurnal cycle. The data are in line, that is 2
   repeated measurements time series are displayed one after the other, separated by a NA value. The 800 genes are a
   subset of the data presented in Smith et al. (2004) selected for periodicity according to the method implemented in
the R package GeneCycle (\url{http://strimmerlab.org/software/genecycle/>}).
}

\usage{data(arth800line)
}

\format{matrix with 800 columns (=genes) and 23 rows (rows 1 to 11 for the
first measumement time series, row 12 is NA and rows 13 to 23 for the second
experiment time series).} 

\source{The microarray experiments were performed in the laboratory of S.
   Smith (Edinburgh). The data are available from the NASCArrays
   database (\url{http://affymetrix.arabidopsis.info/} under 
   experiment reference number NASCARRAYS-60.
}



\references{Smith et al. 2004. Diurnal changes in the transcriptom encoding 
   enzymes of starch metabolism provide evidence for both
   transcriptional and posttranscriptional regulation of starch
   metabolism in Arabidopsis leaves. Plant Physiol. 136: 2687-2699.}

\author{Lebre Sophie (\url{http://stat.genopole.cnrs.fr/~slebre}).}

\keyword{datasets}

\examples{
   # load G1DBN library
   library(G1DBN)

   # load data set
   data(arth800line)

   # plot first ten time series
   plot(1:23,arth800line[,1],type="l",ylim=c(5,12))
   for (i in 2:10)lines(1:23,arth800line[,i,],col=i)

}

 
