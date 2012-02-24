\name{arth800line}
\alias{arth800line}
\title{Arabidopsis Thaliana temporal gene expression data}

\description{This data set describes the temporal log2 transformed
expression of 800 genes of A. thaliana during the diurnal cycle. The
data are in line, that is 2 repeated measurements time series are
displayed one after the other, separated by a 'NA' value. The 800 genes
are a subset of the data presented in Smith et al. (2004) selected for
periodicity according to the method implemented in the R package
GeneCycle (\url{http://strimmerlab.org/software/genecycle/>}).  }

\usage{data(arth800line)
}

\format{matrix with 800 columns (=genes) and 23 rows (rows 1 to 11 contain the
first measumement time series, row 12 contain 'NA' values and rows 13 to 23 contain the second
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

\author{Lebre Sophie (\url{http://newlsiit.u-strasbg.fr/bfo_en/index.php/Sophie_Lebre}), 

Chiquet Julien (\url{http://stat.genopole.cnrs.fr/~jchiquet/}).}

\keyword{datasets}

\examples{
## load G1DBN library
library(G1DBN)

## load data set
data(arth800line)
id<-c(60, 141, 260, 333, 365, 424, 441, 512, 521, 578, 789, 799)

## plot first ten time series
plot(1:23,arth800line[,60],type="l",ylim=c(2,12), xlab="Time",
 ylab="Log2 transformed expression",lwd=2,
 main="Log2 transformed expression of a subset of genes of A. Thaliana")

color=1
for (i in id){
  color=color+1
  lines(1:23,arth800line[,i,],col=color,lwd=2)
}

}

 
