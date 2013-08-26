\name{siRNA}
\docType{data}
\alias{siRNA}

\title{Dataset from siRNA experiments used to reproduce results in Fletcher et al., 2012.}

\usage{
data(siRNA)
}

\description{
The data consists of differentially expressed genes in MCF-7 cells after knockdown experiments.
}

\format{
  A \code{list} object:
\itemize{
\item \code{siRNA$ESR1}: differentially expressed genes in MCF-7 cells after knocking down ESR1 gene. 
\item \code{siRNA$SPDEF}: differentially expressed genes in MCF-7 cells after knocking down SPDEF gene.
\item \code{siRNA$PTTG1}: differentially expressed genes in MCF-7 cells after knocking down PTTG1 gene.
}
}

\source{
Michael NC Fletcher, Mauro AA Castro, Suet-Feung Chin, Oscar Rueda, Xin Wang, Carlos Caldas, Bruce AJ Ponder, Florian Markowetz, Kerstin B Meyer. Master regulators of FGFR2 signalling and breast cancer risk. Journal paper (submitted). 2012.
}

\note{
The differential expression analysis is documented in the package 'Fletcher2013a', and row gene expression data is available at \code{\link[Fletcher2013a:siOTHERS]{siOTHERS}} and \code{\link[Fletcher2013a:siOTHERS]{siESR1}}.
}

\examples{
  data(siRNA)
}
\keyword{datasets}

