\name{rtni.data}
\alias{rtni1st}
\alias{rtni2nd}
\alias{rtniNormals}
\alias{rtniTALL}
\alias{rtniIDs}
\docType{data}

\title{Transcriptional network datasets.}

\usage{
data(rtni1st)
data(rtni2nd)
data(rtniNormals)
data(rtniTALL)
}

\description{
The datasets consist of a transcriptional networks computed by the package RTN.
}

\format{
  A set of \code{TNI} objects:
\itemize{
\item \code{rtni1st}: A TF-centric network based on 2000 breast cancer gene expression profiles - Cohort I (Curtis, C. et al).
\item \code{rtni2nd}: A TF-centric network based on 2000 breast cancer gene expression profiles - Cohort II (Curtis, C. et al).
\item \code{rtniNormals}: A TF-centric network based on normal breast gene expression profiles (Curtis, C. et al).
\item \code{rtniTALL}: A TF-centric network based non-breast cancer gene expression profiles, derived from T-cell acute l
ymphoblastic leukaemia (Van Vlierberghe, P. et al.).
\item \code{rtniIDs}: A \code{data.frame} with gene ids.
}
}

\source{
Michael NC Fletcher, Mauro AA Castro, Suet-Feung Chin, Oscar Rueda, Xin Wang, Carlos Caldas, Bruce AJ Ponder, Florian Markowetz, Kerstin B Meyer. Master regulators of FGFR2 signalling and breast cancer risk. Nature Communications, 4:2464, 2013.

Curtis, C. et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. Nature,  486:346-52, 2012.

Van Vlierberghe, P. et al. ETV6 mutations in early immature human T cell leukemias. J Exp Med, 208:2571-9, 2011.

}


\examples{
  data(rtni1st)
}

\keyword{datasets}
