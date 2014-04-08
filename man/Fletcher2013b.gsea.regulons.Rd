\name{Fletcher2013b.gsea.regulons}
\alias{Fletcher2013gsea.regulons}

\title{
Supplementary GSEA analyses to reproduce results for Fletcher et al. 2013.
}

\description{
Supplementary functions to reproduce results for Fletcher et al. 2013.
}

\usage{
 Fletcher2013gsea.regulons(what = "Exp1", timepoint = 6, verbose = TRUE)
}

\arguments{
\item{what}{
a single character value specifying one of the experimental systems implemented in the study (Options: 'Exp1','Exp2' and 'Exp3').
}
\item{timepoint}{
a single integer value specifying a timepoint of the experimental system (Options: 6 or 24).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE).
}
}

\value{ 
All results will be saved in the current work directory.
}

\source{
Michael NC Fletcher, Mauro AA Castro, Suet-Feung Chin, Oscar Rueda, Xin Wang, Carlos Caldas, Bruce AJ Ponder, Florian Markowetz, Kerstin B Meyer. Master regulators of FGFR2 signalling and breast cancer risk. Nature Communications, 4:2464, 2013.
}

\author{
Mauro Castro \email{mauro.a.castro@gmail.com}
}

\examples{
\dontrun{
Fletcher2013gsea.regulons(what="Exp1")
}
}









