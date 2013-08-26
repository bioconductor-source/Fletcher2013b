\name{Fletcher2013b.pipelines}
\alias{Fletcher2013pipeline.mra1st}
\alias{Fletcher2013pipeline.mra2nd}
\alias{Fletcher2013pipeline.mraNormals}
\alias{Fletcher2013pipeline.mraTALL}
\alias{Fletcher2013pipeline.masters}
\alias{Fletcher2013pipeline.siRNA}
\alias{Fletcher2013pipeline.synergyShadow}
\alias{Fletcher2013pipeline.consensusnet}
\alias{Fletcher2013pipeline.enrichmap}
\alias{Fletcher2013pipeline.chipseq}
\alias{Fletcher2013pipeline.motifs}
\alias{Fletcher2013b}

\title{
A pipeline to reproduce results for Fletcher et al. 2013.
}

\description{
Pipeline functions to reproduce results for Fletcher et al. 2013.
}

\usage{
 Fletcher2013pipeline.mra1st(hits, minRegulonSize=20, idtype="probeid", 
    pAdjustMethod="holm",tnet="dpi", eps=0, pValueCutoff=1e-4, verbose=TRUE, ...)
 Fletcher2013pipeline.mra2nd(hits, minRegulonSize=20, idtype="probeid", 
    pAdjustMethod="holm",tnet="dpi", eps=0, pValueCutoff=1e-4, verbose=TRUE, ...)
 Fletcher2013pipeline.mraNormals(hits, minRegulonSize=20, idtype="probeid", 
    pAdjustMethod="holm",tnet="dpi", eps=0, pValueCutoff=1e-4, verbose=TRUE, ...)
 Fletcher2013pipeline.mraTALL(hits, minRegulonSize=20, idtype="probeid", 
    pAdjustMethod="holm", tnet="dpi", eps=0, pValueCutoff=0.01, verbose=TRUE, ...)
 Fletcher2013pipeline.siRNA(what="siESR1")
 Fletcher2013pipeline.chipseq(what="SPDEF")
 Fletcher2013pipeline.masters()
 Fletcher2013pipeline.synergyShadow()
 Fletcher2013pipeline.consensusnet()
 Fletcher2013pipeline.enrichmap()
 Fletcher2013pipeline.motifs()
}

\arguments{
\item{hits}{
a character vector of gene identifiers for those considered as hits (see \code{\link[RTN:TNA-class]{TNA-class}}).
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of elements in a regulon that must map to elements of the gene universe 
(see \code{\link[RTN:tna.mra]{tna.mra}}).
}
\item{idtype}{
a single character value specifying the inpute gene id (Options: 'probeid' or 'entrez').
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be used (see \code{\link[stats:p.adjust]{p.adjust}} for details).
}
\item{tnet}{
a single character value specifying which transcriptional network should to used to compute the MRA analysis. Options: "dpi" and "ref". 
}
\item{eps}{
a single numeric value specifying the threshold under which Aracne algorithm should apply the dpi filter 
(see \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}}).
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered significant. 
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE). 
}
\item{what}{
a single character value specifying which data should be used in the pipeline. Options for \code{Fletcher2013pipeline.siRNA}: 'siESR1','siSPDEF', and 'siPTTG1'. Options for \code{Fletcher2013pipeline.chipseq}: 'ESR1','SPDEF', 'FOXA1', and 'GATA3'.
}
\item{...}{
other arguments passed to the \code{RTN} package.
}
}

\value{ 
All results will be saved in the current work directory.
}

\source{
Michael NC Fletcher, Mauro AA Castro, Suet-Feung Chin, Oscar Rueda, Xin Wang, Carlos Caldas, Bruce AJ Ponder, Florian Markowetz, Kerstin B Meyer. Master regulators of FGFR2 signalling and breast cancer risk. Journal paper (submitted). 2012.
}

\author{
Mauro Castro \email{mauro.a.castro@gmail.com}
}

\examples{
\dontrun{
hits <- Fletcher2013pipeline.deg(what="Exp1")
mra1 <- Fletcher2013pipeline.mra1st(hits=hits$E2FGF10)
}
}












