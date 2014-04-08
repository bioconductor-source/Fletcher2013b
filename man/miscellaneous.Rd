\name{miscellaneous}
\alias{risksites}
\alias{randsites}
\alias{chromlen}
\alias{ESR1bdsites}
\alias{FOXA1bdsites}
\alias{GATA3bdsites}
\alias{SPDEFbdsites}
\alias{fimoESR1}
\alias{fimoFOXA1}
\alias{fimoGATA3}
\alias{metaPCNA}
\alias{consensus}
\alias{tfs}
\docType{data}
\alias{miscellaneous}

\title{Miscellaneous datasets.}

\usage{
data(miscellaneous)
}

\description{
Different data sets used to produce a variety of analyses and figures in Fletcher et al., 2013.
}

\format{
  A set of \code{miscellaneous} data objects:
\itemize{
\item \code{risksites}: a data.frame with top 1385 risk SNPs derived from UK2 GWAS study for breast cancer (mapped to genome assembly NCBI36/hg18).
\item \code{randsites}: a data.frame with random SNPs derived from Affy SNP-6 array (sites mapped to hg19).
\item \code{chromlen}: a vector listing human chromosome length (genome assembly NCBI36/hg18).
\item \code{ESR1bdsites}: a data.frame listing ChIP-seq ESR1 binding sites in MCF-7 cells (mapped to genome assembly NCBI36/hg18).
\item \code{FOXA1bdsites}: a data.frame listing ChIP-seq FOXA1 binding sites in MCF-7 cells (mapped to genome assembly NCBI36/hg18).
\item \code{GATA3bdsites}: a data.frame listing ChIP-seq GATA3 binding sites in MCF-7 cells (mapped to genome assembly NCBI36/hg18).
\item \code{SPDEFbdsites}: a data.frame listing ChIP-seq SPDEF binding sites in MCF-7 cells (mapped to genome assembly NCBI36/hg18).
\item \code{fimoESR1}: a list with ESR1 motifs mapped across the human genome.TRANSFAC PWM was used as input for the FIMO DNA motif identification tool, Grant et al., 2011 (mapped to hg19).
\item \code{fimoFOXA1}: a list with FOXA1 motifs mapped across the human genome.TRANSFAC PWM was used as input for the FIMO DNA motif identification tool, Grant et al., 2011 (mapped to hg19).
\item \code{fimoGATA3}: a list with GATA3 motifs mapped across the human genome.TRANSFAC PWM was used as input for the FIMO DNA motif identification tool, Grant et al., 2011 (mapped to hg19).
\item \code{metaPCNA}: a vector listing genes from the metaPCNA proliferation-based gene signature (Venet, D. et al., 2011).
\item \code{consensus}: a list with consensus breast cancer master regulators described in Fletcher et a., 2013.
\item \code{tfs}: a vector listing the transcription factors used to compute the transcriptional networks \code{\link[Fletcher2013b:rtni1st]{rtni1st}},\code{\link[Fletcher2013b:rtni2nd]{rtni2nd}}, \code{\link[Fletcher2013b:rtniNormals]{rtniNormals}} and \code{\link[Fletcher2013b:rtniTALL]{rtniTALL}}.
}
}

\details{
ChIP-seq datasets are represetative of 3 independent experiments, with peaks overlapping in at least 2 out of 3 replicates (taking one as reference). All peaks are provided related to the summit positions (+- 35 bp), including peak height and significance (in the form of -10*log10(pvalue)). Additional details about this and the other datasets are provided in the vignette.
}

\source{

Michael NC Fletcher, Mauro AA Castro, Suet-Feung Chin, Oscar Rueda, Xin Wang, Carlos Caldas, Bruce AJ Ponder, Florian Markowetz, Kerstin B Meyer. Master regulators of FGFR2 signalling and breast cancer risk. Nature Communications, 4:2464, 2013.

Grant CE, Bailey TL, Noble WS: FIMO: scanning for occurrences of a given motif. Bioinformatics, 27(7):1017-1018, 2011.

Venet, D., Dumont, J.E. & Detours, V. Most random gene expression signatures are significantly associated with breast cancer outcome. PLoS Comput Biol, 7:e1002240, 2011.
}

\examples{
  data(miscellaneous)
}

\keyword{datasets}

