\name{data sets}

\alias{CiData}
\alias{CiGene}
\alias{DlData}
\alias{DlGene}

\docType{data}

\title{
Example Data Sets
}

\description{
Examples of the types of data to which \code{\link{mixmod}} and \code{\link{mdmixmod}} may be used to fit mixture models.
}

\format{
For \code{CiData} and \code{DlData}, a list with three elements, \code{binding}, \code{expression}, and \code{conservation}, which are numeric vectors, matrices, or data frames.  For \code{CiGene} and \code{DlGene}, a data frame with elements containing the symbol, name, CG ID, FlyBase ID, chromosome, strand, start position, stop position, and target status for each gene in the corresponding \code{*Data} data set.
}

\details{
Both \code{Ci*} and \code{Dl*} contain data relating to identification of transcription factor (TF) target genes in Drosophila involved in embryonic development.  \code{Ci} relates to cubitus interruptus, a TF involved in regulation of almost all Hedgehog-responsive (Hh-responsive) genes (Von Ohlen et al., 1997).  Binding data represents log-ratios of Ci binding in the regulatory regions of genes vs. background binding.  Expression values are a matrix of log-ratios of expression in mutant vs. wild-type embryos, mutants being homozygous null for one of four proteins known to affect Ci's regulatory function.  The proteins are Smoothened (Smo), Patched (Ptc), and Ci and Hh themselves.  Data are preprocessed and scaled from the raw data available at the Gene Expression Omnibus (GEO) accession number GSE24055.

\code{Dl} relates to the Dorsal TF, which controls dorsal-ventral patterning in early embryogenesis.  Binding data represent the log-ratios of binding in regulatory regions vs. background for dorsal and Snail (Sna), a TF which is an early target of Dl and plays an important role in the dorsal-ventral patterning process.  Raw data are available at GEO GSE26285.  Expression data represent log-ratios of gene expression for different mutant strains with varying levels of Dl throughout the embryo (pipe-/pipe- vs. toll10B and pipe-/pipe- vs. tollrm9/tollrm10).  Raw data are available at GEO GSE5434.

For both \code{Ci*} and \code{Dl*}, cross-species gene sequence conservation is calculated from PhastCons using 12 fly species with one species each of mosquito, honeybee, and beetle as outgroups.  The conservation values used in the analysis are a univariate vector calculated from the sums of PhastCons highly conserved element (HCE) scores for HCEs which overlap genes.  These scores are available from the University of California, Santa Cruz (UCSC) Genome Browser.  \dQuote{Known target} status, represented by the \code{target} element of \code{*Gene}, is calculated from previous studies and from Gene Ontology (GO) and Berkeley Drosophila Genome Project (BDGP) annotation.
}

\examples{\dontrun{
data(CiGene)
data(CiData)
fit <- mdmixmod(CiData, c(2,3,2), topology="chained",
    family=c("pvii", "norm", "pvii"))
rocauc(fit, CiGene$target, quasi=TRUE) # 0.8888273
}}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}, \code{\link{rocauc}} for the functions used in Examples; \code{\link{simulation}} for functions to simulate data with similar characteristics to the real data.
}

\keyword{datasets}