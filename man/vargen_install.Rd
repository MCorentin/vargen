% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vargen.R
\name{vargen_install}
\alias{vargen_install}
\title{Download the files needed to run \code{\link{vargen_pipeline}}}
\usage{
vargen_install(
  install_dir = "./",
  gtex_version = "v8",
  timeout = 10000,
  verbose = FALSE
)
}
\arguments{
\item{install_dir}{the path to the installation folder (default = "./")}

\item{gtex_version}{the version of gtex to download, only "v7" and "v8" are
supported (default = "v8")}

\item{timeout}{the timeout set in options(), reading/downloading files online
might fail with the default timeout of 60 seconds.}

\item{verbose}{if TRUE will print progress messages (default = FALSE)}
}
\value{
nothing, download files in "install_dir".
}
\description{
Only need to run it once. Will download the following files in
"install_dir":
\itemize{
  \item the latest gwas catalog, eg: gwas_catalog_v1.0.2-associations_e93_r2019-01-11.tsv
  \item hg19ToHg38.over.chain.gz (will be unzipped)
  \item GTEx_Analysis_v8_eQTL.tar.gz (will be untared)
  \item GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz
  \item enhancer_tss_associations.bed
}
}
\examples{
vargen_install("./", verbose = TRUE)
}
