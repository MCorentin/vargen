% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwas.R
\name{get_gwas_variants}
\alias{get_gwas_variants}
\title{Get the variants from the gwas catalog associated to the traits of interest}
\usage{
get_gwas_variants(gwas_traits, gwas_cat)
}
\arguments{
\item{gwas_traits}{a vector of gwas traits, can be obtained from
\code{\link{list_gwas_traits}}}

\item{gwas_cat}{output from \code{\link{create_gwas}}}
}
\value{
a data.frame contaning the variants linked to the traits in the gwas catalog
The data.frame contains the following columns:
\itemize{
  \item chr (chromosome)
  \item pos (position of the variant)
  \item rsid (variant ID)
  \item ensembl_gene_id ("gene id" of the gene associated with the variants)
  \item hgnc_symbol ("hgnc symbol" of the gene associated with the variants)
  \item source (here "gwas")
}
}
\description{
uses \code{\link[gwascat]{subsetByTraits}} to get the variants.
}
\examples{
# if you need to get the list of gwas traits:
# list_gwas_traits(keywords = c("Obesity"))
obesity_gwas <- c("Obesity (extreme)", "Obesity-related traits", "Obesity")
gwas_cat <- create_gwas()

gwas_variants <- get_gwas_variants(obesity_gwas, gwas_cat)
}
