# ---- VarPhen Pipeline ----

#' @title From keywords, get a list of phenotypes to use with snp mart
#' @description The list of phenotypes obtained can be used with
#' \code{\link{get_variants_from_phenotypes}} to get a list of variants associated
#' with the phenotypes.
#'
#' @param keywords a keyword, all the phenotypes contaning this keyword will be
#' returned (based on \code{\link[base]{grep}} with ignore.case)
#' @param snp_mart optional, a connection to ensembl snp mart, can be created
#' using \code{\link{connect_to_snp_ensembl}} (If missing this function will be
#' used to create the connection).
#' @return a vector of phenotype names
#'
#' @examples
#' # First, connect to snp ensembl
#' snp_mart <- connect_to_snp_ensembl()
#'
#' # Get the phenotypes related to diabetes
#' get_phenotype_terms(keywords = "diabetes", snp_mart = snp_mart)
#' @export
get_phenotype_terms <- function(keywords, snp_mart) {
  if(missing(keywords)) keywords <- ""

  if(missing(snp_mart) || class(snp_mart) != 'Mart'){
    snp_mart <- connect_to_snp_ensembl()
    warning("Snp mart not provided (or not a valid Mart object).",
            "We used one from connect_to_snp_ensembl() instead.")
  }

  # Get all the possible values for the phenotype filter:
  biomart_filters <- biomaRt::listFilterOptions(filter = "phenotype_description",
                                                mart = snp_mart)
  # Unlist the result and use grep to subset the relevant phenotypes
  filters <- unlist(strsplit(x = as.character(biomart_filters), split = ","))
  # paste() is done in case there are more than one keyword, the "|" will act as an OR
  return(filters[grep(paste(keywords,collapse="|"), filters, ignore.case = TRUE)])
}


#' @title Get variants associated with certain phenotypes
#' @description The vector of phenotypes can be obtained with
#' \code{\link{get_phenotype_terms}}
#'
#' @param phenotypes a vector of phenotypes of interest (obtained from
#' \code{\link{get_phenotype_terms}})
#' @param snp_mart a connection to ensembl snp mart, can be generated
#' from \code{\link{connect_to_snp_ensembl}} (which will be used if this argument
#' is not provided)
#' @return a data.frame of variants associated with the phenotypes:
#' \itemize{
#'   \item chr_name (chromosome)
#'   \item chrom_start (variant position)
#'   \item refsnp_id (variant id)
#'   \item refsnp_source (the source of information (eg: "dbSNP"))
#'   \item ensembl_transcript_chrom_strand (strand information (1, -1 or NA))
#'   \item associated_variant_risk_allele (risk allele (eg: "A"))
#'   \item phenotype_description (description of the phenotype (eg: "Diabetes mellitus type 1"))
#'   \item clinical_significance (from clinVar (eg: "pathogenic", "begnin" ...))
#'   \item validated (Variant supporting evidence (1000Genomes, TOPMed, etc...))
#' }
#'
#' @examples
#' # First, connect to ensembl snp mart
#' snp_mart <- connect_to_snp_ensembl()
#'
#' # Get the phenotypes in biomaRt containing the word "Diabetes"
#' DM_phen <- get_phenotype_terms(keywords = "diabetes", snp_mart = snp_mart)
#'
#' # Get all the variants associated with the penotypes in "DM_phen"
#' get_variants_from_phenotypes(phenotypes = DM_phen, snp_mart = snp_mart)
#' @export
get_variants_from_phenotypes <- function(phenotypes, snp_mart) {
  if( missing(snp_mart) || class(snp_mart) != 'Mart'){
    snp_mart <- connect_to_snp_ensembl()
    warning("Snp mart not provided (or not a valid Mart object).",
            "We used one from connect_to_snp_ensembl() instead.")
  }

  pheno_variants <- biomaRt::getBM(attributes = c("chr_name", "chrom_start",
                                                  "refsnp_id", "refsnp_source",
                                                  "associated_variant_risk_allele",
                                                  "phenotype_description"),
                                   filters = c("phenotype_description"),
                                   values = c(phenotypes),
                                   mart = snp_mart
  )
  pheno_variants$chr_name <- format_chr(pheno_variants$chr_name)
  return(pheno_variants)
}
