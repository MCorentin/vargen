# --- Utils ----

#' @title Connect to gene Mart
#' @description Connect to the "hsapiens_gene_ensembl" dataset in the "ENSEMBL_MART_ENSEMBL"
#' BioMart database, using \code{\link[biomaRt]{useEnsembl}}
#'
#' @param mirror (optional) to use an alternative mirror. see
#' \code{\link[biomaRt]{useEnsembl}} for the possible values (default:"www")
#' @return a Mart object with a connection to the "hsapiens_gene_ensembl" dataset
#'
#' @examples
#' # Connect with the default mirror
#' gene_mart <- connect_to_gene_ensembl()
#'
#' # Connect with the "www" mirror
#' gene_mart <- connect_to_gene_ensembl(mirror="www")
#' @export
connect_to_gene_ensembl <- function(mirror = "www"){
  gene_mart <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                                   host = "https://www.ensembl.org",
                                   mirror = mirror,
                                   dataset = "hsapiens_gene_ensembl")
  return(gene_mart)
}


#' @title Connect to snp Mart
#' @description Connect to the "hsapiens_snp" dataset in the "snp" BioMart
#' database, using \code{\link[biomaRt]{useEnsembl}}
#'
#' @param mirror (optional) to use an alternative mirror. see
#' \code{\link[biomaRt]{useEnsembl}} for the possible values (default:"www")
#' @return a Mart object with a connection to the "hsapiens_snp" dataset
#'
#' @examples
#' # Connect with the default mirror
#' ensembl <- connect_to_snp_ensembl()
#'
#' # Connect with the "www" mirror
#' www_ensembl <- connect_to_snp_ensembl(mirror="www")
#' @export
connect_to_snp_ensembl <- function(mirror = "www"){
  snp_mart <- biomaRt::useEnsembl(biomart = "snp",
                                  host = "https://www.ensembl.org",
                                  mirror = mirror,
                                  dataset = "hsapiens_snp")
  return(snp_mart)
}


#' @title Using ensembl API to get the variants at certain locations
#' @description Use \code{\link[httr]{GET}} to query the ensembl REST API to
#' gather variants located on specific locations. The content is fetched as a
#' JSON and then translated as a data.frame, using \code{\link[jsonlite]{toJSON}}
#' /!\ the column order is sometimes different.
#'
#' @param locations a vector of locations (eg: c("chr1:1000:10000", "chr2:10-100"))
#' @param verbose if TRUE, will print progress messages (default: FALSE)
#' @return a data.frame containing the following columns:
#' \itemize{
#'   \item clinical_significance
#'   \item start
#'   \item consequence_type
#'   \item alleles
#'   \item seq_region_name
#'   \item feature_type
#'   \item source
#'   \item strand
#'   \item assembly_name
#'   \item end
#'   \item id
#' }
#'
#' @examples
#' get_variants_from_locations(c("chr1:100:1000", "chr10:2004:12042"))
#' @export
get_variants_from_locations <- function(locations, verbose = FALSE) {
  # We use: https://rest.ensembl.org/
  variants_at_loc <- data.frame()

  # On top of "x-ratelimit-remaining", we also try to do less than 14 request
  # per seconds. higher rate seems to throw a 429 error (undocumented rate limit?)
  reqs_per_sec_limit <- 14
  n_reqs <- 1

  for(curr_loc in locations){

    # The 14 request per seconds limit
    if(n_reqs >= reqs_per_sec_limit){
      Sys.sleep(1)
      n_reqs <- 1
    }

    query <- paste0("/overlap/region/human/",curr_loc,
                    "?content-type=application/json;feature=variation")

    get_output <- httr::GET(paste("https://rest.ensembl.org", query, sep = ""),
                            httr::content_type("application/json"))

    # We get the header to check for query limits
    header <- httr::headers(get_output)

    while(get_output$status_code == "429"){
      print(paste0("Hit 429 error, waiting for ", header$`Retry-After`, " seconds"))
      Sys.sleep(header$`Retry-After`)

      get_output <- httr::GET(paste("https://rest.ensembl.org", query, sep = ""),
                              httr::content_type("application/json"))

      # We get the header to check for query limits
      header <- httr::headers(get_output)
    }

    # Transform status to R error
    httr::stop_for_status(get_output)

    # If the number of allowed queries is lower than 2:
    if(header$`x-ratelimit-remaining` < 2){
      # We wait until the rate limit is reset (in seconds)
      if(verbose) print(paste0("Waiting ", header$`x-ratelimit-reset`,
                               " seconds for ensembl rate limit to reset"))

      Sys.sleep(header$`x-ratelimit-reset`)
    }

    variants_at_loc <- rbind(variants_at_loc, jsonlite::fromJSON(
      jsonlite::toJSON(httr::content(get_output))))

    n_reqs <- n_reqs + 1
  }

  return(variants_at_loc)
}


#' @title Format chromosome names
#' @description Transform 1 to 22, X, Y and MT into chr1 to chr22, chrX, chrY and chMT.
#' Regex needed because some chr are named CHR17_... and we do not want chrCHR17
#'
#' @param chr a vector of chromosome names
#' @return the formatted vector of chromosome names (eg: "chr1" instead of "1")
format_chr <- function(chr){
  # This should match 1 to 22 and X, Y, MT. We need to replace them as
  # chr1, chrX etc....
  chr_to_rename <- grep(pattern = "^\\d{1,2}$|^X$|^Y$|^MT$", chr,
                        ignore.case = TRUE)

  chr[chr_to_rename] <- paste0("chr", chr[chr_to_rename])

  return(chr)
}


#' @title Format the variants data.frame
#' @description Used to have consistent output for the variants from OMIM, FANTOM5,
#' GTEx and GWAS.
#'
#' @param chr list of variant chromosome
#' @param pos list of variant positions
#' @param rsid list of variant rsid
#' @param ensembl_gene_id gene id associated to the variants
#' @param hgnc_symbol hgnc_symbol associated to the variants
#'
#' @return a data.frame of variants with the following columns
#' \itemize{
#'   \item chr (chromosome)
#'   \item pos (position of the variant)
#'   \item rsid (variant ID)
#'   \item ensembl_gene_id ("gene id" of the gene associated with the variants)
#'   \item hgnc_symbol ("hgnc symbol" of the gene associated with the variants)
#' }
format_output <- function(chr, pos, rsid, ensembl_gene_id, hgnc_symbol) {
  chr <- format_chr(chr)

  return(data.frame(chr = chr,
                    pos = pos,
                    rsid = rsid,
                    ensembl_gene_id = ensembl_gene_id,
                    hgnc_symbol = hgnc_symbol, stringsAsFactors = FALSE))
}


#' @title Annotate variants from a vector of variant IDs
#' @description use \code{\link[myvariant]{getVariants}} and process the output
#' to return a data.frame instead of nested lists.
#'
#' @param rsid a vector of variant IDs
#' @param verbose if TRUE will print progress messages (default = FALSE)
#'
#' @return a data.frame with the following columns:
#' \itemize{
#'   \item rsid (variant id)
#'   \item cadd phred (phred score from CADD, higher more deleterious)
#'   \item fathmm_xf_score (from 0 to 1, higher is more deleterious)
#'   \item fathmm_xf_pred ("D"(DAMAGING) if score > 0.5, "N"(NEUTRAL) otherwise)
#'   \item annot_type (the annotation type from CADD, eg: "Intergenic")
#'   \item consequence (the consequence from CADD, eg: "DOWNSTREAM")
#'   \item clinical_significance (risk factor from clinvar, eg: "Benign", "risk factor")
#'   \item snpeff_ann (impact annotation from snpeff eg: "MODIFIER")
#' }
#'
#' @examples
#' annotate_variants(rsid = c("rs1225680362", "rs12395043", "rs746318172"))
#' @export
annotate_variants <- function(rsid, verbose = FALSE) {

  myvariant <- myvariant::MyVariant()
  myvariant@step  <- 700

  # First, use "getVariants" to annotate the snps
  rsid_annotated <- myvariant::getVariants(hgvsids = rsid, verbose = verbose,
                                           myvariant = myvariant,
                                           fields =  c("cadd", "dbnsfp", "clinvar", "snpeff", "vcf"))

  # Then, format the output as a data.frame
  # Checking "is.null" is needed to avoid errors when there is a list of variants
  # without some features (eg: no cadd.phred scores). Can happen for small lists.
  if(is.null(rsid_annotated$cadd.phred)){
    cadd_phred <- NA
  } else {
    cadd_phred <- unlist(rsid_annotated$cadd.phred)
  }

  if(is.null(rsid_annotated$dbnsfp.fathmm.xf.coding_score)){
    fathmm_score <- NA
  } else {
    fathmm_score <- rsid_annotated$dbnsfp.fathmm.xf.coding_score
  }

  if(is.null(rsid_annotated$dbnsfp.fathmm.xf.coding_pred)){
    fathmm_pred <- NA
  } else {
    fathmm_pred <- rsid_annotated$dbnsfp.fathmm.xf.coding_pred
  }

  if(is.null(rsid_annotated$cadd.annotype)){
    annot_type <- NA
  } else {
    annot_type <- sapply(rsid_annotated$cadd.annotype,
                         paste, collapse = ";", USE.NAMES = FALSE)
    annot_type <- gsub("\n", "", annot_type)
  }

  if(is.null(rsid_annotated$cadd.consequence)){
    consequence <- NA
  }  else {
    consequence <- sapply(rsid_annotated$cadd.consequence,
                          paste, collapse = ";", USE.NAMES = FALSE)
    consequence <- gsub("\n", "", consequence)
  }

  if(is.null(rsid_annotated$snpeff.ann)){
    snpeff <- NA
  }  else {
    snpeff <- sapply(sapply(rsid_annotated$snpeff.ann, "[", "putative_impact"),
                     paste, collapse = ";", USE.NAMES = FALSE)
    snpeff <- gsub("\n", "", snpeff)
  }

  # For this one, getVariants returns a different format depending on the data
  # fetched, so we also need to check if this is a list to extract only relevant information.
  if(is.list(rsid_annotated$clinvar.rcv)) {
    clinical_significance <- sapply(sapply(rsid_annotated$clinvar.rcv, "[", "clinical_significance"),
                                    paste, collapse = ";")
  } else if(is.null(rsid_annotated$clinvar.rcv)) {
    clinical_significance <- NA
  } else {
    clinical_significance <- sapply(rsid_annotated$clinvar.rcv.clinical_significance,
                                    paste, collapse = ";", USE.NAMES = FALSE)
  }

  if(is.null(rsid_annotated$vcf.ref)){
    vcf.ref <- NA
  } else{
    vcf.ref <- rsid_annotated$vcf.ref
  }
  if(is.null(rsid_annotated$vcf.alt)){
    vcf.alt <- NA
  } else {
    vcf.alt <- rsid_annotated$vcf.alt
  }

  rsid_annotated_df <- data.frame(rsid = unlist(rsid_annotated$query),
                                  ref = vcf.ref,
                                  alt = vcf.alt,
                                  cadd_phred = cadd_phred,
                                  fathmm_xf_score = fathmm_score,
                                  fathmm_xf_pred = fathmm_pred,
                                  annot_type = annot_type,
                                  consequence = consequence,
                                  clinical_significance = clinical_significance,
                                  snpeff_ann = snpeff,
                                  stringsAsFactors = FALSE)

  #Replacing the snpeff annotation from "c(\"MODIFIER\", \"MODIFIER\")" to "MODIFIER;MODIFIER"
  rsid_annotated_df$snpeff_ann <- gsub(pattern = "\"|c\\(|\\)", replacement = "",
                                       x = rsid_annotated_df$snpeff_ann)
  rsid_annotated_df$snpeff_ann <- noquote(gsub(pattern = ", ", replacement = ";",
                                               x = rsid_annotated_df$snpeff_ann))

  # Same as above for clinical significance:
  rsid_annotated_df$clinical_significance <- gsub(pattern = "\"|c\\(|\\)", replacement = "",
                                                  x = rsid_annotated_df$clinical_significance)
  rsid_annotated_df$clinical_significance <- noquote(gsub(pattern = ", ", replacement = ";",
                                                          x = rsid_annotated_df$clinical_significance))

  return(unique(rsid_annotated_df))
}
