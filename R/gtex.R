# ---- GTEx ----

#' @title Generates a vector containing the available GTEx tissue files
#' @description The directory can be downloaded with:
#' download.file("http://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
#               destfile = "gwas_catalog_v1.0.2-associations_e93_r2019-01-11.tsv")
#'
#' @param gtex_dir is the directory contaning the files from
#' "GTEx_Analysis_v8_eQTL.tar.gz" (cf: https://gtexportal.org/home/datasets)
#' @return data.frame containing the full file path of each tissue and keywords to select them
#'
#' @examples
#' vargen_install(install_dir = "./vargen_data/")
#' list_gtex_tissues("./vargen_data/GTEx_Analysis_v8_eQTL/")
#' @export
list_gtex_tissues <- function(gtex_dir){
  if(missing(gtex_dir)) stop("Please specify the path to the
                                        GTEx directory (eg: 'GTEx_Analysis_v8_eQTL')")

  file_paths <- list.files(path = paste0(gtex_dir), full.names = TRUE)
  filenames <- list.files(path = paste0(gtex_dir), full.names = FALSE)

  # Selecting only the variant gene pairs files
  var_gene_pairs_paths <- file_paths[grep(".signif_variant_gene_pairs.txt.gz", file_paths)]
  var_gene_pairs_names <- filenames[grep(".signif_variant_gene_pairs.txt.gz", filenames)]

  # We get the tissue name, this will serve as a keyword to select the tissue
  keywords <- sub(pattern = ".v[78].*",  replacement = "", x = var_gene_pairs_names)

  return(data.frame(keywords = keywords, filepaths = var_gene_pairs_paths))
}


#' @title Get the GTEx files corresponding to the tissues of interest
#' @description Use \code{\link{list_gtex_tissues}} to get the list of GTEx files
#' and their associated keywords. Grep the user's query against the list of
#' keywords (case is ignored) to return all the corresponding tissue files.
#' eg: "adipose" will return the files corresponding to "Adipose_Subcutaneous" AND
#' "Adipose_Visceral_Omentum". The tissues in the query are grepped one after the
#' other.
#' If one of the keywords does not correspond to any tissues, \code{\link[base]{stop}}
#' is called and an error message is displayed.
#'
#' @param gtex_dir is the directory contaning the files from
#' "GTEx_Analysis_v8_eQTL.tar.gz" (cf: https://gtexportal.org/home/datasets)
#' @param tissues_query is a vector of keywords to select the tissues of interest
#' (eg: "adipose", "brain", "Heart_Left_Ventricle"). NOT case sensitive. (default:
#' "", return everything)
#' @return the list of files corresponding to the tissues of interest
#'
#' @examples
#' vargen_install(install_dir = "./vargen_data/")
#' select_gtex_tissues("./vargen_data/GTEx_Analysis_v8_eQTL/",
#'                     c("adipose", "Heart_Left_Ventricle", "leg"))
#'
#' select_gtex_tissues("./vargen_data/GTEx_Analysis_v8_eQTL/",
#'                     "brain")
#' @export
select_gtex_tissues <- function(gtex_dir, tissues_query = ""){
  if(missing(gtex_dir)) stop("Please specify the path to the
                              GTEx directory (eg: 'GTEx_Analysis_v8_eQTL')")

  selected_files <- c()
  gtex_tissues <- list_gtex_tissues(gtex_dir)

  for(tissue_keyword in tissues_query){
    corr_tissues <- gtex_tissues[grep(tissue_keyword, gtex_tissues[,1],
                                      ignore.case = TRUE),]

    # If the grep had no rows: tissue does not exists, we exit
    if(dim(corr_tissues)[1] == 0){
      stop(paste0("Tissue '", tissue_keyword,"' not recognized, aborting...",
                  "(please check the keyword and the GTEx folder)"))
    }

    selected_files <- c(selected_files, as.character(corr_tissues[,2]))
  }

  return(selected_files)
}


#' @title Convert GTEx IDs to rsid
#' @description GTEx IDs are in the format "chr_pos_ref_alt_build". To make it
#' consistent with the rest of the package, we are converting them to rsids.
#' We us the "gtex_lookup" which is a data.frame created by reading the gtex lookup file.
#' (see: get_gtex_variants)
#'
#' @param gtex_variants a data.frame of gtex variants obtained during "get_gtex_variants"
#' should contain at least one column with the gtex_id as "chr_pos_alt_ref_build"
#' @param gtex_lookup_file the lookup file, GTEx to rsids. "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"
#' Can be obtained using \code{\link{vargen_install}}.
#' @param verbose if true, will print information about the conversion (default: FALSE)
#' @return a vector of rsids (as some gtex ids do not have a corresponding rsid
#' the output can be smalled than the input)
convert_gtex_to_rsids <- function(gtex_variants, gtex_lookup_file, verbose = FALSE) {

  # Check the gtex build (b37 or b38)
  gtex_build <- unique(sub(".*_", "", gtex_variants$variant_id))
  if(length(gtex_build) > 1){
    stop(paste0("All the GTEx variants are not in the same build, builds detected: ",
                paste(gtex_build, collapse = ", ")))
  }

  if(verbose) print("Loading GTEx lookup table... Please be patient")
  # Column 1, 7 and 8 contains the GTEx ids v8, rsids and GTEx ids v7 respectively
  gtex_lookup <- data.table::fread(select = c(1,7,8), sep = "\t", header = TRUE,
                                   file = gtex_lookup_file, stringsAsFactors = FALSE)

  # In the GTEx lookup table:
  # "variant_id" correspond to the variants from b38
  # "variant_id_b37" correspond to the variants from b37
  # However, from the tissue files, the name is always "variant_id" regardless
  # of the GTEx version
  if(gtex_build == "b38"){
    # Subsetting the lookup table to make merge faster
    gtex_lookup <- gtex_lookup[gtex_lookup$variant_id %in% gtex_variants$variant_id,]

    gtex_variants_rsids <- merge(x = gtex_variants, y = gtex_lookup,
                                 by = "variant_id", all.x = TRUE)
  }

  if(gtex_build == "b37"){
    # Subsetting the lookup table to make merge faster
    gtex_lookup <- gtex_lookup[gtex_lookup$variant_id_b37 %in% gtex_variants$variant_id,]

    gtex_variants_rsids <- merge(x = gtex_variants, y = gtex_lookup,
                                 by.x = "variant_id", by.y = "variant_id_b37",
                                 all.x = TRUE)
  }

  # If verbose on, tell the user how many gtex snps have no corresponding rsids
  if(verbose){
    n_removed <- nrow(gtex_variants_rsids[gtex_variants_rsids$rs_id_dbSNP151_GRCh38p7 == ".",])
    print(paste0("Number of GTEx ids removed (no corresponding rsid): ", n_removed))
  }

  rm(gtex_lookup)

  # We remove gtex ids without a rsid:
  return(gtex_variants_rsids[gtex_variants_rsids$rs_id_dbSNP151_GRCh38p7 != ".",])
}


#' @title Get variants from GTEx linked to the given ensembl genes
#' @description Take as input one or more tissue files from
#'  "GTEx_Analysis_v8_eQTL" (or v7), as well as a vector of ensembl gene ids.
#'  This function will return the variants that are associated with changes
#'  in the expression of the selected genes in the selected tissues. The
#'  assocations are based on the "signif_variant_gene_pairs" files.
#'  The ensembl ids from the GTEx file will be converted to stable ids.
#'
#' @param tissue_files a vector containing the name of the "signif_variant_gene_pairs.txt.gz"
#' files. This will be read using \code{\link[base]{gzfile}}. Output from
#' \code{\link{select_gtex_tissues}} can be used.
#' @param omim_genes output from \code{\link{get_omim_genes}}
#' @param gtex_lookup_file the lookup file, GTEx to rsids. "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"
#' Can be obtained using \code{\link{vargen_install}}.
#' @param snp_mart optional, a connection to ensembl snp mart, can be created
#' using \code{\link{connect_to_snp_ensembl}} (If missing this function will be
#' used to create the connection).
#' @param verbose will be given to subsequent functions to print progress.
#' @return a data.frame with information about the variants associated with a
#' change in expression on the gene of interest.
#' The data.frame will contain the following columns:
#' \itemize{
#'   \item chr (chromosome)
#'   \item pos (position of the variant)
#'   \item rsid (variant ID)
#'   \item ensembl_gene_id ("gene id" of the gene associated with the variant)
#'   \item hgnc_symbol ("hgnc symbol" of the gene associated with the variant)
#'   \item source (here the value will be "gtex")
#' }
get_gtex_variants <- function(tissue_files, omim_genes, gtex_lookup_file,
                              snp_mart, verbose = FALSE){

  gtex_variants_final <- data.frame()

  if(missing(snp_mart) || class(snp_mart) != 'Mart'){
    if(verbose) print("Connecting to the snp mart...")
    snp_mart <- connect_to_snp_ensembl()
    warning("Snp mart not provided (or not a valid Mart object).",
            "We used one from connect_to_snp_ensembl() instead.")
  }

  # First, get the significant variants from all the tissues
  list.variants.tissues <- vector('list', length(tissue_files))
  i <- 1
  for(file in tissue_files){
    tissue_variants <- data.table::fread(select = c(1,2), sep = "\t", header = TRUE,
                                         file = file, stringsAsFactors = FALSE)
    # We get the tissue name from the filename
    tissue <- sub(pattern = ".v[78].*",  replacement = "", x = basename(file))
    tissue_variants$tissue <- paste0("gtex (", tissue, ")")

    # Tranforming the "ensembl ID" from GTEx to "stable ensembl gene id"
    # eg: ENSG00000135100.17 to ENSG00000135100 (the ".17" correspond to the version number)
    tissue_variants$stable_gene_id <- stringr::str_replace(tissue_variants$gene_id,
                                                           pattern = ".[0-9]+$",
                                                           replacement = "")
    list.variants.tissues[[i]] <- tissue_variants
    i <- i + 1
  }
  # Removing the NULL elements of the list if exists and transform the list to a DF
  list.variants.tissues <- list.variants.tissues[!sapply(list.variants.tissues, is.null)]
  list.variants.tissues <- do.call('rbind', list.variants.tissues)


  # Select the variants that are affecting our genes of interest (omim_genes)
  gtex_variants <- list.variants.tissues[list.variants.tissues$stable_gene_id %in%
                                           omim_genes$ensembl_gene_id,]

  if(nrow(gtex_variants) > 0){
    gtex_variants <- convert_gtex_to_rsids(gtex_variants = gtex_variants,
                                           gtex_lookup_file = gtex_lookup_file,
                                           verbose = verbose)
    if(nrow(gtex_variants) > 0){

      colnames(gtex_variants)[which(names(gtex_variants) == "rs_id_dbSNP151_GRCh38p7")] <- "rsid"

      # GTEx use dbsnp v151. Some rsids have been merged. We can use biomart with
      # "synonym" to get the corresponding rsid in the new version:
      rsid_updated <- biomaRt::getBM(attributes = c("synonym_name", "refsnp_id"),
                                     filters = "snp_synonym_filter",
                                     values = gtex_variants$rsid,
                                     mart = snp_mart)

      # If we do not have synonyms to merge we skip the rsid update from db151
      # since they are all valid
      if(nrow(rsid_updated) > 0){
        # This will add a column to gtex_variants with "<NA>" except where there is
        # a synonym SNP, if so we will change the rsid by the "ref_snp" from rsid_updated
        gtex_variants_update <- merge(x = gtex_variants, y = rsid_updated,
                                      by.x = "rsid", by.y = "synonym_name", all.x = TRUE)

        gtex_variants_update[!is.na(gtex_variants_update$refsnp_id),"rsid"] <-
          gtex_variants_update[!is.na(gtex_variants_update$refsnp_id),"refsnp_id"]

      }

      # get rsid positions with biomaRt
      variants_pos <- biomaRt::getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                                     filters = "snp_filter",
                                     values = gtex_variants$rsid,
                                     mart = snp_mart)

      gtex_variants_pos <- merge(x = gtex_variants, y = variants_pos, all.x = TRUE,
                                 by.x = "rsid", by.y = "refsnp_id")

      # To add the hgnc symbol:
      gtex_variants_hgnc <- merge(x = gtex_variants_pos,
                                  y = omim_genes[,c("ensembl_gene_id", "hgnc_symbol")],
                                  by.x = "stable_gene_id",
                                  by.y = "ensembl_gene_id",
                                  all.x = TRUE)

      # Use format output and return the variants
      gtex_variants_formatted <- format_output(chr = gtex_variants_hgnc$chr_name,
                                               pos =  gtex_variants_hgnc$chrom_start,
                                               rsid = gtex_variants_hgnc$rsid,
                                               ensembl_gene_id = gtex_variants_hgnc$stable_gene_id,
                                               hgnc_symbol = gtex_variants_hgnc$hgnc_symbol)

      # Merge back the tissue after the formatting
      gtex_variants_final <- merge(x = gtex_variants_formatted,
                                   y = gtex_variants_hgnc[,c("rsid", "stable_gene_id", "tissue")],
                                   by.x = c("rsid","ensembl_gene_id"),
                                   by.y = c("rsid", "stable_gene_id"),
                                   all.x = TRUE)

      # Reorder and rename the columns to match "master_variants"
      gtex_variants_final <- gtex_variants_final[c("chr", "pos", "rsid", "ensembl_gene_id",
                                                   "hgnc_symbol", "tissue")]
      colnames(gtex_variants_final)[which(names(gtex_variants_final) == "tissue")] <- "source"
    }
  }

  return(unique(gtex_variants_final))
}
