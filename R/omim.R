#---- OMIM ----

#' @title List the available OMIM morbid IDs and descriptions
#' @description List the "OMIM IDs" available in biomaRt, will also
#' print the corresponding description. Use \code{\link[biomaRt]{getBM}} with the
#' following filter: "with_mim_morbid".
#' This function is useful to get OMIM IDs for \code{\link{get_omim_genes}}
#'
#' @param gene_mart a connection to hsapiens gene mart, can be generated
#' from \code{\link{connect_to_gene_ensembl}} (which will be used if this argument
#' if not specified).
#' @param keywords a vector containing a list of keywords. This will be grepped
#' (case insensitive) against the list of omim descriptions to retrieve a list of
#' omim IDs.
#'
#' @return a dataframe with two columns ("mim_morbid_accession" and "mim_morbid_description")
#'
#' @examples
#' gene_mart <- connect_to_gene_ensembl()
#' list_omim_accessions(gene_mart)
#' list_omim_accessions(gene_mart, "alzheimer")
#' list_omim_accessions(gene_mart, c("alzheimer", "obesity"))
#' @export
list_omim_accessions <- function(gene_mart, keywords){
  if(missing(gene_mart) || class(gene_mart) != 'Mart'){
    gene_mart <- connect_to_gene_ensembl()
    warning("Gene mart not provided (or not a valid Mart object).",
            "We used one from connect_to_gene_ensembl() instead.")
  }

  omim_list <- biomaRt::getBM(mart = gene_mart,
                              attributes = c("mim_morbid_accession",
                                             "mim_morbid_description"),
                              filters = "with_mim_morbid",
                              values = c(TRUE),
                              uniqueRows = TRUE)

  if(!missing(keywords)){
    # use "|" as an OR for the grep search
    keys <- paste(keywords, collapse = "|")
    omim_list <- omim_list[grep(keys,
                                omim_list$mim_morbid_description,
                                ignore.case = TRUE),]

    if(nrow(omim_list) == 0) print(paste0("No OMIM ID found for '", keywords, "'"))
  }

  return(omim_list)
}


#' @title Get ensembl gene IDs from "mim_morbid_accession"
#' @description Use \code{\link[biomaRt]{getBM}} to filter results by
#' "mim_morbid_accession". One or more accessions can be given as input.
#'
#' @param omim_ids is a vector of morbity accession
#' eg: c("125853","222100") for Diabetes
#' @param gene_mart a connection to hsapiens gene mart, can be generated
#' from \code{\link{connect_to_gene_ensembl}} (which will be used if this argument
#' if not specified).
#' @return a dataframe with
#' \itemize{
#'  \item ensembl_gene_id (ensembl stable gene ID, eg:"ENSG00000134640")
#'  \item chromosome_name (chromosome name)
#'  \item start_position (gene start)
#'  \item end_position (gene end)
#'  \item hgnc_symbol (gene id from the Hugo Gene Nomenclature Committee, eg:"MTNR1B")
#'  \item mim_morbid_accession (the OMIM morbid ID, eg:"125853")
#'  \item mim_morbid_accession (the OMIM description, eg: "NONINSULIN-DEPENDENT DIABETES MELLITUS")
#' }
#'
#' @examples
#' gene_mart <- connect_to_gene_ensembl()
#' DM_genes <- get_omim_genes(c("125853","222100"), gene_mart)
#' @export
get_omim_genes <- function(omim_ids, gene_mart) {

  if(missing(gene_mart) || class(gene_mart) != 'Mart'){
    gene_mart <- connect_to_gene_ensembl()
    warning("Gene mart not provided (or not a valid Mart object).",
            "We used one from connect_to_gene_ensembl() instead.")
  }

  mim_genes <- data.frame()

  mim_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                             "start_position", "end_position",
                                             "hgnc_symbol", "mim_morbid_accession",
                                             "mim_morbid_description"),
                              filters = c("mim_morbid_accession"),
                              values = c(omim_ids),
                              mart = gene_mart, uniqueRows = TRUE)

  if(nrow(mim_genes) == 0) print(paste0("No genes found for: ", paste(omim_ids, collapse = ", ")))

  if(nrow(mim_genes) > 0) mim_genes$chromosome_name <- format_chr(chr = mim_genes$chromosome_name)

  return(mim_genes)
}


#' @title Get the variants located on the genes
#' @description Get the variants located between the start and stop positions
#' of the genes given as input.
#'
#' @param genes list of genes (can be obtained from \code{\link{get_omim_genes}})
#' @param verbose if true, will print progress information (default: FALSE)
#'
#' @return a data.frame of variants with the following columns (see: \code{\link{format_output}})
#' \itemize{
#'   \item chr (chromosome)
#'   \item pos (position of the variant)
#'   \item rsid (variant ID)
#'   \item ensembl_gene_id ("gene id" of the gene associated with the variant)
#'   \item hgnc_symbol ("hgnc symbol" of the gene associated with the variant)
#'   \item source (here the value will be "omim")
#' }
get_genes_variants <- function(genes, verbose = FALSE){
  # Get the variants on the OMIM genes:
  genes_variants <- data.frame()

  # Used to avoid growing the data frame for every gene, we only do it at the end
  # with: "do.call('rbind', list.variants)"
  list.variants <- vector('list', nrow(genes))

  for(current_gene in 1:nrow(genes)){
    variants_loc <- get_variants_from_locations(locations = paste0(genes[current_gene,"chromosome_name"], ":",
                                                                   genes[current_gene,"start_position"], ":",
                                                                   genes[current_gene,"end_position"]),
                                                verbose = verbose)
    if(length(variants_loc) != 0){
      current_gene_variants <- cbind(ensembl_gene_id = genes[current_gene, "ensembl_gene_id"],
                                     hgnc_symbol = genes[current_gene, "hgnc_symbol"],
                                     variants_loc)

      current_gene_variants_df <- format_output(chr = unlist(current_gene_variants$seq_region_name),
                                                pos =  unlist(current_gene_variants$start),
                                                rsid = unlist(current_gene_variants$id),
                                                ensembl_gene_id = unique(current_gene_variants$ensembl_gene_id),
                                                hgnc_symbol = unique(current_gene_variants$hgnc_symbol))

      list.variants[[current_gene]] <- current_gene_variants_df
    }
  }

  # Removing the NULL elements of the list if exists
  list.variants <- list.variants[!sapply(list.variants, is.null)]
  # Then concatenate the list into a data.frame
  genes_variants <- do.call('rbind', list.variants)
  genes_variants$source <- "omim"

  return(genes_variants)
}
