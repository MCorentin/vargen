# ---- vargen Pipeline ----

#' @title Download the files needed to run \code{\link{vargen_pipeline}}
#' @description Only need to run it once. Will download the following files in
#' "install_dir":
#' \itemize{
#'   \item the latest gwas catalog, eg: gwas_catalog_v1.0.2-associations_e93_r2019-01-11.tsv
#'   \item hg19ToHg38.over.chain.gz (will be unzipped)
#'   \item GTEx_Analysis_v8_eQTL.tar.gz (will be untared)
#'   \item GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz
#'   \item enhancer_tss_associations.bed
#' }
#'
#' @param install_dir the path to the installation folder (default = "./")
#' @param gtex_version the version of gtex to download, only "v7" and "v8" are
#' supported (default = "v8")
#' @param verbose if TRUE will print progress messages (default = FALSE)
#'
#' @return nothing, download files in "install_dir".
#'
#' @examples
#' vargen_install("./", verbose = TRUE)
#' @export
vargen_install <- function(install_dir = "./", gtex_version = "v8", verbose = FALSE){

  if (!file.exists(install_dir)){
    if(verbose) print(paste0("Creating folder '", install_dir, "'"))
    dir.create(install_dir)
  }

  if(verbose) print("Dowloading FANTOM's enhancer tss associations")
  # FANTOM5 TSS associations
  utils::download.file("http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed",
                       destfile = paste0(install_dir, "/enhancer_tss_associations.bed"))

  cat("\n")

  if(verbose) print("Downloading liftOver chain file from ucsc")
  # liftOver file
  utils::download.file(url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                       destfile = paste0(install_dir, "/hg19ToHg38.over.chain.gz"))
  R.utils::gunzip(filename = paste0(install_dir, "/hg19ToHg38.over.chain.gz"),
                  skip = TRUE, remove = TRUE)


  # download gwas file:
  if(verbose) print("Downloading the gwas catalog file from ebi")
  gwasurl <- "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
  # gwasheaders contains the http header from the URL
  gwasheaders <- curlGetHeaders(url = gwasurl, redirect = TRUE, verify = TRUE)
  # From the header, we can extract the line where the filename is
  gwasheaders <- gwasheaders[grep("filename", gwasheaders)]
  # Then we extract the filname from the header line, which should resemble:
  # "Content-Disposition: attachement; filename=gwas_catalog_v1.0-associations_e96_r2019-10-14.tsv\r\n"
  gwasfilename <- unlist(strsplit(x = gwasheaders, split = "filename="))[2]
  # This should remove new lines "\r", "\r\n" or "\n"
  gwasfilename <- gsub("\r?\n|\r", "", gwasfilename)

  # Now that we have the gwas catalog filename, we can download it:
  utils::download.file(url = gwasurl,
                       destfile = paste0(install_dir, "/", gwasfilename))

  cat("\n")

  # Add gtex folder
  if(verbose) print("Downloading GTEx variant association file... This may take a while")
  if(gtex_version == "v7") {
    gtex_filename = "GTEx_Analysis_v7_eQTL.tar.gz"
    gtex_url = "https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz"
  } else if(gtex_version == "v8") {
    gtex_filename = "GTEx_Analysis_v8_eQTL.tar"
    gtex_url = "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar"
  } else {
    stop(paste0("Please set 'gtex_version' as v7 or v8, current value: '", gtex_version, "'"))
  }

  # The tar file is considered as a binary file, so the "mode = wb" option is needed
  # or else there is a "corrupt archive" error during untar.
  utils::download.file(url = gtex_url,
                       destfile = paste0(install_dir, "/", gtex_filename),
                       mode = "wb")
  utils::untar(tarfile = paste0(install_dir, "/", gtex_filename),
               exdir = paste0(install_dir))

  if(file.remove(paste0(install_dir, "/", gtex_filename))){
    if(verbose) print(paste0(install_dir, "/", gtex_filename, " untared and removed succesfully"))
  } else{
    print(paste0("Error while removing ", gtex_filename))
  }

  # Installing GTEx lookup table
  if(verbose) print("Downloading GTEx lookup table... This may take a while")
  gtex_lookup_filename <- "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"
  gtex_lookup_url <- paste0("https://storage.googleapis.com/gtex_analysis_v8/reference/", gtex_lookup_filename)
  utils::download.file(url = gtex_lookup_url,
                       destfile = paste0(install_dir, "/", gtex_lookup_filename),
                       mode = "wb")
}


#' @title Main vargen function, to get the list of variants
#' @description Will get a list of variants related to certain OMIM morbid IDs.
#' Be aware that some of these variants will not be necessarily associated to the
#' phenotype. We advise to filter the results by annotation ("CADD phred score",
#' "snpEff impact" etc...). If you want a smaller list of variants that are all
#' associated with the disease, then run \code{\link{get_variants_from_phenotypes}}
#' The variants are fetched from the following sources:
#' \itemize{
#'   \item OMIM: get variants on the genes related to the disease
#'   \item FANTOM5: get the variants on the enhancers of the OMIM genes
#'   \item GTEx: get variants impacting the expression of the OMIM genes in specific tissues.
#'   \item GWAS: get variants related to the phenotype of interest from the gwas catalog
#' }
#' The pipeline will also annotate the variants using \code{\link[myvariant]{getVariants}}
#'
#' @param vargen_dir directory with the following file (can be generated with
#' \code{\link{vargen_install}})
#' @param omim_morbid_ids a vector containing the omim morbid id(s) of the phenotype(s)
#' of interest. You can search on the Online Mendelian Inheritance in Man website
#' (https://www.omim.org/) or use \code{\link{list_omim_accessions}}
#' @param fantom_corr the minimum correlation (z-score) to consider a FANTOM5
#' enhancer/gene association valid (default: 0.25).
#' A z-score greater than 0 represents an element greater than the mean, this
#' means that this association has more correlation than random motifs.
#' @param outdir the output directory, some files will be written during the
#' running of this function
#' @param gtex_tissues a vector containing the name of the "signif_variant_gene_pairs.txt.gz"
#' files. Output from \code{\link{select_gtex_tissues}} can be used.
#' @param gwas_traits a vector with the trait of interest (as characters). The list
#' of available traits can be obtained with \code{\link{list_gwas_traits}}
#' @param gene_mart optional, a connection to ensembl gene mart, can be created
#' using \code{\link{connect_to_gene_ensembl}} (If missing this function will be
#' used to create the connection).
#' @param snp_mart optional, a connection to ensembl snp mart, can be created
#' using \code{\link{connect_to_snp_ensembl}} (If missing this function will be
#' used to create the connection).
#' @param verbose if TRUE, will print progress messages (default: FALSE)
#' @return a data.frame with the variants fetched from OMIM, FANTOM5, GTEx and GWAS.
#' The data.frame will contain the following columns:
#' \itemize{
#'   \item chr (chromosome)
#'   \item pos (position of the variant)
#'   \item rsid (variant ID)
#'   \item ensembl_gene_id ("gene id" of the gene associated with the variant)
#'   \item hgnc_symbol ("hgnc symbol" of the gene associated with the variant)
#'   \item source ("omim", "fantom5", "gtex" or "gwas")
#'   \item trait (the "omim ids" seperated by ';' for omim,fantom and gtex variants and the gwas trait
#'   for the gwas variants).
#' }
#'
#' @examples
#' vargen_install("./vargen_data/")
#'
#' # Simple query
#' DM1_simple <- vargen_pipeline(vargen_dir = "./vargen_data/", omim_morbid_ids = "222100",
#'                               fantom_corr = 0.25, outdir = "./", verbose = TRUE)
#'
#'
#' # Query with gtex and gwas
#' pancreas_tissues <- select_gtex_tissues(gtex_dir = "./vargen_data/GTEx_Analysis_v8_eQTL/",
#'                                         tissues_query = "pancreas")
#'
#' # list_gwas_traits("diabetes")
#'
#' DM1 <- vargen_pipeline(vargen_dir = "./vargen_data/", omim_morbid_ids = "222100",
#'                        fantom_corr = 0.25, outdir = "./",
#'                        gtex_tissues = pancreas_tissues,
#'                        gwas_traits = "Type 1 diabetes", verbose = TRUE)
#' @export
vargen_pipeline <- function(vargen_dir, omim_morbid_ids, fantom_corr = 0.25,
                            outdir = "./", gtex_tissues, gwas_traits,
                            gene_mart, snp_mart, verbose = FALSE) {

  if(missing(omim_morbid_ids)){
    stop("Please provide at least one OMIM morbid id. Stopping now.")
  }

  if(!dir.exists(outdir)){
    if(verbose) print(paste0("Creating folder '", outdir, "'"))
    dir.create(outdir)
  }

  #_____________________________________________________________________________
  # Loading the necessary resources
  #_____________________________________________________________________________
  if(missing(gene_mart) || class(gene_mart) != 'Mart'){
    if(verbose) print("Connecting to the gene mart...")
    gene_mart <- connect_to_gene_ensembl()
    warning("Gene mart not provided (or not a valid Mart object).",
            "We used one from connect_to_gene_ensembl() instead.")
  }

  if(missing(snp_mart) || class(snp_mart) != 'Mart'){
    if(verbose) print("Connecting to the snp mart...")
    snp_mart <- connect_to_snp_ensembl()
    warning("Snp mart not provided (or not a valid Mart object).",
            "We used one from connect_to_snp_ensembl() instead.")
  }

  # no gwas traits = no need to generate the gwas_cat object
  if(!missing(gwas_traits)){
    if(verbose) print("Building the gwascat object...")
      gwas_cat <- create_gwas(vargen_dir)
    # Check if the gwas traits are in the gwas catalog:
    for(trait in gwas_traits){
      if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)){
        stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
      }
    }
  }

  if(verbose) print(paste0("Reading the enhancer tss association file for FANTOM5... '" ,
                           vargen_dir, "/enhancer_tss_associations.bed'"))
  fantom_df <- prepare_fantom(enhancer_tss_association = paste0(vargen_dir,
                                                                "/enhancer_tss_associations.bed"))

  hg19ToHg38.over.chain <- paste0(vargen_dir, "/hg19ToHg38.over.chain")
  if(!file.exists(hg19ToHg38.over.chain)){
    stop(paste0("Can not read: ", hg19ToHg38.over.chain, ", stopping now."))
  }

  gtex_lookup_file <- paste0(vargen_dir, "/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")
  if(!file.exists(gtex_lookup_file)){
    stop(paste0("Can not read: ", gtex_lookup_file, ", stopping now."))
  }

  #_____________________________________________________________________________
  # Getting variants from genes related to OMIM disease
  #_____________________________________________________________________________
  if(verbose) print("Starting the pipeline...")

  omim_all_genes <- data.frame()
  master_variants <- data.frame()

  for(omim_morbid in omim_morbid_ids){
    if(verbose) print(paste0("Getting genes for OMIM: ", omim_morbid))
    # First: get all the genes linked to the different omim ids:
    omim_genes <- get_omim_genes(omim_morbid, gene_mart)

    if(nrow(omim_genes) == 0){
      if(verbose) print(paste0("warning: no genes found for omim id: ", omim_morbid))
    } else {
      omim_all_genes <- rbind(omim_all_genes, omim_genes)
    }
  }

  if(nrow(omim_all_genes) > 0){
    omim_all_genes <- unique(omim_all_genes)

    # We get the variants on the genes:
    genes_variants <- get_genes_variants(genes = omim_all_genes, verbose = verbose)

    if(length(genes_variants) != 0) master_variants <- rbind(master_variants,
                                                             genes_variants)

    # We get the variants on the enhancers of the genes:
    fantom_variants <- get_fantom5_variants(fantom_df = fantom_df,
                                            omim_genes = omim_all_genes,
                                            corr_threshold = fantom_corr,
                                            hg19ToHg38.over.chain = hg19ToHg38.over.chain,
                                            verbose = verbose)

    if(length(fantom_variants) != 0) master_variants <- rbind(master_variants,
                                                              fantom_variants)

    # master_variants$trait <- omim_morbid
    # print(head(master_variants))

    if(verbose) print(paste0("Writing the list of genes to: ", outdir, "/genes_info.tsv"))
    # We write the list of genes in a file.
    utils::write.table(x = omim_all_genes, quote = FALSE, sep = "\t", row.names = FALSE,
                       file = paste0(outdir, "/genes_info.tsv"))

    #_____________________________________________________________________________
    # Getting variants associated with change of expression in GTEx (need tissues as input)
    #_____________________________________________________________________________
    if(!missing(gtex_tissues)){
      if(verbose) print("Getting the GTEx variants...")
      gtex_variants <- get_gtex_variants(tissue_files = gtex_tissues,
                                         omim_genes = omim_all_genes,
                                         gtex_lookup_file = gtex_lookup_file,
                                         snp_mart = snp_mart,
                                         verbose = verbose)

      if(length(gtex_variants) != 0) master_variants <- rbind(master_variants,
                                                              gtex_variants)

    } else{
      if(verbose) print("No values for 'gtex_tissues', skipping GTEx step...")
    }

      # Adding the traits for each variant.
      master_variants$trait <- ""

      for(gene in unique(master_variants$hgnc_symbol)){
        # First get the list of omim ids for this gene:
        gene_omims <- paste0(unique(omim_all_genes[omim_all_genes$hgnc_symbol == gene,]$mim_morbid_accession),
                             collapse = ";")

        # Then add it to the "trait" column
        master_variants[master_variants$hgnc_symbol == gene,]$trait <- gene_omims
      }
  }

  #_____________________________________________________________________________
  # Getting variants associated to the disease in the gwas catalog
  #_____________________________________________________________________________
  # GWAS variants (only if list of gwas traits were given)
  if(!missing(gwas_traits)){
    if(verbose) print("Getting the gwas variants,,,")
    master_variants <- rbind(master_variants, get_gwas_variants(gwas_traits = gwas_traits,
                                                                gwas_cat = gwas_cat))
  } else{
    if(verbose) print("No values for 'gwas_traits', skipping gwas step...")
  }

  # writing the variants data.frame to a file
  if(verbose) print(paste0("Writing the variants to ",
                           paste0(outdir, "/vargen_variants.tsv")))
  utils::write.table(x = unique(master_variants), append = FALSE, quote = FALSE,
                     sep = "\t", row.names = FALSE,
                     file = paste0(outdir, "/vargen_variants.tsv"))

  return(unique(master_variants))
}


#' @title Get variants related to specific genes.
#' @description Alternative to vargen_pipeline for a set of custom genes.
#' The variants are fetched from the following sources:
#' \itemize{
#'   \item GENE: get variants on the genes
#'   \item FANTOM5: get the variants on the enhancers / promoters of the genes
#'   \item GTEx: get variants impacting the expression of the genes in specific tissues.
#'   \item GWAS: get variants related to the phenotype of interest from the gwas catalog
#' }
#'
#' @param vargen_dir directory with the following file (can be generated with
#' \code{\link{vargen_install}})
#' @param gene_ids list of ensembl gene IDs of interest
#' @param fantom_corr the minimum correlation (z-score) to consider a FANTOM5
#' enhancer/gene association valid (default: 0.25).
#' A z-score greater than 0 represents an element greater than the mean, this
#' means that this association has more correlation than random motifs.
#' @param outdir the output directory, some files will be written during the
#' running of this function
#' @param gtex_tissues a vector containing the name of the "signif_variant_gene_pairs.txt.gz"
#' files. Output from \code{\link{select_gtex_tissues}} can be used.
#' @param gwas_traits a vector with the trait of interest (as characters). The list
#' of available traits can be obtained with \code{\link{list_gwas_traits}}
#' @param gene_mart optional, a connection to ensembl gene mart, can be created
#' using \code{\link{connect_to_gene_ensembl}} (If missing this function will be
#' used to create the connection).
#' @param snp_mart optional, a connection to ensembl snp mart, can be created
#' using \code{\link{connect_to_snp_ensembl}} (If missing this function will be
#' used to create the connection).
#' @param verbose if TRUE, will print progress messages (default: FALSE)
#'
#' @return a data.frame with the variants fetched from OMIM, FANTOM5, GTEx and GWAS.
#' The data.frame will contain the following columns:
#' \itemize{
#'   \item chr (chromosome)
#'   \item pos (position of the variant)
#'   \item rsid (variant ID)
#'   \item ensembl_gene_id ("gene id" of the gene associated with the variant)
#'   \item hgnc_symbol ("hgnc symbol" of the gene associated with the variant)
#'   \item source ("omim", "fantom5", "gtex" or "gwas")
#'   \item trait (an empty string for the gene, fantom and gtex variants, since
#'   this comes from a list of genes, no omim id is associated with them.
#'   Contains the gwas trait for the gwas variants. )
#' }
#'
#' @examples
#' # Simple query
#' vargen_install(install_dir = "./vargen_data/")
#'
#' vargen_custom(vargen_dir = "./vargen_data/",
#'               gene_ids = c("ENSG00000166603", "ENSG00000155846"),
#'               outdir = "./", verbose = TRUE)
#'
#' # With gwas and gtex
#' adipose_tissues <- select_gtex_tissues(gtex_dir = "./vargen_data/GTEx_Analysis_v8_eQTL/",
#'                                        tissues_query = "adipose")
#' vargen_custom(vargen_dir = "./vargen_data/",
#'               gene_ids = c("ENSG00000166603", "ENSG00000155846"),
#'               outdir = "./", verbose = TRUE,
#'               gtex_tissues = adipose_tissues,
#'               gwas_traits = "Obesity")
#' @export
vargen_custom <- function(vargen_dir, gene_ids, fantom_corr = 0.25, outdir = "./",
                          gtex_tissues, gwas_traits, gene_mart, snp_mart,
                          verbose = FALSE) {

  if (!dir.exists(outdir)){
    if(verbose) print(paste0("Creating folder '", outdir, "'"))
    dir.create(outdir)
  }

  #_____________________________________________________________________________
  # Loading the necessary resources
  #_____________________________________________________________________________
  if(missing(gene_mart) || class(gene_mart) != 'Mart'){
    if(verbose) print("Connecting to the gene mart...")
    gene_mart <- connect_to_gene_ensembl()
    warning("Gene mart not provided (or not a valid Mart object).",
            "We used one from connect_to_gene_ensembl() instead.")
  }

  if(missing(snp_mart) || class(snp_mart) != 'Mart'){
    if(verbose) print("Connecting to the snp mart...")
    snp_mart <- connect_to_snp_ensembl()
    warning("Snp mart not provided (or not a valid Mart object).",
            "We used one from connect_to_snp_ensembl() instead.")
  }

  # no gwas traits = no need to generate the gwas object
  if(!missing(gwas_traits)){
    if(verbose) print("Building the gwascat object...")
    gwas_cat <- create_gwas(vargen_dir)
    # Check if the gwas traits are in the gwas catalog:
    for(trait in gwas_traits){
      if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)) stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
    }
  }

  if(verbose) print(paste0("Reading the enhancer tss association file for FANTOM5... '" ,
                           vargen_dir, "/enhancer_tss_associations.bed'"))
  fantom_df <- prepare_fantom(enhancer_tss_association = paste0(vargen_dir, "/enhancer_tss_associations.bed"))


  hg19ToHg38.over.chain <- paste0(vargen_dir, "/hg19ToHg38.over.chain")
  if(!file.exists(hg19ToHg38.over.chain)){
    stop(paste0("Can not read: ", hg19ToHg38.over.chain, ", stopping now."))
  }

  gtex_lookup_file <- paste0(vargen_dir, "/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")
  if(!file.exists(gtex_lookup_file)){
    stop(paste0("Can not read: ", gtex_lookup_file, ", stopping now."))
  }

  master_variants <- data.frame()
  #_____________________________________________________________________________
  # Getting variants from genes related to OMIM disease
  #_____________________________________________________________________________
  # Get the genes related to the phenotype:
  genes_info <- biomaRt::getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                              "start_position", "end_position",
                                              "hgnc_symbol"),
                              filters = c("ensembl_gene_id"),
                              values = c(gene_ids),
                              mart = gene_mart, uniqueRows = TRUE)

  if(nrow(genes_info) > 0){

    if(verbose) print(paste0("Writing the list of genes to: ",
                             outdir, "/custom_genes_info.tsv"))

    utils::write.table(x = genes_info, quote = FALSE, sep = "\t",
                       row.names = FALSE,
                       file = paste0(outdir, "/custom_genes_info.tsv"))

    # Get the variants on the genes:
    master_variants <- get_genes_variants(genes = genes_info, verbose = verbose)
    # Replacing "omim" by "gene"
    master_variants$source <- "gene"

    #_____________________________________________________________________________
    # Getting variants on the enhancers of the OMIM genes, using FANTOM5
    #_____________________________________________________________________________
    fantom_variants <- get_fantom5_variants(fantom_df = fantom_df,
                                            omim_genes = genes_info,
                                            corr_threshold = fantom_corr,
                                            hg19ToHg38.over.chain = hg19ToHg38.over.chain,
                                            verbose = verbose)

    if(length(fantom_variants) != 0) master_variants <- rbind(master_variants,
                                                              fantom_variants)

    #_____________________________________________________________________________
    # Getting variants associated with change of expression in GTEx (need tissues as input)
    #_____________________________________________________________________________
    if(!missing(gtex_tissues)){
      if(verbose) print("Getting the GTEx variants...")
      gtex_variants <- get_gtex_variants(tissue_files = gtex_tissues,
                                         omim_genes = genes_info,
                                         gtex_lookup_file = gtex_lookup_file,
                                         snp_mart = snp_mart,
                                         verbose = verbose)

      if(length(gtex_variants) != 0) master_variants <- rbind(master_variants,
                                                              gtex_variants)

    } else{
      print("No values for 'gtex_tissues', skipping GTEx step...")
    }

    master_variants$trait <- ""

  } else {
    print(paste0("No genes found for: ", paste(gene_ids, collapse = ", ")))
  }

  #_____________________________________________________________________________
  # Getting variants associated to the disease in the gwas catalog
  #_____________________________________________________________________________
  # GWAS variants (only if list of gwas traits were given)
  if(!missing(gwas_traits)){
    if(verbose) print("Getting the gwas variants...")
    master_variants <- rbind(master_variants,
                             get_gwas_variants(gwas_traits = gwas_traits,
                                               gwas_cat = gwas_cat))
  } else{
    if(verbose) print("No values for 'gwas_traits', skipping gwas step...")
  }

  # writing the variants data.frame to a file
  if(verbose) print(paste0("Writing the variants to ",
                           paste0(outdir, "/custom_vargen_variants.tsv")))

  utils::write.table(x = unique(master_variants), append = FALSE, quote = FALSE,
                     sep = "\t", row.names = FALSE,
                     file = paste0(outdir, "/custom_vargen_variants.tsv"))

  return(unique(master_variants))
}


#' @title Generate a plot of the variants on the OMIM genes
#' @description The plot contains 4 tracks:
#' \itemize{
#'   \item  The chromosome, with a red marker on the gene location
#'   \item  The ensembl transcripts
#'   \item  The variant consequences, grouped by type (eg: INTRONIC, STOP LOST etc...),
#' each green bar represent a variant
#'   \item  The cadd phred score, each dot represent a variant
#' }
#'
#' @param annotated_snps a data.frame from the merging of the outputs of
#' \code{\link{vargen_pipeline}} and \code{\link{annotate_variants}}
#' @param outdir the directory that will contain the plots
#' @param rsid_highlight optional, a vector of rsids, which will be plotted in red
#' on the cadd track.
#' @param device only "pdf" and "png" are supported now.
#' @param verbose if TRUE, will display progress messages
#' @param gene_mart optional, a connection to ensembl gene mart, can be created
#' using \code{\link{connect_to_gene_ensembl}}. If missing this function will be
#' used to create the connection.
#'
#' @return nothing, create the plots in "outdir". One file per gene, the name format
#' is <hgnc_symbol>_<ensembl_id>_GVIZ
#'
#' @examples
#' vargen_install("./vargen_data/")
#'
#' # Simple query
#' gene_mart <- connect_to_gene_ensembl()
#' DM1_simple <- vargen_pipeline(vargen_dir = "./vargen_data/", omim_morbid_ids = "222100",
#'                               fantom_corr = 0.25, outdir = "./", verbose = TRUE)
#'
#' vargen_visualisation(annotated_snps = DM1_simple, verbose = TRUE,
#'                      outdir = "./DM1_gviz", device = "png",
#'                      gene_mart = gene_mart)
#' @export
vargen_visualisation <- function(annotated_snps, outdir = "./", rsid_highlight,
                                 device = "pdf", verbose = FALSE, gene_mart){
  if(device != "png" && device != "pdf"){
    stop("Please specify device as png or pdf")
  }

  if (!file.exists(outdir)){
    if(verbose) print(paste0("Creating folder '", outdir, "'"))
    dir.create(outdir)
  }

  if(missing(gene_mart) || class(gene_mart) != 'Mart'){
    gene_mart <- connect_to_gene_ensembl()
    warning("Gene mart not provided (or not a valid Mart object).",
            "We used one from connect_to_gene_ensembl() instead.")
  }

  bckg_col <- "#D95F02"

  # Get genes position for omim genes (not gwas / gtex genes as we currently limit
  # the plot to start and stop of the gene, and gwas snps often lies outside these coordinates)
  omim_genes <- annotated_snps[annotated_snps$source == "omim","ensembl_gene_id"]
  genes_list <- biomaRt::getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                              "start_position", "end_position",
                                              "hgnc_symbol"),
                               filters = c("ensembl_gene_id"),
                               values = unique(omim_genes),
                               mart = gene_mart, uniqueRows = TRUE)

  # UCSC chr format
  genes_list$chromosome_name <- paste0("chr", genes_list$chromosome_name)

  # From get_genes_from_omim output:
  genes <- genes_list[,2:4]
  colnames(genes) <- c("chromosome", "start", "end")

  regions <- GenomicRanges::makeGRangesFromDataFrame(genes)
  GenomeInfoDb::genome(regions) <- "hg38"
  GenomeInfoDb::seqlevelsStyle(regions) <- "UCSC"

  for(i in 1:nrow(genes_list)){
    if(verbose) print(genes_list[i,]$hgnc_symbol)

    # Subselect the snps on the gene locus only:
    track_variants <- subset(annotated_snps, annotated_snps$chr == genes_list[i,]$chromosome_name)
    track_variants <- subset(track_variants,track_variants$pos >= genes_list[i,]$start_position)
    track_variants <- subset(track_variants,track_variants$pos <= genes_list[i,]$end_position)
    track_variants <- track_variants[!(is.na(track_variants$consequence)),]
    track_variants <- track_variants[!(track_variants$consequence == ""),]
    track_variants$chr <- as.vector(track_variants$chr)

    if(nrow(track_variants) != 0){
      # Name of consequences are displayed on the left of the track, so we leave 10%
      # of the gene length as a margin on the left to avoid cutting names.
      gene_length <- genes_list[i,]$end_position - genes_list[i,]$start_position
      start_plot <- min(track_variants$pos) - (15/100) * gene_length
      end_plot <- max(track_variants$pos) + 10

      # Track of variants consequences
      conseq_track <- Gviz::AnnotationTrack(background.title = bckg_col,
                                            start = track_variants$pos,
                                            end = track_variants$pos,
                                            chromosome = unique(track_variants$chr),
                                            group = track_variants$consequence,
                                            genome = "hg38", col = "#1B9E77",
                                            fill = "#1B9E77", fontsize = 18,
                                            name = "Variant Consequence")

      # Need to have the variants in position order for the data track
      cadd_track_variants <- track_variants[!(is.na(track_variants$cadd_phred)),]
      cadd_track_variants <- cadd_track_variants[order(cadd_track_variants$pos),]

      if(nrow(cadd_track_variants) != 0){
        # Track for cadd phrad score, width of 1 because of SNP
        cadd_track <- Gviz::DataTrack(background.title = bckg_col,
                                      data = cadd_track_variants$cadd_phred,
                                      start = cadd_track_variants$pos,
                                      width = 1, fontsize = 16,
                                      chromosome = unique(cadd_track_variants$chr),
                                      genome = "hg38", name = "CADD Phred score",
                                      ylim = c(0,40), group = cadd_track_variants$group)

        # Overlaying a second data track to highlight some snps:
        # Maybe it would be better to use the "group" paramater of the
        # DataTrack rather than overlapping data tracks.
        cadd_highlights <- data.frame()
        if(!missing(rsid_highlight)){
          cadd_highlights <- cadd_track_variants[cadd_track_variants$rsid %in% rsid_highlight,]
          if(nrow(cadd_highlights) != 0){
            cadd_highlight_track <- Gviz::DataTrack(background.title = bckg_col,
                                                    data = cadd_highlights$cadd_phred,
                                                    start = cadd_highlights$pos,
                                                    width = 1, fontsize = 16,
                                                    group = cadd_highlights$rsid,
                                                    chromosome = unique(cadd_highlights$chr),
                                                    genome = "hg38", name = "CADD Phred score",
                                                    col = "red", ylim = c(0,40))
          }
        }
      }

      gene_track <- Gviz::UcscTrack(background.title = bckg_col, fontsize = 18,
                                    genome = "hg38", track = "knownGene",
                                    chromosome = genes_list[i,]$chromosome_name,
                                    from = start_plot, to = end_plot,
                                    trackType = "GeneRegionTrack",
                                    rstarts = "exonStarts", rends = "exonEnds",
                                    gene = "name", symbol = "name", transcript = "name",
                                    strand = "strand", fill = "#8282d2",
                                    transcriptAnnotation = "transcript",
                                    name = genes_list[i,"hgnc_symbol"])

      itrack <- Gviz::IdeogramTrack(genome="hg38", chromosome = genes_list[i,]$chromosome_name)
      axis_track <- Gviz::GenomeAxisTrack()

      plot_name <- paste0(outdir,"/",genes_list[i,"hgnc_symbol"], "_", genes_list[i,"ensembl_gene_id"], "_GVIZ")
      if(device == "pdf") grDevices::pdf(paste0(plot_name,".pdf"), width = 15, height = 10)
      if(device == "png") grDevices::png(paste0(plot_name,".png"), width = 15, height = 10, units = "in", res = 300)

      tracklist <- list(itrack, axis_track, gene_track, conseq_track)

      # If some variants are to be highlighted, we combine the cadd track and the
      # highlight track into one "overlay track". Else, if we have cadd scores
      # we just plot the "cadd_track" on its own.
      if(nrow(cadd_highlights) != 0){
        highlight_track <- Gviz::OverlayTrack(background.title = bckg_col,
                                              trackList = list(cadd_track,
                                                               cadd_highlight_track))
        tracklist <- append(tracklist, highlight_track)
      } else {
        if(nrow(cadd_track_variants) != 0){
          tracklist <- append(tracklist, cadd_track)
        }
      }

      Gviz::plotTracks(tracklist,
                       from = start_plot,
                       to = end_plot,
                       groupAnnotation = "group")
      grDevices::dev.off()
    }
  }
}
