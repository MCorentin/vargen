#---- Utils ----

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
#' # Connect with the "uswest" mirror
#' uswest_gene_mart <- connect_to_gene_ensembl(mirror="uswest")
#' @export
connect_to_gene_ensembl <- function(mirror = "www"){
  gene_mart <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                                   host = "www.ensembl.org",
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
#' # Connect with the "uswest" mirror
#' uswest_ensembl <- connect_to_snp_ensembl(mirror="uswest")
#' @export
connect_to_snp_ensembl <- function(mirror = "www"){
  snp_mart <- biomaRt::useEnsembl(biomart = "snp",
                                host = "www.ensembl.org",
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
  variants_at_loc <- data.frame()
  ensembl_rest_limit <- 0

  for(curr_loc in locations){
    # ensembl has a limit of 15 requests / seconds
    if(ensembl_rest_limit >= 14){
      if(verbose) print("waiting... (to not reach the limit of queries on ensembl REST API)")
      Sys.sleep(1)
      ensembl_rest_limit <- 0
    }

    query <- paste0("/overlap/region/human/",curr_loc, "?content-type=application/json;feature=variation")
    get_output <- httr::GET(paste("https://rest.ensembl.org", query, sep = ""),
                            httr::content_type("application/json"))
    httr::stop_for_status(get_output)
    variants_at_loc <- rbind(variants_at_loc, jsonlite::fromJSON(
                          jsonlite::toJSON(httr::content(get_output))))

    #if(verbose) print(paste0("location '", curr_loc, "' done !"))

    ensembl_rest_limit <- ensembl_rest_limit + 1
  }

  return(variants_at_loc)
}


#' @title Format chromosome names
#' @description Transform 1 to 22, X, Y and MT into chr1 to chr22, chrX, chrY and chMT.
#' Regex needed because some region are named CHR17_... and we do not want
#' chrCHR17...
#'
#' @param chr a vector of chromosome names
#' @return the formatted vector of chromosome names (eg: "chr1" instead of "1")
format_chr <- function(chr){
  # This should match 1 to 22 and X, Y, MT. We need to replace them as
  # chr1, chrX etc....
  chr_to_rename <- grep(pattern = "^\\d{1,2}$|^X$|^Y$|^MT$", chr, ignore.case = T)
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
  # First, use "getVariants" to annotate the snps
  rsid_annotated <- myvariant::getVariants(hgvsids = rsid, verbose = verbose,
                                           fields =  c("cadd", "dbnsfp", "clinvar", "snpeff"))

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
  
  
  rsid_annotated_df <- data.frame(rsid = unlist(rsid_annotated$query),
                                  cadd_phred = cadd_phred,
                                  fathmm_xf_score = fathmm_score,
                                  fathmm_xf_pred = fathmm_pred,
                                  annot_type = annot_type,
                                  consequence = consequence,
                                  clinical_significance = clinical_significance,
                                  snpeff_ann = snpeff,
                                  stringsAsFactors = FALSE)
  

  #Replacing the snp eff annotation from "c(\"MODIFIER\", \"MODIFIER\")" to "MODIFIER;MODIFIER"
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


#' @title Generate a plot of the variants on the OMIM genes
#' @description The plot contains 4 tracks:
#' \itemize{
#'   \item
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
#' DM1_simple <- vargen_pipeline(vargen_dir = "./vargen_data/", omim_morbid = "222100",
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
  if(missing(gene_mart)) gene_mart <- connect_to_gene_ensembl()

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


# ---- varphen Pipeline ----

#' @title From keywords, get a list of phenotypes to use with snp mart
#' @description The list of phenotypes obtained can be used with
#' \code{\link{get_variants_from_phenotypes}} to get a list of variants associated
#' with the phenotypes.
#'
#' @param keywords a keyword, all the phenotypes contaning this keyword will be
#' returned (based on \code{\link[base]{grep}} with ignore.case)
#' @param snp_mart a connection to ensembl snp mart, can be generated from
#' \code{\link{connect_to_snp_ensembl}}
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
  if(missing(snp_mart)) snp_mart <- connect_to_snp_ensembl()

  # Get all the possible values for the phenotype filter:
  biomart_filters <- biomaRt::filterOptions(filter = "phenotype_description",
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

  if(missing(snp_mart)) snp_mart <- connect_to_snp_ensembl()

  pheno_variants <- biomaRt::getBM(attributes = c("chr_name", "chrom_start",
                                                  "refsnp_id", "refsnp_source",
                                                  "associated_variant_risk_allele",
                                                  "phenotype_description",
                                                  "clinical_significance", "validated"),
                                   filters = c("phenotype_description"),
                                   values = c(phenotypes),
                                   mart = snp_mart
  )
  pheno_variants$chr_name <- format_chr(pheno_variants$chr_name)
  return(pheno_variants)
}


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
  if(missing(gene_mart)) gene_mart <- connect_to_gene_ensembl()

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
                                ignore.case = T),]

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
  if(missing(gene_mart)) gene_mart <- connect_to_gene_ensembl()

  mim_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                             "start_position", "end_position",
                                             "hgnc_symbol", "mim_morbid_accession",
                                             "mim_morbid_description"),
                              filters = c("mim_morbid_accession"),
                              values = c(omim_ids),
                              mart = gene_mart, uniqueRows = TRUE)

  if(is.null(gene_mart)) paste0("/!\\ no genes found for:", paste(omim_ids, collapse = ", "))

  if(!is.null(gene_mart)) mim_genes$chromosome_name <- format_chr(chr = mim_genes$chromosome_name)
  
  return(mim_genes)
}


#' @title Get the variants on the OMIM genes
#' @description Get the variants on the gene locations obtained from
#' \code{\link{get_omim_genes}}
#'
#' @param omim_genes list of omim genes from \code{\link{get_omim_genes}}
#' @param verbose if true, will print progress information (default: FALSE)
#'
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
get_omim_variants <- function(omim_genes, verbose = FALSE){
  # Get the variants on the OMIM genes:
  omim_variants <- data.frame()

  # Used to avoid growing the data frame for every gene, we only do it at the end
  # with: "do.call('rbind', list.variants)"
  list.variants <- vector('list', nrow(omim_genes))

  for(gene in 1:nrow(omim_genes)){
    variants_loc <- get_variants_from_locations(locations = paste0(omim_genes[gene,"chromosome_name"], ":",
                                                                   omim_genes[gene,"start_position"], ":",
                                                                   omim_genes[gene,"end_position"]),
                                                verbose = verbose)
    if(length(variants_loc) != 0){
      gene_variants <- cbind(ensembl_gene_id = omim_genes[gene, "ensembl_gene_id"],
                             hgnc_symbol = omim_genes[gene, "hgnc_symbol"],
                             variants_loc)
      
      gene_variants_df <- format_output(chr = unlist(gene_variants$seq_region_name),
                                        pos =  unlist(gene_variants$start),
                                        rsid = unlist(gene_variants$id),
                                        ensembl_gene_id = gene_variants$ensembl_gene_id,
                                        hgnc_symbol = gene_variants$hgnc_symbol)
      
      list.variants[[gene]] <- gene_variants_df
    }
  }

  # Removing the NULL elements of the list if exists
  list.variants <- list.variants[!sapply(list.variants, is.null)]
  # Then concatenate the list into a data.frame
  omim_variants <- do.call('rbind', list.variants)
  omim_variants$source <- "omim"

  return(omim_variants)
}


# ---- FANTOM ----

#' @title Generate data.frame from FANTOM5 enhancer file
#' @description Prepare FANTOM5 for \code{\link{get_fantom5_enhancers_from_hgnc}}
#'
#' @param enhancer_tss_association the "enhancer_tss_associations.bed" file from FANTOM5.
#' The file can be downloaded here:
#' http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed
#' @return a data.frame containing the file information
prepare_fantom <- function(enhancer_tss_association) {
  # TSS = Transcripton Start Site
  fantom <- utils::read.delim(enhancer_tss_association, skip = 1, stringsAsFactors = FALSE)

  fantom_df <- as.data.frame(splitstackshape::cSplit(fantom, splitCols="name",
                                                     sep=";", direction="wide"))

  locs <- strsplit(as.character(fantom_df$name_1),"[:-]")
  fantom_df$chr <- sapply(locs,"[",1)
  fantom_df$start <- as.numeric(sapply(locs,"[",2))
  fantom_df$end <- as.numeric(sapply(locs,"[",3))
  fantom_df$symbol <- fantom_df$name_3
  fantom_df$corr <- sub("R:","",fantom_df$name_4)
  fantom_df$fdr <- sub("FDR:","",fantom_df$name_5)

  return(fantom_df)
}


#' @title Get the enhancers associated to certain genes from FANTOM5
#' @description from the FANTOM5 dataset, get the enhancers associated to the
#' genes using the HGNC symbols, the association must pass the correlation threshold given
#' by the user. Is used internally by \code{\link{get_fantom5_variants}}
#'
#' @param fantom_df the output of \code{\link{prepare_fantom}}
#' @param hgnc_symbols vector of HUGO ids for the genes of interest
#' @param corr_threshold the minimum correlation (z-score) to consider a
#' enhancer/gene assocation valid (default: 0.25).
#' A z-score greater than 0 represents an element greater than the mean, this means
#' that this association has more correlation than random motifs.
#' @return a subset of the FANTOM5 data.frame containing the information about the
#' enhancers. The data.frame contains the following columns
#' \itemize{
#'   \item chr (chromosome)
#'   \item start (start of the enhancers)
#'   \item end (end of the enhancers)
#'   \item symbol (HGNC symbol of the gene associated to the enhancer)
#'   \item corr (correlation z-score from FANTOM5)
#'   \item fdr (False Discovery Rate)
#' }
#'

# vargen_install(install_dir = "./vargen_data/")
# fantom_df <- prepare_fantom("./vargen_data/enhancer_tss_associations.bed")
# get_fantom5_enhancers_from_hgnc(fantom_df = fantom_df,
#                                 hgnc_symbols = "FOXP3",
#                                 corr_threshold = 0.30)
get_fantom5_enhancers_from_hgnc <- function(fantom_df, hgnc_symbols, corr_threshold = 0.25){ #fdr_threshold = 1e-5) {
  # Get the enhancers with a decent level of correlations & significance
  # Correlation is Pearson as a z-score = "(pearson corr - (mean of random motifs)) / std(pearson of random motifs)"
  # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0560-6
  # A z-score greater than 0 represents an element greater than the mean = more correlation than random motifs ?
  enhancers_corr <- unique(subset(fantom_df,
                                  corr >= corr_threshold,# & fdr < fdr_threshold,
                                  select = c("chr", "start", "end", "symbol", "corr", "fdr"))
  )
  # Get list of enhancer related to list of genes
  return(enhancers_corr[(enhancers_corr$symbol %in% hgnc_symbols),])
}


#' @title Get variants on the enhancers of the list of genes given as input
#' @description FANTOM5 is used to get the enhancers of the genes, then the variants
#' located on the enhancers are fetched with \code{\link{get_variants_from_locations}}
#' It is also possible to specify a correlation threshold to limit the number
#' of association to output (default: 0.25). The enhancer tss association file
#' is based on hg19, so liftOver is needed.
#'
#' @param fantom_df the output of \code{\link{prepare_fantom}}
#' @param omim_genes output from \code{\link{get_omim_genes}}
#' @param corr_threshold the minimum correlation (z-score) to consider a
#' enhancer/gene association valid (default: 0.25). A z-score greater than 0
#' represents an element greater than the mean, this means that this association
#'  has more correlation than random motifs.
#' @param hg19ToHg38.over.chain the chain file to liftOver locations from
#' hg19 to hg38.
#' @param verbose if true, will print progress information (default: FALSE)
#'
#' @return a data.frame with information about the variants located on the enhancers.
#' The data.frame will contain the following columns:
#' \itemize{
#'   \item chr (chromosome)
#'   \item pos (position of the variant)
#'   \item rsid (variant ID)
#'   \item ensembl_gene_id ("gene id" of the gene associated with the variant)
#'   \item hgnc_symbol ("hgnc symbol" of the gene associated with the variant)
#'   \item source (here the value will be "fantom5")
#' }
#
# gene_mart <- connect_to_gene_ensembl()
# DM1_genes <- get_omim_genes(omim_ids = "222100", gene_mart = gene_mart)
#
# vargen_install(install_dir = "./vargen_data/")
# fantom_df <- prepare_fantom("./vargen_data/enhancer_tss_associations.bed")
# get_fantom5_variants(fantom_df, DM1_genes, 0.25, "hg19ToHg38.over.chain")
get_fantom5_variants <- function(fantom_df, omim_genes, corr_threshold = 0.25,
                                 hg19ToHg38.over.chain, verbose = FALSE) {
  fantom_variants <- data.frame()
  list.variants <- vector('list', nrow(omim_genes))

  for(gene in 1:nrow(omim_genes)){
    enhancers_df <- get_fantom5_enhancers_from_hgnc(fantom_df = fantom_df,
                                                    hgnc_symbols = omim_genes[gene, "hgnc_symbol"],
                                                    corr_threshold = corr_threshold)
    if(nrow(enhancers_df) != 0) {
      enhancers_df <- GenomicRanges::makeGRangesFromDataFrame(enhancers_df,
                                                              keep.extra.columns = TRUE)
      enhancers_df <- unlist(rtracklayer::liftOver(enhancers_df, rtracklayer::import.chain(hg19ToHg38.over.chain)))
      fantom_locs <- paste0(GenomeInfoDb::seqnames(enhancers_df), ":",
                            BiocGenerics::start(enhancers_df)-1, ":",
                            BiocGenerics::end(enhancers_df))
      fantom_locs <- sub("^chr", "", fantom_locs)

      
      variants_loc <- get_variants_from_locations(fantom_locs,
                                                  verbose = verbose)
      if(length(variants_loc) != 0){
        enhancer_variants <- cbind(ensembl_gene_id = omim_genes[gene, "ensembl_gene_id"],
                                   hgnc_symbol = omim_genes[gene, "hgnc_symbol"],
                                   variants_loc)
  
        enhancer_variants_df <- format_output(chr = unlist(enhancer_variants$seq_region_name),
                                              pos =  unlist(enhancer_variants$start),
                                              rsid = unlist(enhancer_variants$id),
                                              ensembl_gene_id = enhancer_variants$ensembl_gene_id,
                                              hgnc_symbol = enhancer_variants$hgnc_symbol)
  
        list.variants[[gene]] <- enhancer_variants_df
      }
    }
  }
  
  # Removing the NULL elements of the list if exists
  list.variants <- list.variants[!sapply(list.variants, is.null)]
  # Then concatenate the list into a data.frame
  fantom_variants <- do.call('rbind', list.variants)
  if(length(fantom_variants) != 0){
    fantom_variants$source <- "fantom5"
  }

  return(fantom_variants)
}


# ---- GWAS ----

#' @title Create a gwaswloc object
#' @description Create a gwaswloc object from \code{\link[gwascat]{makeCurrentGwascat}}
#' or by reading a local file (if "vargen_dir" specified as a parameter). If 
#' "vargen_dir" contains more than one gwas catalog, the user will be prompted to 
#' choose one.
#' The function will use the filename to determine the extract date, please have
#' it in the format: \[filename\]_r**YYYYY**-**MM**-**DD**.tsv
#' eg: "gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv"
#'
#' @param vargen_dir (optional) a path to vargen data directory, created 
#' during \code{\link{vargen_install}}. If not specified, the gwas object will
#' be created with \code{\link[gwascat]{makeCurrentGwascat}}. If specified, this
#' function will look for files that begin with "gwas_catalog". If more than one 
#' is found, the user will have to choose one via a text menu 
#' @param verbose if true, will print progress information (default: FALSE)
#'
#' @examples
#' gwas_cat <- create_gwas()
#' @export
create_gwas <- function(vargen_dir, verbose = FALSE){
  # If usemakegwas is true (user choice or no local gwas file found, then we 
  # download the catalog with gwascat::makeCurrentGwascat())
  usemakegwas <- FALSE
  
  if(missing(vargen_dir)){
    usemakegwas <- TRUE
  } else{
    # Listing the list of gwas catalogs in "vargen_dir"
    gwasfiles <- list.files(path = vargen_dir, pattern = "gwas_catalog*", 
                            full.names = TRUE, include.dirs = FALSE)
    
    # If more than one catalog found, we let the user choose
    if(length(gwasfiles) > 1){
      print(paste0("More than 1 gwas catalog found in ", vargen_dir))
      print("Please choose one of the file (type 0 to use gwascat::makeCurrentGwascat")
      # Choice can only get indexes belonging to the menu, so no need to check for
      # out of array indexes.
      choice <- utils::menu(choices = gwasfiles)
      if(choice == 0){
        print("No correct choice made, downloading gwas catalog with gwascat::makeCurrentGwascat")
        usemakegwas <- TRUE
      } else {
        gwascat_file <- gwasfiles[choice]
      }
      
    } else if(length(gwasfiles) == 1){
      gwascat_file <- gwasfiles
    } else {
      # Any other option we use "gwascat::makeCurrentGwascat()"
      if(verbose) print(paste0("No gwas catalogs detected in ", vargen_dir, 
                               ". Downloading gwas catalog with gwascat::makeCurrentGwascat"))
      usemakegwas <- TRUE
    }
  }
  
  if(usemakegwas == TRUE){
    require("gwascat")
    utils::data("ebicat38")
    # The following line downloads file and build the gwasloc object all in one step
    gwas_cat <- gwascat::makeCurrentGwascat(table.url = "http://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
                                            fixNonASCII = TRUE,
                                            genome = "GRCh38",
                                            withOnt = TRUE)
    extract_date <- as.character(Sys.Date())
  } else {
    # File from "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
    # quote = "" is important to avoid error about EOF in quoted string.
    gwas_cat <- utils::read.delim(file = gwascat_file, check.names = FALSE, quote = "",
                                  stringsAsFactors = FALSE, header = TRUE, sep = "\t")
    # This suppose the file is under the format: gwas_catalog_v1.0.2-associations_e93_r2019-01-11.tsv
    extract_date <- sub(".*_r", "", gwascat_file)
    extract_date <- sub(".tsv", "", extract_date)
  }
  
  # We transform the gwas cat object into a GRanges object:
  gwas_cat <- gwascat:::gwdf2GRanges(df = gwas_cat, extractDate = extract_date)
  rtracklayer::genome(gwas_cat) <- "GRCh38"

  return(gwas_cat)
}


#' @title List the available gwas traits
#' @description Return the gwas traits available in the gwas catalog, based on keywords.
#' The traits are found using \code{\link[base]{grep}} on the `DISEASE/TRAIT`
#' column from the gwas object produced by \code{\link[gwascat]{makeCurrentGwascat}}.
#' Output can be used as parameter for \code{\link{vargen_pipeline}}
#'
#' @param keywords a vector of keywords to grep the traits from the gwas catalog. (default: "")
#' @param vargen_dir (optional) a path to vargen data directory, created 
#' during \code{\link{vargen_install}}. If not specified, the gwas object will
#' be created with \code{\link[gwascat]{makeCurrentGwascat}}. If specified, this
#' function will look for files that begin with "gwas_catalog". If more than one 
#' is found, the user will have to choose one via a text menu 
#' @return a vector of traits
#'
#' @examples
#' list_gwas_traits(keywords = c("type 1 diabetes", "Obesity"))
#' @export
list_gwas_traits <- function(keywords = "", vargen_dir) {
  traits <- c()
  
  # If the user prefers to create a gwasloc object from gwascat::makeCurrentGwascat()
  # instead of using the local catalog in "vargen_dir" he can use "makegwas = TRUE"
  if(missing(vargen_dir)){
    gwas_cat <- create_gwas()
  } else{
    gwas_cat <- create_gwas(vargen_dir = vargen_dir)
  }

  # Collapsing the values with a "|" for "OR" in grep search
  keys <- paste(keywords, collapse = "|")
  traits <- unique(gwas_cat$`DISEASE/TRAIT`[grep(keys,
                                                 gwas_cat$`DISEASE/TRAIT`,
                                                 ignore.case=T)])
  if(length(traits) == 0) print(paste0("No gwas traits found for '", keys, "'"))

  return(unique(traits))
}


#' @title Get the variants from the gwas catalog associated to the traits of interest
#' @description uses \code{\link[gwascat]{subsetByTraits}} to get the variants.
#'
#' @param gwas_cat output from \code{\link{create_gwas}}
#' @param gwas_traits a vector of gwas traits, can be obtained from
#' \code{\link{list_gwas_traits}}
#'
#' @return a data.frame contaning the variants linked to the traits in the gwas catalog
#' The data.frame contains the following columns:
#' \itemize{
#'   \item chr (chromosome)
#'   \item pos (position of the variant)
#'   \item rsid (variant ID)
#'   \item ensembl_gene_id ("gene id" of the gene associated with the variants)
#'   \item hgnc_symbol ("hgnc symbol" of the gene associated with the variants)
#'   \item source (here "gwas")
#' }
#'
#' @examples
#' gwas_cat <- create_gwas()
#' # if you need to get the list of gwas traits:
#' # list_gwas_traits(keywords = c("Obesity"))
#' obesity_gwas <- c("Obesity (extreme)", "Obesity-related traits", "Obesity")
#' gwas_cat <- create_gwas()
#'
#' gwas_variants <- get_gwas_variants(gwas_cat, obesity_gwas)
#' @export
get_gwas_variants <- function(gwas_cat, gwas_traits){
  require("gwascat")
  gwas_variants <- gwascat::subsetByTraits(x = gwas_cat, tr = gwas_traits)
  GenomeInfoDb::seqlevelsStyle(gwas_variants) <- "UCSC"

  # "gwas_variants@ranges" contains "start" "stop" "width", that is why we only
  # select columns c(1,2,5,6,7,8).
  gwas_variants_df <- unique(data.frame(chr = gwas_variants@seqnames, 
                                        pos = gwas_variants@ranges,
                                        rsid = gwas_variants$SNPS, 
                                        ensembl_gene_id = gwas_variants$SNP_GENE_IDS,
                                        hgnc_symbol = gwas_variants$MAPPED_GENE, 
                                        source = "gwas")[,c(1,2,5,6,7,8)])
  colnames(gwas_variants_df)[2] <- "pos"

  return(gwas_variants_df)
}



#' @title Manhattan plot for variants found in GWAS
#' @description Display a manhattan plot. Only the variants related to the traits
#' given as parameter will be displayed. The two horizontal lines on the plot
#' correspond to the "suggestive" and "significant" threshold in genome wide
#' studies.
#'
#' @param gwas_cat a gwaswloc object obtained from \code{\link[gwascat]{makeCurrentGwascat}}
#' @param traits a vector with the trait of interest (as characters). The list
#' of available traits can be obtained with \code{\link{list_gwas_traits}}
#' @return nothing, diplay the plot
#'
#' @references
#' Lander E, Kruglyak L. Genetic dissection of complex traits: guidelines for
#' interpreting and reporting linkage results. Nat Genet. 1995;11(3):241-7.
#'
#' @examples
#' gwas_cat <- create_gwas()
#'
#' plot_manhattan_gwas(gwas_cat = gwas_cat, traits = c("Type 1 diabetes", "Type 2 diabetes"))
#' @export
plot_manhattan_gwas <- function(gwas_cat, traits) {
  # Check if the gwas traits are in the gwas catalog:
  for(trait in traits){
    if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)) stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
  }

  # Need ggbio for manhattan plot in gwascat
  require("ggbio")
  # These are the two standard gwas thresholds, they will be plotted as lines
  # on the manhattan plot.
  suggestive  <- 1*(10^-5)
  significant <- 5*(10^-8)

  # Just select a subset of all the variants in the gwas catalog (those related
  # to the traits of interest)
  variants_traits <- gwascat::subsetByTraits(x = gwas_cat, tr = traits)

  gwascat::traitsManh(gwr = variants_traits, sel = gwas_cat, traits = traits) +
    # genome-wide significant threshold (p-value < 1 x 10-8)
    # because : -log10(5*10^-8) = 7.30
    ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(suggestive),
                                     linetype = paste0("Suggestive: ",  suggestive)),
                        color = "blue",  size = 0.3) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(significant),
                                     linetype = paste0("Significant: ", significant)),
                        color = "red", size = 0.3) +

    ggplot2::xlab("Genomic Coordinates") + ggplot2::ylab("-log10(p-value)") +

    ggplot2::theme(strip.text.x = ggplot2::element_text(size=6),
                   axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +

    # Overriding the guide legend for the linetype
    # To allow for a legend for "geom_hline"
    ggplot2::scale_linetype_manual(name = "Thresholds p-values\n", values = c(2,2),
                                   guide = ggplot2::guide_legend(override.aes = list(color = c("red", "blue"))))
}


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
      stop(paste0("Tissue '", tissue_keyword,"' not recognized, aborting... (please check the keyword and the GTEx folder)"))
    }

    selected_files <- c(selected_files, as.character(corr_tissues[,2]))
  }

  return(selected_files)
}


#' @title Convert GTEx IDs to GRanges
#' @description GTEx IDs are in the format "chr_pos_ref_alt_build".
#' Transform the IDs to a GRanges object and get the rsids based on the position.
#' If the build is equals to "b37", then liftOver is performed.
#'
#' @param gtex_ids a vector of gtex_ids in the format "chr_pos_ref_alt_build"
#' @param hg19ToHg38.over.chain the chain file to liftOver locations from
#' hg19 to hg38. GTEx reports positions based on hg19.
#' @param verbose if true, will print information about the conversion (default: FALSE)
#'
#' @return a GRanges object, with the variant locations on hg38.
convert_gtex_id_to_granges <- function(gtex_ids, hg19ToHg38.over.chain, verbose = FALSE) {
  gtex_build <- unique(sub(".*_", "", gtex_ids))

  if(length(gtex_build) > 1){
    stop(paste0("All the variants are not in the same build, builds detected: ",
                paste(gtex_build, collapse = ", ")))
  }

  locs <- strsplit(as.character(gtex_ids), "_")
  # This is a variants, so start and stop are the same (even for Indels)
  # Lines from https://f1000researchdata.s3.amazonaws.com/manuscripts/15281/
  # 05161f5a-754b-4c56-acf5-8d3da7ee1583_13577_-_enrico_ferrero_v2.pdf?doi=10.12688/f1000research.13577.2
  # &numberOfBrowsableCollections=17&numberOfBrowsableGateways=23
  gtex_df <- data.frame(chr = sapply(locs, "[",1),
                        start = sapply(locs, "[", 2),
                        stop = sapply(locs, "[", 2))

  gtex_GRanges <- GenomicRanges::makeGRangesFromDataFrame(gtex_df,
                                                          keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(gtex_GRanges) <- "UCSC"

  # If build if hg19, liftOver to hg38
  if(gtex_build == "b37"){
    if(verbose) print("b37 build detected, performing liftOver")

    gtex_GRanges <- unlist(rtracklayer::liftOver(gtex_GRanges,
                                                 rtracklayer::import.chain(hg19ToHg38.over.chain)))
  }

  return(gtex_GRanges)
}


#' @title Convert GTEx IDs to rsid
#' @description GTEx IDs are in the format "chr_pos_ref_alt_build".
#' Transform the IDs to a dataframe and get the rsids based on the position.
#' If the build is equals to "b37", then liftOver is performed.
#' Use \code{\link{convert_gtex_id_to_granges}}
#'
#' @param gtex_ids a vector of gtex_ids in the format "chr_pos_ref_alt_build"
#' @param hg19ToHg38.over.chain the chain file to liftOver locations from
#' hg19 to hg38. GTEx reports positions based on hg19.
#' @param verbose if true, will print information about the conversion (default: FALSE)
#' @return the same as \code{\link{get_variants_from_locations}}
#'
#' @examples
#' vargen_install(install_dir = "./vargen_data/")
#'
#' convert_gtex_to_rsids(c("1_98929_A_G_b37", "10_61803062_T_G_b37"),
#'                         "./vargen_data/hg19ToHg38.over.chain", FALSE)
#'
#' @export
convert_gtex_to_rsids <- function(gtex_ids, hg19ToHg38.over.chain, verbose = FALSE) {
  gtex_GRanges <- convert_gtex_id_to_granges(gtex_ids, hg19ToHg38.over.chain, verbose)

  return(get_variants_from_locations(paste0(gtex_GRanges@seqnames, ":",
                                            gtex_GRanges@ranges, ":",
                                            gtex_GRanges@ranges),
                                     verbose))
}


#' @title Get variants from GTEx linked to the given ensembl genes
#' @description Take as input one or more tissue files from
#'  "GTEx_Analysis_v8_eQTL", as well as a vector of ensembl gene ids.
#'  This function will return the variants that are associated with changes
#'  in the expression of the selected genes in the selected tissues. The
#'  assocations are based on the "signif_variant_gene_pairs" files.
#'  The ensembl ids from the GTEx file will be converted to stable ids, and the
#'  locations from GTEx translated from hg19 to hg38.
#'
#' @param tissue_files a vector containing the name of the "signif_variant_gene_pairs.txt.gz"
#' files. This will be read using \code{\link[base]{gzfile}}. Output from
#' \code{\link{select_gtex_tissues}} can be used.
#' @param omim_genes output from \code{\link{get_omim_genes}}
#' @param hg19ToHg38.over.chain the chain file to liftOver locations from
#' hg19 to hg38. GTEx reports positions based on hg19.
#' @param verbose will be given to subsequent functions to print progress.
#' #' @return a data.frame with information about the variants associated with a
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

# vargen_install(install_dir = "./vargen_data/")
#
# adipose_tissues <- select_gtex_tissues(gtex_dir = "./vargen_data/GTEx_Analysis_v8_eQTL/",
#                                     tissues_query = "adipose")
#
# gene_mart <-  connect_to_gene_ensembl()
# DM1_genes <- get_omim_genes(omim_ids = "222100",
#                            gene_mart = gene_mart)
#
# get_gtex_variants(tissue_files = adipose_tissues,
#                                      omim_genes = DM1_genes,
#                                      hg19ToHg38.over.chain = "./vargen_data/hg19ToHg38.over.chain",
#                                      verbose = TRUE)
get_gtex_variants <- function(tissue_files, omim_genes, hg19ToHg38.over.chain,
                              verbose = FALSE){
  # This will contain the master vector of variants obtained from the different tissues
  gtex_master <- c()

  list.tissues <- vector('list', length(tissue_files))
  i <- 1
  for(file in tissue_files){
    list.tissues[[i]] <- utils::read.delim(gzfile(file), stringsAsFactors = FALSE)
    i <- i + 1
  }
  gtex_master <- do.call('rbind', list.tissues)

  # Tranforming the "ensembl ID" from GTEx to "stable ensembl gene id"
  # eg: ENSG00000135100.17 to ENSG00000135100 (the ".17" correspond to the version number)
  gtex_master$stable_gene_id <- stringr::str_replace(gtex_master$gene_id,
                                                     pattern = ".[0-9]+$",
                                                     replacement = "")

  gtex_variants <- data.frame()
  list.variants <- vector('list', nrow(omim_genes))
  i <- 1
  for(gene in omim_genes$ensembl_gene_id){
    # Get gtex variants that alter the expression of a specific gene
    gene_variants <- gtex_master[gtex_master$stable_gene_id %in% gene, ]

    if(nrow(gene_variants) > 0){
      # GTEx ids are =/= to rsid, we need to convert them to have consistency
      gene_variants <- convert_gtex_to_rsids(gtex_ids = gene_variants$variant_id,
                                             hg19ToHg38.over.chain = hg19ToHg38.over.chain,
                                             verbose = verbose)

      gene_variants_df <- format_output(chr = unlist(gene_variants$seq_region_name),
                                        pos =  unlist(gene_variants$start),
                                        rsid = unlist(gene_variants$id),
                                        ensembl_gene_id = omim_genes[omim_genes$ensembl_gene_id == gene, "ensembl_gene_id"],
                                        hgnc_symbol = omim_genes[omim_genes$ensembl_gene_id == gene, "hgnc_symbol"])

      list.variants[[i]] <- gene_variants_df
      i <- i + 1
    }
  }
  
  # Removing the NULL elements of the list if exists
  list.variants <- list.variants[!sapply(list.variants, is.null)]
  # Then concatenate the list into a data.frame
  gtex_variants <- do.call('rbind', list.variants)

  # If no variants are found, we do not add the source
  if(length(gtex_variants) != 0)  gtex_variants$source <- "gtex"

  return(gtex_variants)
}


# ---- vargen Pipeline ----

#' @title Download the files needed to run \code{\link{vargen_pipeline}}
#' @description Only need to run it once. Will download the following files in
#' "install_dir":
#' \itemize{
#'   \item the latest gwas catalog, eg: gwas_catalog_v1.0.2-associations_e93_r2019-01-11.tsv
#'   \item hg19ToHg38.over.chain.gz (will be unzipped)
#'   \item GTEx_Analysis_v8_eQTL.tar.gz (will be untared)
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
  if(verbose) print("Download GTEx variant association file... This may take a while")
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
#' }
#'
#' @examples
#' vargen_install("./vargen_data/")
#'
#' # Simple query
#' DM1_simple <- vargen_pipeline(vargen_dir = "./vargen_data/", omim_morbid = "222100",
#'                               fantom_corr = 0.25, outdir = "./", verbose = TRUE)
#'
#'
#' # Query with gtex and gwas
#' pancreas_tissues <- select_gtex_tissues(gtex_dir = "./vargen_data/GTEx_Analysis_v8_eQTL/",
#'                                         tissues_query = "pancreas")
#'
#' # list_gwas_traits("diabetes")
#'
#' DM1 <- vargen_pipeline(vargen_dir = "./vargen_data/", omim_morbid = "222100",
#'                        fantom_corr = 0.25, outdir = "./",
#'                        gtex_tissues = pancreas_tissues,
#'                        gwas_traits = "Type 1 diabetes", verbose = TRUE)
#' @export
vargen_pipeline <- function(vargen_dir, omim_morbid_ids, fantom_corr = 0.25, 
                            outdir = "./", gtex_tissues, gwas_traits, verbose = FALSE) {
  if(missing(omim_morbid_ids)){
    stop("Please provide OMIM morbid ids. Stopping now")
  }

  if(!file.exists(outdir)){
    if(verbose) print(paste0("Creating folder '", outdir, "'"))
    dir.create(outdir)
  }

  #_____________________________________________________________________________
  # Loading the necessary resources
  #_____________________________________________________________________________
  if(verbose) print("Connection to gene and snp marts...")
  gene_mart <- connect_to_gene_ensembl()
  snp_mart <- connect_to_snp_ensembl()

  # no gwas traits = no need to generate the gwas object
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

  if(verbose) print(paste0("Reading the liftOver chain file... '",
                           vargen_dir, "/hg19ToHg38.over.chain'"))
  hg19ToHg38.over.chain <- paste0(vargen_dir, "/hg19ToHg38.over.chain")


  #_____________________________________________________________________________
  # Getting variants from genes related to OMIM disease
  #_____________________________________________________________________________
  if(verbose) print("Starting the pipeline...")
  omim_all_genes <- data.frame()
  master_variants <- data.frame()
  for(omim_morbid in omim_morbid_ids){
    
    # First: get all the genes linked to the different omim ids:
    omim_genes <- get_omim_genes(omim_morbid, gene_mart)
    if(nrow(omim_genes) == 0){
      if(verbose) print(paste0("warning: no genes found for omim id: ", omim_morbid))
    } else {
      
      omim_all_genes <- rbind(omim_all_genes, omim_genes)
      
      # We get the variants on the genes:
      genes_variants <- get_omim_variants(omim_genes = omim_genes, 
                                          verbose = verbose)
      
      if(length(genes_variants) != 0) master_variants <- rbind(master_variants, genes_variants)
      
      # We get the variants on the enhancers of the genes:
      fantom_variants <- get_fantom5_variants(fantom_df = fantom_df,
                                              omim_genes = omim_genes,
                                              corr_threshold = fantom_corr,
                                              hg19ToHg38.over.chain = hg19ToHg38.over.chain,
                                              verbose = verbose)

      if(length(fantom_variants) != 0) master_variants <- rbind(master_variants, fantom_variants)
    }
  }
  
  # We write the list of genes in a file.
  utils::write.table(x = omim_all_genes, quote = FALSE, sep = "\t", row.names = FALSE,
                     file = paste0(outdir, "/genes_info.tsv"))

  #_____________________________________________________________________________
  # Getting variants associated with change of expression in GTEx (need tissues as input)
  #_____________________________________________________________________________
  if(!missing(gtex_tissues)){
    gtex_variants <- get_gtex_variants(tissue_files = gtex_tissues,
                                       omim_genes = omim_all_genes,
                                       hg19ToHg38.over.chain = hg19ToHg38.over.chain,
                                       verbose = verbose)

    if(length(gtex_variants) != 0) master_variants <- rbind(master_variants, gtex_variants)

  } else{
    if(verbose) print("No values for 'gtex_tissues', skipping GTEx step...")
  }

  #_____________________________________________________________________________
  # Getting variants associated to the disease in the gwas catalog
  #_____________________________________________________________________________
  # GWAS variants (only if list of gwas traits were given)
  if(!missing(gwas_traits)){
    master_variants <- rbind(master_variants, get_gwas_variants(gwas_cat, gwas_traits))
  } else{
    if(verbose) print("No values for 'gwas_traits', skipping gwas step...")
  }

  # writing the variants data.frame to a file
  if(verbose) print(paste0("Writing the variants to ",
                           paste0(outdir, "/vargen_variants.tsv")))
  utils::write.table(x = master_variants, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
                     file = paste0(outdir, "/vargen_variants.tsv"))

  return(unique(master_variants))
}


#' @title Get variants related to specific genes.
#' @description Alternative to vargen_pipeline for a set of custom genes.
#' The variants are fetched from the following sources:
#' \itemize{
#'   \item get variants on the genes
#'   \item get the variants on the enhancers / promoters of the genes
#'   \item GTEx: get variants impacting the expression of the genes in specific tissues.
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
                          gtex_tissues, gwas_traits, verbose = FALSE) {
  if (!file.exists(outdir)){
    if(verbose) print(paste0("Creating folder '", outdir, "'"))
    dir.create(outdir)
  }

  #_____________________________________________________________________________
  # Loading the necessary resources
  #_____________________________________________________________________________
  if(verbose) print("Connection to gene and snp marts...")
  gene_mart <- connect_to_gene_ensembl()
  snp_mart <- connect_to_snp_ensembl()

  # no gwas traits = no need to generate the gwas object
  if(!missing(gwas_traits)){
    if(verbose) print("Building the gwascat object...")
    #if(missing(gwascat_file)){
    #  gwas_cat <- create_gwas()
    #} else {
    gwas_cat <- create_gwas(vargen_dir)
    #}
    # Check if the gwas traits are in the gwas catalog:
    for(trait in gwas_traits){
      if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)) stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
    }
  }

  if(verbose) print(paste0("Reading the enhancer tss association file for FANTOM5... '" ,
                           vargen_dir, "/enhancer_tss_associations.bed'"))
  fantom_df <- prepare_fantom(enhancer_tss_association = paste0(vargen_dir, "/enhancer_tss_associations.bed"))


  if(verbose) print(paste0("Reading the liftOver chain file... '", vargen_dir, "/hg19ToHg38.over.chain'"))
  hg19ToHg38.over.chain <- paste0(vargen_dir, "/hg19ToHg38.over.chain")

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

  utils::write.table(x = genes_info, quote = FALSE, sep = "\t", row.names = FALSE,
                     file = paste0(outdir, "/custom_genes_info.tsv"))

  # Get the variants on the genes: (we are reusing the get_omim_variants function)
  master_variants <- get_omim_variants(omim_genes = genes_info, verbose = verbose)
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

  if(length(fantom_variants) != 0) master_variants <- rbind(master_variants, fantom_variants)

  #_____________________________________________________________________________
  # Getting variants associated with change of expression in GTEx (need tissues as input)
  #_____________________________________________________________________________
  if(!missing(gtex_tissues)){
    gtex_variants <- get_gtex_variants(tissue_files = gtex_tissues,
                                       omim_genes = genes_info,
                                       hg19ToHg38.over.chain = hg19ToHg38.over.chain,
                                       verbose = verbose)

    if(length(gtex_variants) != 0) master_variants <- rbind(master_variants, gtex_variants)

  } else{
    print("No values for 'gtex_tissues', skipping GTEx step...")
  }


  #_____________________________________________________________________________
  # Getting variants associated to the disease in the gwas catalog
  #_____________________________________________________________________________
  # GWAS variants (only if list of gwas traits were given)
  if(!missing(gwas_traits)){
    master_variants <- rbind(master_variants, get_gwas_variants(gwas_cat, gwas_traits))
  } else{
    if(verbose) print("No values for 'gwas_traits', skipping gwas step...")
  }

  # writing the variants data.frame to a file
  if(verbose) print(paste0("Writing the variants to ", paste0(outdir, "/custom_vargen_variants.tsv")))
  utils::write.table(x = master_variants, append = FALSE, quote = FALSE, sep = "\t",
                     row.names = FALSE, file = paste0(outdir, "/custom_vargen_variants.tsv"))

  return(master_variants)
}
