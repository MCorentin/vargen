# ---- GWAS ----

#' @title Create a gwaswloc object
#' @description Create a gwaswloc object by reading the file at
#' "http://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
#' or by reading a local gwas catalog (if "vargen_dir" specified as a parameter).
#' If "vargen_dir" contains more than one gwas catalog, the user will be prompted to
#' choose one. The function will use the filename to determine the extract date,
#' please have it in the format: \[filename\]_r**YYYYY**-**MM**-**DD**.tsv
#' eg: "gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv"
#'
#' @param vargen_dir (optional) a path to vargen data directory, created
#' during \code{\link{vargen_install}}. If not specified, the gwas object will
#' be created by reading the file at "http://www.ebi.ac.uk/gwas/api/search/downloads/alternative".
#' If specified, this function will look for files that begin with "gwas_catalog".
#' If more than one is found, the user will have to choose one via a text menu
#' @param verbose if true, will print progress information (default: FALSE)
#' @param timeout the timeout set in options(), reading/downloading files online
#' might fail with the default timeout of 60 seconds.
#'
#' @return a gwaswloc object containing the information from the gwas catalog.
#'
#' @examples
#' gwas_cat <- create_gwas()
#' @export
create_gwas <- function(vargen_dir, verbose = FALSE, timeout = 1000){
  # If useURL is TRUE (user choice or no local gwas file found, then we download
  # the catalog from the URL)
  original_timeout <- getOption("timeout")
  options(timeout = timeout)
  if(verbose) print(paste0("Setting the timeout to value '", timeout, "'"))

  useURL <- FALSE
  table.url <- "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
  #table.url <- "ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations.tsv"

  if(missing(vargen_dir)){
    # If the vargen_dir is not specified, we use the url instead
    useURL <- TRUE
  } else{
    # Listing the list of gwas catalogs in "vargen_dir"
    gwasfiles <- list.files(path = vargen_dir, pattern = "gwas_catalog*",
                            full.names = TRUE, include.dirs = FALSE)

    # If more than one catalog found, we let the user choose
    if(length(gwasfiles) > 1){
      print(paste0("More than 1 gwas catalog found in ", vargen_dir))
      print(paste0("Please choose one of the file (type 0 to read the latest gwas catalog from: ", table.url))
      # Choice can only get indexes belonging to the menu (no need to check out of array)
      choice <- utils::menu(choices = gwasfiles)
      if(choice == 0){
        print(paste0("No correct choice made, reading the gwas catalog from: ", table.url))
        useURL <- TRUE
      } else {
        gwascat_file <- gwasfiles[choice]
      }
    } else if(length(gwasfiles) == 1) {
      # If there is only one gwas catalog in the directory, we use it directly
      gwascat_file <- gwasfiles
    } else {
      # Any other option we get the gwas cat from the URL
      if(verbose) print(paste0("No gwas catalogs detected in ", vargen_dir,
                               ". Downloading the gwas catalog..."))
      useURL <- TRUE
    }
  }

  if(useURL == TRUE){
    gwas_cat <- utils::read.delim(file = url(table.url), check.names = FALSE,
                                  quote = "", stringsAsFactors = FALSE,
                                  header = TRUE, sep = "\t")
    extract_date <- as.character(Sys.Date())
  } else {
    # File from "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
    # quote = "" is important to avoid error about EOF in quoted string.
    gwas_cat <- utils::read.delim(file = gwascat_file, check.names = FALSE,
                                  quote = "", stringsAsFactors = FALSE,
                                  header = TRUE, sep = "\t")
    # This suppose the file is under the format:
    # gwas_catalog_v1.0.2-associations_e93_r2019-01-11.tsv
    extract_date <- sub(".*_r", "", gwascat_file)
    extract_date <- sub(".tsv", "", extract_date)
  }

  # We transform the gwas cat object into a GRanges object:
  gwas_cat <- gwascat:::gwdf2GRanges(df = gwas_cat, extractDate = extract_date)
  if(typeof(gwas_cat) == "list"){ gwas_cat <- gwas_cat$okrngs }
  rtracklayer::genome(gwas_cat) <- "GRCh38"

  options(timeout = original_timeout)
  if(verbose) print(paste0("Resetting the timeout to previous value '", original_timeout, "'"))


  return(gwas_cat)
}


#' @title List the available gwas traits
#' @description Return the gwas traits available in the gwas catalog, based on keywords.
#' The traits are found using \code{\link[base]{grep}} on the `DISEASE/TRAIT`
#' column from the gwas object produced by \code{\link{create_gwas}}.
#' Output can be used as parameter for \code{\link{vargen_pipeline}}
#'
#' @param keywords a vector of keywords to grep the traits from the gwas catalog. (default: "")
#' @param gwas_cat output from \code{\link{create_gwas}}
#' @return a vector of traits
#'
#' @examples
#' list_gwas_traits(keywords = c("type 1 diabetes", "Obesity"))
#' @export
list_gwas_traits <- function(keywords = "", gwas_cat) {

  if(missing(gwas_cat) || class(gwas_cat) != 'gwaswloc'){
    gwas_cat <- create_gwas()
    warning("gwas_cat not provided (or not a gwaswloc object).",
            "Created one with create_gwas().")
  }

  traits <- c()

  # Collapsing the values with a "|" for "OR" in grep search
  keys <- paste(keywords, collapse = "|")
  traits <- unique(gwas_cat$`DISEASE/TRAIT`[grep(keys,
                                                 gwas_cat$`DISEASE/TRAIT`,
                                                 ignore.case = TRUE)])
  if(length(traits) == 0) print(paste0("No gwas traits found for '", keys, "'"))

  return(unique(traits))
}


#' @title Get the variants from the gwas catalog associated to the traits of interest
#' @description uses \code{\link[gwascat]{subsetByTraits}} to get the variants.
#'
#' @param gwas_traits a vector of gwas traits, can be obtained from
#' \code{\link{list_gwas_traits}}
#' @param gwas_cat output from \code{\link{create_gwas}}
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
#' # if you need to get the list of gwas traits:
#' # list_gwas_traits(keywords = c("Obesity"))
#' obesity_gwas <- c("Obesity (extreme)", "Obesity-related traits", "Obesity")
#' gwas_cat <- create_gwas()
#'
#' gwas_variants <- get_gwas_variants(obesity_gwas, gwas_cat)
#' @export
get_gwas_variants <- function(gwas_traits, gwas_cat){

  if(missing(gwas_cat) || class(gwas_cat) != 'gwaswloc'){
    gwas_cat <- create_gwas()
    warning("gwas_cat not provided (or not a gwaswloc object).",
            "Created one with create_gwas().")
  }

  gwas_variants <- gwascat::subsetByTraits(x = gwas_cat, tr = gwas_traits)
  GenomeInfoDb::seqlevelsStyle(gwas_variants) <- "UCSC"

  # "gwas_variants@ranges" contains "start" "stop" "width", we only take the start
  gwas_variants_df <- unique(data.frame(chr = gwas_variants@seqnames,
                                        pos = gwas_variants@ranges@start,
                                        rsid = gwas_variants$SNPS,
                                        ensembl_gene_id = gwas_variants$SNP_GENE_IDS,
                                        hgnc_symbol = gwas_variants$MAPPED_GENE,
                                        source = "gwas",
                                        trait = gwas_variants$`DISEASE/TRAIT`))

  return(gwas_variants_df)
}



#' @title Manhattan plot for variants found in GWAS
#' @description Display a manhattan plot. Only the variants related to the traits
#' given as parameter will be displayed. The two horizontal lines on the plot
#' correspond to the "suggestive" and "significant" thresholds in genome wide
#' studies.
#'
#' @param traits a vector with the trait of interest (as characters). The list
#' of available traits can be obtained with \code{\link{list_gwas_traits}}
#' @param gwas_cat a gwaswloc object obtained from \code{\link{create_gwas}}
#' @param list_chr (optional) a list of chromosome to plot. (format: "chr{1-22-X-Y})
#' @return nothing (just display the plot)
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
plot_manhattan_gwas <- function(traits, gwas_cat, list_chr) {

  if(missing(gwas_cat) || class(gwas_cat) != 'gwaswloc'){
    gwas_cat <- create_gwas()
    warning("gwas_cat not provided (or not a gwaswloc object).",
            "Created one with create_gwas().")
  }

  # Check if the gwas traits are in the gwas catalog:
  for(trait in traits){
    if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)) stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
  }

  if(!missing(list_chr)){
    # Checking if the list of chromosomes has a valid format
    valid_chr <- unique(gwas_cat@seqnames@values)
    for(chr in list_chr){
      if(!(chr %in% valid_chr)){
        stop(paste0("Chromosome name: '", chr,
                    "' not found in gwas catalog, stopping now. ",
                    "Valid chromosome names are: ",
                    paste0(valid_chr, collapse = ", ")))
      }
    }
  }

  # Just select the variants  related to the traits of interest
  variants_traits <- gwascat::subsetByTraits(x = gwas_cat, tr = traits)

  # Select the variants on the selected chromosomes
  if(!missing(list_chr)){
    # We check if the chromsome has variants for the traits
    if(length(list_chr[list_chr %in% variants_traits@seqnames@values]) > 0){
      variants_traits <- gwascat::subsetByChromosome(x = variants_traits, ch = list_chr)
    } else {
      stop(paste0("No variants found on chromosomes: ", paste0(list_chr, collapse = ", "),
                  " for the traits: ", paste0(traits, collapse = ", ")))
    }
  }

  # These are the two standard gwas thresholds, they will be plotted as lines
  # on the manhattan plot.
  suggestive  <- 1*(10^-5)
  significant <- 5*(10^-8)

  # Chromosomes names and lengths to have correct lengths in the manhattan plot
  chrNames <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                "chrX", "chrY")
  chrNames <- as.factor(chrNames)

  chrLengths <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                  159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                  114364328, 107043718, 101991189, 90338345, 83257441,  80373285,
                  58617616,  64444167, 46709983,  50818468, 156040895,  57227415)

  names(chrLengths) <- chrNames
  GenomeInfoDb::seqlengths(variants_traits) <- chrLengths[names(GenomeInfoDb::seqlengths(variants_traits))]

  # The next 9 commands below are from the traitsManh function from the gwascat package.
  # I do not use the function directly because it throws an error if ggplot is
  # required after ggbio. (because ggplot2::autoplot will be called instead of
  # ggbio::autoplot, and ggplot2::autoplot does not handle gwaswloc objects).
  #   gwascat::traitsManh(gwr = variants_traits, selr = gwas_cat, traits = traits)
  variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                gwas_cat))]

  availtr <- as.character(variants_traits$`DISEASE/TRAIT`)
  oth <- which(!(availtr %in% traits))
  availtr[oth] <- "Other"
  variants_traits$Trait <- availtr
  pv <- variants_traits$PVALUE_MLOG
  # Capping the values at 25
  variants_traits$PVALUE_MLOG = ifelse(pv > 25, 25, pv)


  # Generate the plot
  ggbio::autoplot(object = variants_traits, geom = "point",
                  ggplot2::aes(y = variants_traits$PVALUE_MLOG,
                               color = variants_traits$Trait)) +
    # To add the chromosome name on top
    ggplot2::facet_grid(. ~factor(seqnames, levels = chrNames), scales = "free_x", switch="both") +

    # genome-wide significant threshold (p-value < 1 x 10-8)
    # because : -log10(5*10^-8) = 7.30
    ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(suggestive),
                                     linetype = paste0("Suggestive: ",  suggestive)),
                        color = "blue",  size = 0.3) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(significant),
                                     linetype = paste0("Significant: ", significant)),
                        color = "red", size = 0.3) +

    ggplot2::xlab("Genomic Coordinates") + ggplot2::ylab("-log10(p-value)") +

    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10),
                   axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +

    # Overriding the guide legend for the linetype to allow for geom_hline's legend
    ggplot2::scale_linetype_manual(name = "Thresholds p-values\n", values = c(2,2),
                                   guide = ggplot2::guide_legend(override.aes = list(color = c("red", "blue"))))
}
