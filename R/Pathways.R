#---- Utils ----

#' @title Checks if the params given for the kegg_graph function are valid.
#' @description Checks if any given param for othe kegg_graph function are
#' valid and if it isn't, then the .cleanParam is called to fix the paramaters
#' if possible. A list of checks is completed as follows:
#' 1. Check that the minimum information is correct and given and if not
#' given a proper warning to the user. Checks:
#' A. Directory vargen_dir is checked to see if it exists.
#' B. The Gene Mode is checked to see if it is valid and if it isn't then
#' it is made to be valid and set to the default mode of 2.
#' C. Checks that the threshold is an integer if the integer can not be
#' converted to a double from a string then it throws an error.
#' D. Checks that the chromosomes are in the valid format and fixes them
#' if they're able to be corrected.
#' E. Checks if output directory exists and if it doesn't then it creates
#' the directory in the current directory.
#' F. Checks that the rsids are not NA.
#' G. Checks the colors are not null or na.
#' @param: A vector containing all values possible even if not passed into
#' the function by the user.
#' @return If passed: A vector containing all values possible even if not
#' passed into the function by the user and fixed values if something was
#' not right and they could be corrected. If falsed: then is FALSE
.check_kegg_validity <- function(params){
  # Check the directory exists for the VarGen data.
  if(dir.exists(params[[1]]) == FALSE){
    stop(paste0("The directory used for vargen_dir: ", params[[1]],
                " was not found in the current directory ", getwd()))
  }
  
  # Check that the genes are valid.
  if(!(params[[7]] %in% c(0,1,2)) && is.na(params[[7]]) == FALSE) {
    print(paste0("Gene_mode: ", params[[7]],
                 " does not exist so the default mode of 2 is being used instead."))
    params[[7]] <- 2
  }
  
  # Check that the p-value is valid.
  if(typeof(params[[9]]) != "double" && is.na(params[[9]]) == FALSE) {
    if(anyNA(as.double(params[[9]]))){
      print(paste0("p-value threshold could not be converted to a double from the given value: ",
                   params[[9]], ". No pvalue threshold will be used."))
      params[[9]] <- NA
    }
  }
  
  # Check if the chromosome/s are valid and if they are not attempt to fix them.
  if(length(grep("chr[0-9]", x = params[[4]])) != length(params[[4]]) && is.na(params[[4]]) == FALSE){
    #Separate the bad and good chromosomes
    bad_chrs <- params[[4]][!(params[[4]] %in% params[[4]][grep("chr[0-9]", x = params[[4]])])]
    good_chrs <- params[[4]][grep("chr[0-9]", x = params[[4]])]
    
    print(paste0("Chromosomes must be given as a value 'chr3'. Currently given as: ",
                 bad_chrs, ". Cleaning up bad entries to fit entry format as possible."))
    
    for(i in 1:length(bad_chrs)){
      if(anyNA(as.double(bad_chrs[i])) == FALSE) {
        bad_chrs[i] <- paste("chr", bad_chrs[i], sep = "")
      } else {
        bad_chrs <- bad_chrs[bad_chrs!= bad_chrs[i]]
      }
    }
    if(length(good_chrs) != 0){
      params[[4]] <- c(good_chrs, bad_chrs)
      print(paste("These are the remaining chromosomes after fixing: ", params[[4]]))
    } else {
      print(paste0("Chromosomes must be given as a value 'chr3'. Currently given as: ",
                   bad_chrs, ". All chromosomes to be used by default"))
    }
  }
  
  # Check if the output directory exists and if it doesn't then make it.
  if(is.na(params[[2]]) == FALSE) {
    if(dir.exists(params[[2]]) == FALSE) {
      # Make the proper output directory if one isn't given.
      if(.Platform$OS.type == "windows"){
        if(dir.exists(paste0(getwd(), "\\", "KEGG_images\\")) == FALSE){
          dir.create("KEGG_images")
        }
        params[[2]] <- paste0(getwd(), "\\", "KEGG_images\\")
      } else {
        if(dir.exists(paste0(getwd(), "/", "KEGG_images/")) == FALSE){
          dir.create("KEGG_images")
        }
        params[[2]] <- paste0(getwd(), "/", "KEGG_images/")
      }
    }
  }
  
  # If a rsid is given then search just by the rsid and not by any other params even
  # if they have been given. A default title should made for each rsid as well and
  # not just the default.
  if(is.null(params[[8]])) {
    params[[8]] = NA
  } else if(is.na(params[[8]]) == FALSE) {
    traits = NA
    genes = NA
    chrs = NA
  }
  
  # Check that the colors are properly set.
  if(is.null(params[[10]]) || anyNA(params[[10]])){
    params[[10]] = "#da8cde"
  }
  if(is.null(params[[12]]) || anyNA(params[[12]])){
    params[[12]] = "#1254c7"
  }
  if(is.null(params[[11]]) || anyNA(params[[11]])){
    params[[11]] = "#09873e"
  }
  
  return(params)
}

#' @title Restricts variants by for the kegg_graph function.
#' @description Takes the variant list after it has been checked for
#' validity and restricts the GWAS catelog and restricts by the given
#' parameters from the user.
#' @param: A vector containing all values possible even if not passed into
#' the function by the user after they have been deemed valid.
#' @return a dataframe for restricted variants from the GWAS catelog
.restrict_snp_gwas <- function(params){
  if(is.na(params[[1]]) == FALSE && is.na(params[[3]]) == FALSE
     && is.na(params[[4]]) && is.na(params[[6]])) {
    # Case 1: Restrict the variants by trait only.
    gwas_cat <- create_gwas(params[[1]])
    for(trait in params[[3]]){
      if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)){
        stop(paste0("gwas trait '", trait,
                    "' not found in gwas catalog, stopping now."))
      }
    }
    variants_traits <- gwascat::subsetByTraits(x = gwas_cat, tr = params[[3]])
    variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                  gwas_cat))]
  } else if(is.na(params[[1]]) == FALSE && is.na(params[[3]])
            && is.na(params[[4]]) == FALSE && is.na(params[[6]])) {
    # Case 2: Restrict the variants by chromosome only.
    gwas_cat <- create_gwas(params[[1]])
    variants_traits <- gwascat::subsetByChromosome(x = gwas_cat, ch = chrs)
    variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                  gwas_cat))]
  } else if(is.na(params[[1]]) == FALSE && is.na(params[[3]])
            && is.na(params[[4]]) && is.na(params[[6]]) == FALSE) {
    # Case 3: Restrict by gene only.
    gwas_cat <- create_gwas(params[[1]])
    if(params[[7]] == 1){
      variants_traits <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                               gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                               gwas_cat$"REPORTED GENE(S)")
      colnames(variants_traits) <- cbind("rsid", "ensembl_gene_id",
                                         "p-value", "mapped_genes",
                                         "reported_genes")
      
      variants_traits <- data.frame(variants_traits)
      variants_single <- subset(variants_traits,
                                variants_traits[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_traits)){
        print(paste0("searching ", i, " out of ", nrow(variants_traits), " special cases."))
        inner_genes <- strsplit(variants_traits[["mapped_genes"]], ", ")[[i]]
        if(params[[6]] %in% inner_genes){
          variants_single <- merge(x = variants_single, y = variants_traits[i, ],
                                   by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                          "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_traits <- variants_single
      
    } else if(params[[7]] == 0) {
      variants_traits <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                               gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                               gwas_cat$"REPORTED GENE(S)")
      colnames(variants_traits) <- cbind("rsid", "ensembl_gene_id",
                                         "p-value", "mapped_genes",
                                         "reported_genes")
      
      variants_traits <- data.frame(variants_traits)
      variants_multi <- subset(variants_traits,
                               variants_traits[["reported_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_traits)){
        inner_genes <- strsplit(variants_traits[["reported_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_multi <- merge(x = variants_single, y = variants_traits[i, ],
                                  by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                         "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_traits <- variants_multi
    } else {
      variants_traits <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                               gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                               gwas_cat$"REPORTED GENE(S)")
      colnames(variants_traits) <- cbind("rsid", "ensembl_gene_id",
                                         "p-value", "mapped_genes",
                                         "reported_genes")
      
      variants_traits <- data.frame(variants_traits)
      variants_single_reported <- subset(variants_traits,
                                         variants_traits[["reported_genes"]] == params[[6]])
      variants_single_mapped <- subset(variants_traits,
                                       variants_traits[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_traits)){
        print(paste0("searching ", i, " out of ", nrow(variants_traits), " special cases."))
        inner_genes_reported <- strsplit(variants_traits[["reported_genes"]], ", ")[[i]]
        inner_genes_mapped <- strsplit(variants_traits[["mapped_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_single_reported <- merge(x = variants_single_reported, y = variants_traits[i, ],
                                            by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                                   "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_traits <- merge(x = variants_traits, y = variants_single_mapped,
                               by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                      "mapped_genes", "reported_genes"), all = TRUE)
    }
  } else if (is.na(params[[1]]) == FALSE && is.na(params[[3]]) == FALSE
             && is.na(params[[4]]) == FALSE && is.na(params[[6]])) {
    # Case 4: Restrict by both trait and chromosome.
    gwas_cat <- create_gwas(params[[1]])
    for(trait in params[[3]]){
      if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)){
        stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
      }
    }
    variants_traits_chrs <- gwascat::subsetByChromosome(x = gwas_cat, ch = chrs)
    variants_traits_chrs <- variants_traits_chrs[which(IRanges::overlapsAny(variants_traits_chrs,
                                                                            gwas_cat))]
    variants_traits <- gwascat::subsetByTraits(x = variants_traits_chrs, tr = params[[3]])
    variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                  gwas_cat))]
  } else if (is.na(params[[1]]) == FALSE && is.na(params[[3]]) == FALSE
             && is.na(params[[4]]) && is.na(params[[6]]) == FALSE) {
    # Case 5: Restrict by both trait and gene.
    gwas_cat <- create_gwas(params[[1]])
    for(trait in params[[3]]){
      if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)){
        stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
      }
    }
    variants_traits <- gwascat::subsetByTraits(x = gwas_cat, tr = params[[3]])
    variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                  gwas_cat))]
    if(gene_mode == 1){
      variants_info <- cbind(variants_traits$"SNPS", variants_traits$"SNP_GENE_IDS",
                             variants_traits$"P-VALUE", variants_traits$"MAPPED_GENE",
                             variants_traits$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_single <- subset(variants_info,
                                variants_info[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes <- strsplit(variants_info[["mapped_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_single <- merge(x = variants_single, y = variants_info[i, ],
                                   by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                          "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_traits <- variants_single
      
    } else if(gene_mode == 0) {
      variants_info <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                             gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                             gwas_cat$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_multi <- subset(variants_info,
                               variants_info[["reported_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes <- strsplit(variants_info[["reported_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_multi <- merge(x = variants_single, y = variants_info[i, ],
                                  by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                         "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_info <- variants_multi
    } else {
      variants_info <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                             gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                             gwas_cat$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_single_reported <- subset(variants_info,
                                         variants_info[["reported_genes"]] == params[[6]])
      variants_single_mapped <- subset(variants_info,
                                       variants_info[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes_reported <- strsplit(variants_info[["reported_genes"]], ", ")[[i]]
        inner_genes_mapped <- strsplit(variants_info[["mapped_genes"]], ", ")[[i]]
        if(genes %in% inner_genes_mapped || genes %in% inner_genes_reported){
          variants_single_reported <- merge(x = variants_single_reported, y = variants_info[i, ],
                                            by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                                   "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_info <- merge(x = variants_info, y = variants_single_mapped,
                             by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                    "mapped_genes", "reported_genes"), all = TRUE)
    }
  }  else if (is.na(params[[1]]) == FALSE && is.na(params[[3]]) == FALSE
              && is.na(params[[4]]) && is.na(params[[6]]) == FALSE) {
    # Case 6: Restrict by both chromosome and gene.
    gwas_cat <- create_gwas(params[[1]])
    for(trait in params[[3]]){
      if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)){
        stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
      }
    }
    variants_traits <- gwascat::subsetByChromosome(x = gwas_cat, ch = chrs)
    variants_traits <- variants_traits_chrs[which(IRanges::overlapsAny(variants_traits,
                                                                       gwas_cat))]
    if(gene_mode == 1){
      variants_info <- cbind(variants_traits$"SNPS", variants_traits$"SNP_GENE_IDS",
                             variants_traits$"P-VALUE", variants_traits$"MAPPED_GENE",
                             variants_traits$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_single <- subset(variants_info,
                                variants_info[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes <- strsplit(variants_info[["mapped_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_single <- merge(x = variants_single, y = variants_info[i, ],
                                   by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                          "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_traits <- variants_single
      
    }
    else if(gene_mode == 0) {
      variants_info <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                             gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                             gwas_cat$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_multi <- subset(variants_info,
                               variants_info[["reported_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes <- strsplit(variants_info[["reported_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_multi <- merge(x = variants_single, y = variants_info[i, ],
                                  by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                         "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_info <- variants_multi
    } else {
      variants_info <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                             gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                             gwas_cat$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_single_reported <- subset(variants_info,
                                         variants_info[["reported_genes"]] == params[[6]])
      variants_single_mapped <- subset(variants_info,
                                       variants_info[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes_reported <- strsplit(variants_info[["reported_genes"]], ", ")[[i]]
        inner_genes_mapped <- strsplit(variants_info[["mapped_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_single_reported <- merge(x = variants_single_reported, y = variants_info[i, ],
                                            by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                                   "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_info <- merge(x = variants_info, y = variants_single_mapped,
                             by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                    "mapped_genes", "reported_genes"), all = TRUE)
    }
  } else if (is.na(params[[1]]) == FALSE && is.na(params[[3]]) == FALSE
             && is.na(params[[4]]) && is.na(params[[6]]) == FALSE) {
    # Case 6: Restrict by both chromosome and gene.
    gwas_cat <- create_gwas(params[[1]])
    for(trait in params[[3]]){
      if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)){
        stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
      }
    }
    variants_traits <- gwascat::subsetByChromosome(x = gwas_cat, ch = chrs)
    variants_traits <- variants_traits_chrs[which(IRanges::overlapsAny(variants_traits,
                                                                       gwas_cat))]
    if(gene_mode == 1){
      variants_info <- cbind(variants_traits$"SNPS", variants_traits$"SNP_GENE_IDS",
                             variants_traits$"P-VALUE", variants_traits$"MAPPED_GENE",
                             variants_traits$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_single <- subset(variants_info,
                                variants_info[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes <- strsplit(variants_info[["mapped_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_single <- merge(x = variants_single, y = variants_info[i, ],
                                   by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                          "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_traits <- variants_single
      
    }
    else if(gene_mode == 0) {
      variants_info <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                             gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                             gwas_cat$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_multi <- subset(variants_info,
                               variants_info[["reported_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes <- strsplit(variants_info[["reported_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_multi <- merge(x = variants_single, y = variants_info[i, ],
                                  by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                         "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_info <- variants_multi
    } else {
      variants_info <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                             gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                             gwas_cat$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_single_reported <- subset(variants_info,
                                         variants_info[["reported_genes"]] == params[[6]])
      variants_single_mapped <- subset(variants_info,
                                       variants_info[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes_reported <- strsplit(variants_info[["reported_genes"]], ", ")[[i]]
        inner_genes_mapped <- strsplit(variants_info[["mapped_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_single_reported <- merge(x = variants_single_reported, y = variants_info[i, ],
                                            by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                                   "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_info <- merge(x = variants_info, y = variants_single_mapped,
                             by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                    "mapped_genes", "reported_genes"), all = TRUE)
    }
  }  else if (is.na(params[[1]]) == FALSE && is.na(params[[3]]) == FALSE
              && is.na(params[[4]]) && is.na(params[[6]]) == FALSE) {
    # Case 7: Restrict by both chromosome and gene.
    gwas_cat <- create_gwas(params[[1]])
    for(trait in params[[3]]){
      if(!(trait %in% gwas_cat$`DISEASE/TRAIT`)){
        stop(paste0("gwas trait '", trait, "' not found in gwas catalog, stopping now."))
      }
    }
    variants_traits_chrs <- gwascat::subsetByChromosome(x = gwas_cat, ch = chrs)
    variants_traits_chrs <- variants_traits_chrs[which(IRanges::overlapsAny(variants_traits_chrs,
                                                                            gwas_cat))]
    variants_traits <- gwascat::subsetByTraits(x = variants_traits_chrs, tr = params[[3]])
    variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                  gwas_cat))]
    if(gene_mode == 1){
      variants_info <- cbind(variants_traits$"SNPS", variants_traits$"SNP_GENE_IDS",
                             variants_traits$"P-VALUE", variants_traits$"MAPPED_GENE",
                             variants_traits$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_single <- subset(variants_info,
                                variants_info[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes <- strsplit(variants_info[["mapped_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_single <- merge(x = variants_single, y = variants_info[i, ],
                                   by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                          "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_traits <- variants_single
      
    } else if(gene_mode == 0) {
      variants_info <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                             gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                             gwas_cat$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_multi <- subset(variants_info,
                               variants_info[["reported_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes <- strsplit(variants_info[["reported_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_multi <- merge(x = variants_single, y = variants_info[i, ],
                                  by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                         "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_info <- variants_multi
    } else {
      variants_info <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                             gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                             gwas_cat$"REPORTED GENE(S)")
      colnames(variants_info) <- cbind("rsid", "ensembl_gene_id",
                                       "p-value", "mapped_genes",
                                       "reported_genes")
      
      variants_info <- data.frame(variants_info)
      variants_single_reported <- subset(variants_info,
                                         variants_info[["reported_genes"]] == params[[6]])
      variants_single_mapped <- subset(variants_info,
                                       variants_info[["mapped_genes"]] == params[[6]])
      
      for(i in 1:nrow(variants_info)){
        print(paste0("searching ", i, " out of ", nrow(variants_info), " special cases."))
        inner_genes_reported <- strsplit(variants_info[["reported_genes"]], ", ")[[i]]
        inner_genes_mapped <- strsplit(variants_info[["mapped_genes"]], ", ")[[i]]
        if(genes %in% inner_genes){
          variants_single_reported <- merge(x = variants_single_reported, y = variants_info[i, ],
                                            by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                                   "mapped_genes", "reported_genes"), all = TRUE)
        }
      }
      variants_info <- merge(x = variants_info, y = variants_single_mapped,
                             by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                    "mapped_genes", "reported_genes"), all = TRUE)
    }
  } else if(   is.na(params[[1]])) {
    # Case 8: This case should never happen as it should already be caught.
    stop(paste0("No VarGen directory was found."))
  } else {
    # Case 9: No restrictions given.
    gwas_cat <- create_gwas(params[[1]])
    variants_traits <- gwas_cat
  }
}

#' @title Makes output directory for Windows, Linux, or Unix system.
#' @description Makes a directory if none is given or returns the
#' directory of the computer system so that figures can be placed
#' in a given directory.
#' @param: The directory name that should be made.
#' @return a the directory or make a directory is none is given.
.make_dir <- function(params){
  if(is.na(params[[1]])  || params[[1]] == " " || is.double(params[[1]])) {
    # Make the proper output directory if one isn't given.
    if(.Platform$OS.type == "windows") {
      dir.create("KEGG_images")
      output_dir <- paste0(getwd(), "\\", "KEGG_images\\")
    } else {
      dir.create("KEGG_images")
      output_dir <- paste0(getwd(), "/", "KEGG_images/")
    }
  } else {
    if(.Platform$OS.type == "windows") {
      output_dir <- paste0(params[[1]], "/", "KEGG_images/")
    } else {
      output_dir <- paste0(params[[1]], "/", "KEGG_images/")
    }
  }
  return(output_dir)
}

#' @title Get genes from a kegg pathway.
#' @description Gets all the genes from a given pathway.
#' @param: The kegg pathway.
#' @return A list of containing two values, one with
#' the genes from a given pathway and the other with the
#' kegg path ids, and the kegg id.
.make_genes <- function(kegg_pathway){
  # Separate pathway ID from the enzyme/s.
  kegg_path_id <- unlist(strsplit(kegg_pathway, "\\+"))[1]
  #Get the list of numbers, gene symbols and gene description
  kegg_path_id <- paste0("hsa", kegg_path_id)
  #1. Get genes.
  kegg_genes <- tryCatch({
    KEGGREST::keggGet(kegg_path_id)[[1]][["GENE"]]
  }, error = function(e) {
    NULL
  })
  # 2. Check the values and get the id split form the actual genes
  if(is.null(kegg_genes) == FALSE){
    kegg_id <- kegg_genes[seq(1, length(kegg_genes), by = 2)]
    kegg_genes <- kegg_genes[seq(0, length(kegg_genes), 2)]
    kegg_genes <- gsub("\\;.*","",kegg_genes)
    return(list(kegg_genes, gsub("hsa", "map", kegg_path_id), kegg_id))
  } else {
    return(list(kegg_genes, gsub("hsa", "map", kegg_path_id), kegg_id = NULL))
  }
}

#---- Visualizations ----

#' @title Makes the KEGG graphs for a selected variant.
#' @description Allows for the KEGG graphs to be output into a directory, but
#' also allows for KEGG graphs to be generated based upon a seach of a gene/genes,
#' or by the fathmm score threshold. The default name of the image will be the
#' name of the kegg_id unless other wise specfied.  They will be
#' placed into a directory made labeled KEGG_images.
#' @param vargen_dir: The directory that the VarGen information has been
#' download to or is being stored in.
#' @param output_dir: The directory where the files will be put as desired
#' by the user.
#' @param traits: A vector containing the traits that can be seached for. Can be
#' given as either a vector or as a single string.
#' @param chrs: A vector containing the chromosomes that can be seached for.
#' Can be given as either a vector or as a single string. Ex: "chr1"
#' @param title: The title that each of the graphics will be given with a
#' counter added.
#' @param genes: The mapped and reported genes that are found in the file.
#' @param gene_mode: The mode for selecting if the reported genes and the
#' mapped genes should both be search for the given genes or if only one or the
#' other should be.  0 is for mapped genes only, 1 is for only reported, 2 is
#' for both. The defualt is 2.
#' @param pval_thresh: The cut off threshold for the p-values of the variants.
#' @return Nothing; The KEGG pathways figures in a given or default directory.
#' @title Makes the KEGG graphs for a selected variant.
#' @description Allows for the KEGG graphs to be output into a directory, but
#' also allows for KEGG graphs to be generated based upon a seach of a gene/genes,
#' or by the fathmm score threshold. The default name of the image will be the
#' name of the kegg_id unless other wise specfied.  They will be
#' placed into a directory made labeled KEGG_images.
#' @param vargen_dir: The directory that the VarGen information has been
#' download to or is being stored in.
#' @param output_dir: The directory where the files will be put as desired
#' by the user.
#' @param traits: A vector containing the traits that can be seached for. Can be
#' given as either a vector or as a single string.
#' @param chrs: A vector containing the chromosomes that can be seached for.
#' Can be given as either a vector or as a single string. Ex: "chr1"
#' @param title: The title that each of the graphics will be given with a
#' counter added.
#' @param genes: The mapped and reported genes that are found in the file.
#' @param gene_mode: The mode for selecting if the reported genes and the
#' mapped genes should both be search for the given genes or if only one or the
#' other should be.  0 is for mapped genes only, 1 is for only reported, 2 is
#' for both. The defualt is 2.
#' @param pval_thresh: The cut off threshold for the p-values of the variants.
#' @return Nothing; The KEGG pathways figures in a given or default directory.
kegg_graph <- function(vargen_dir, output_dir = NA, traits = NA,
                       chrs = NA, title = NA, genes = NA, gene_mode = 2,
                       rsids = NA, pval_thresh = NA, color = "#da8cde",
                       pval_low_col = "#09873e", pval_high_col = "#1254c7") {
  # 1. Check if the parameters are valid.
  params = .check_kegg_validity(list(vargen_dir, output_dir, traits, chrs, title,
                                     genes, gene_mode, rsids, pval_thresh, color,
                                     pval_low_col, pval_high_col))
  
  # 2. Restrict the gwas catelog by the given paramaters
  variants_traits <- .restrict_snp_gwas(params)
  
  # 3. Get KEGG pathways for the filtered variants above.
  kegg_paths <- biomaRt::getBM(
    attributes = c("kegg_enzyme", "ensembl_gene_id"),
    filters = c("ensembl_gene_id"),
    values = variants_traits$"SNP_GENE_IDS",
    mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
  
  # 4. Make the directory and title as fit given the parameter inputs.
  if(nrow(kegg_paths) == 0 || anyNA(unique(kegg_paths$"kegg_enzyme")[1])){
    stop(paste0("No KEGG pathways found.  Try broading the search."))
  } else {
    variant_info <- cbind(variants_traits$"SNPS", variants_traits$"SNP_GENE_IDS",
                          variants_traits$"P-VALUE", variants_traits$"MAPPED_GENE",
                          variants_traits$"REPORTED GENE(S)")
    colnames(variant_info) <- cbind("rsid", "ensembl_gene_id", "p-value",
                                    "mapped_genes", "reported_genes")
    variant_info <- merge(x = variant_info, y = kegg_paths,
                          by = "ensembl_gene_id", all.x = TRUE)
    
    # Make output directory.
    output_dir <- .make_dir(params[[2]])
    # Check if the title is given.
    titleKey <- FALSE
    if(is.na(params[[5]])){
      titleKey <- TRUE
    }
    
    # 5. Subset and restrict the variants by the pathways and the p-values.
    # In this section the genes are also gone through for each mode (either
    # reported, mapped or both) as well.
    cnt <- 1
    # 5.1 Remove values without a KEGG pathway.
    kegg_paths <- subset(variant_info, is.na(variant_info[["kegg_enzyme"]]) == FALSE)
    kegg_paths <- subset(kegg_paths, kegg_paths[["kegg_enzyme"]] != "")
    
    # 5.2 Get the values that are missing KEGG pathways.
    kegg_paths_na <- subset(variant_info, is.na(variant_info[["kegg_enzyme"]]) == TRUE)
    if(is.na(params[[9]]) == FALSE){
      kegg_paths_na <- subset(kegg_paths_na, kegg_paths_na[["p-value"]] >= is.na(params[[9]]))
    }
    # 5.3 Write output to text file in same directory of values not found with KEGG pathways
    # in the next for loop.
    file_out<-file(paste(output_dir, "variants_without_kegg.txt"))
    
    fg <- rep()
    bg <- rep()
    # 5.4 Restrict it by p-value threshold
    if(is.na(pval_thresh) == FALSE){
      kegg_paths <- subset(kegg_paths, kegg_paths[["p-value"]] >= is.na(params[[9]]))
    }
    
    if(nrow(kegg_paths) == 0){
      writeLines(c(variant_info[["rsid"]][cnt]), file_out)
      print("No KEGG pathways available for this search request. Try broading search if possible.")
    }
    
    # 6. Output KEGG pathways for each pathway in the output firectory
    # assuming it doesn't already exist. Also write the file containing
    # the rsids that were aquired by filtered out due to the p-value
    # threshold.
    index <- c()
    for (kegg_pathway in kegg_paths[["kegg_enzyme"]]) {
      # 6.1 Get the genes from a given pathway.
      kegg_info <- .make_genes(kegg_pathway)
      kegg_genes <- kegg_info[1]
      kegg_path_id <- kegg_info[2]
      # 6.2 Get the color vectors.
      #fg <- replicate(length(unlist(kegg_genes, use.names=FALSE)), "#ffffff")
      #bg <- replicate(length(unlist(kegg_genes, use.names=FALSE)), "#000000")
      # 6.3 Find the unique values that match in the mapped or the reported genes
      # depending on the mode.
      matched <- c()
      if (anyNA(unique(kegg_paths[["mapped_gene"]])) == FALSE) {
        if (params[[7]] == 1) {
          for (i in 1:length(unique(kegg_paths[["mapped_genes"]]))) {
            if (anyNA(match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes))) == FALSE) {
              index <- c(index, match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes)))
              for (j in 1:length(index)) {
                matched <- c(matched, unlist(kegg_genes)[index[j]])
              }
            } else {
              index <- c(index, NA)
              matched <- c(matched, NA)
            }
          }
        } else if (params[[7]] == 0) {
          for (i in 1:length(unique(kegg_paths[["reported_genes"]]))) {
            if (anyNA(match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes))) == FALSE) {
              index <-
                c(index, match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes)))
              for (j in 1:length(index)) {
                matched <- c(matched, unlist(kegg_genes)[index[j]])
              }
            } else {
              index <- c(index, NA)
              matched <- c(matched, NA)
            }
          }
        } else if (params[[7]] == 2) {
          for (i in 1:length(unique(kegg_paths[["mapped_genes"]]))) {
            if (anyNA(match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes))) == FALSE) {
              index <- c(index, match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes)))
              for (j in 1:length(index)) {
                matched <- c(matched, unlist(kegg_genes)[index[j]])
              }
            } else {
              index <- c(index, NA)
              matched <- c(matched, NA)
            }
          }
          for (i in 1:length(unique(kegg_paths[["reported_genes"]]))) {
            if (is.na(match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes))) == FALSE) {
              # This might cause an error... needs testing.
              append(index, c(index, match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes))), length(index))
              for (j in 1:length(index)) {
                matched <- c(matched, unlist(kegg_genes)[index[j]])
              }
            } else {
              index <- c(index, NA)
              matched <- c(matched, NA)
            }
          }
        }
      }
      
      matched <- na.omit(matched)
      print(matched)
      # 6.4 Make sure that some value is found so that a pathway can be made.
      if (is.null(unlist(kegg_genes)) == FALSE
          && identical(unlist(kegg_genes), character(0)) == FALSE
          && length(unlist(kegg_genes, use.names=FALSE)) != 0
          && length(matched) != 0) {
        
        matched <- unlist(unique(matched))
        
        df <- data.frame(unlist(kegg_genes), unlist(kegg_info[3]))
        df <- df[is.element(df$"unlist.kegg_genes.", matched),]
        
        temp <- df[["unlist.kegg_info.3.."]]
        for(i in 1:length(df[["unlist.kegg_info.3.."]])) {
          temp[i] <- paste("hsa:", df[["unlist.kegg_info.3.."]][i], sep = "")
        }
        kegg_id <- temp
        
        kegg_genes <- df[["unlist.kegg_genes."]]
        
        fg <- replicate(length(kegg_genes), color)
        bg <- replicate(length(kegg_genes), "#000000")
        
        url <- KEGGREST::color.pathway.by.objects(paste("path:", kegg_path_id, sep = ""),
                                                  unlist(kegg_genes, use.names=FALSE),
                                                  fg.color.list = unlist(strsplit(fg[!anyNA(fg)], ",")),
                                                  bg.color.list = unlist(strsplit(bg[!anyNA(bg)], ",")))
        
        fg <- c()
        bg <- c()
        # 6.6 Make the correct titles for each of the pathways.
        if (titleKey == TRUE) {
          title <- paste0("hsa_", kegg_path_id, "_kegg.png")
        } else {
          if (file.exists(paste0(output_dir, cnt, "_", title)) == FALSE) {
            title <- paste0(cnt, "_", title)
          }
        }
        # 6.7 Write the output files for the paths not found.
        writeLines(c(kegg_paths[["rsid"]][cnt]), file_out)
        if (file.exists(paste(output_dir, title, sep = ""))) {
          file_cnt <- sapply(output_dir, function(output_dir) {
            length(list.files(output_dir, pattern = title))
          })
          title <- paste0(file_cnt, "_", title, sep = "")
        }
        # 6.8 Download the files from the website.
        download.file(url, paste(output_dir, title, sep = ""), mode = "wb")
        
        cnt <- cnt + 1
      }
    }
    close(file_out)
  }
}
