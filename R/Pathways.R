#---- Universal Utils ----

#' @title Checks Dataframe.
#' @description Checks if a dataframe is valid and if not then stops the
#' program. This was meant for the VarGen pipeline dataset but can be
#' used for any dataframe.
#' @param data: Meant for data given from the VarGen pipeline.
#' @return None; stops the program if the dataframe is not valid.
.check_dataframe <- function(data) {
  # Check if the dataset is null, na, or contains no data.
  if(is.na(data) || is.null(data) || nrow(data) == 0) {
    stop(paste0("The given data from the VarGen pipeline is either null,",
                "na, or contains no data."), sep = "")
  }
}

#' @title Makes output directory for Windows, Linux, or Unix system.
#' @description Checks if a directory is valid and if it exists.
#' If it doesn't exist, then it makes the directory. The function works
#' on Linux, Mac, and windows OS.
#' @example output_dir = ./this_is_an_exmaple_dir/
#' @example output_dir = /this_is_an_exmaple_dir
#' @example output_dir = /this_is_an_exmaple_dir/to_my_nonExisting_dir/
#' @example this will cause an error: output_dir = help_is_an_exmaple_dir/
#' @param output_dir: The user given directory from function input.
#' @param name: The name of the function using this util function.
#' @return The output directory.
.check_output_dir <- function(output_dir, name = NULL) {
  if(is.na(output_dir) || is.null(output_dir)) {
    if(.Platform$OS.type == "windows") {
      if(dir.exists(paste0(getwd(), "\\", name, "\\")) == FALSE) {
        output_dir <- paste0(".", output_dir)
        dir.create(output_dir)
        output_dir <- gsub(".", "", output_dir)
      }
      output_dir <- paste0(getwd(), "\\", output_dir, "\\", sep = "")
    } else {
      if(dir.exists(paste0(getwd(), "/", name, "/")) == FALSE){
        output_dir <- paste0(".", output_dir)
        dir.create(output_dir)
        output_dir <- gsub(".", "", output_dir)
      }
      output_dir <- paste0(getwd(), "/", output_dir, "/", sep = "")
      print(paste0("Using output directory ", output_dir, "/"))
    }
  } else {
    output_dir <- gsub("[.\\/]+", "VARGENrEPLACEsTRING", output_dir)
    if(.Platform$OS.type == "windows") {
      output_dir <- gsub("VARGENrEPLACEsTRING", "[\\]", output_dir)
      if(dir.exists(paste0(getwd(), output_dir)) == FALSE) {
        output_dir <- paste0(".", output_dir)
        dir.create(output_dir)
        output_dir <- gsub(".", "", output_dir)
      }
      output_dir <- paste0(getwd(), output_dir, "\\", sep = "")
    } else {
      output_dir <- gsub("VARGENrEPLACEsTRING", "/", output_dir)
      if(dir.exists(paste0(getwd(), output_dir)) == FALSE){
        output_dir <- paste0(".", output_dir)
        dir.create(output_dir)
        output_dir <- gsub(".", "", output_dir)
      }
      output_dir <- paste0(getwd(), output_dir, "/", sep = "")
      print(paste0("Using output directory ", output_dir))
    }
  }
  
  return(output_dir)
}


#---- KEGG Utils ----

.check_vargen_dir <- function(vargen_dir) {
  # 1. Check the directory exists for the VarGen data.
  if(dir.exists(vargen_dir) == FALSE) {
    stop(paste0("The directory used for vargen_dir: ", vargen_dir,
                " was not found in the current directory ", getwd()))
  }
}

.check_gene_mode <- function(gene_mode, genes) {
  # 2. Check that the gene mode is valid.
  if((!(gene_mode %in% c(0,1,2)) && is.na(gene_mode) == FALSE && is.na(gene_mode)) ||
     (!(gene_mode %in% c(0,1,2)) && is.na(gene_mode) == FALSE && is.na(gene_mode) == FALSE)) {
    # Case 1: If the given gene mode is not in one of the valid selection types.
    print(paste0("Gene_mode: ", gene_mode,
                 " does not exist so the default mode of 2 is being used instead."))
    gene_mode <- 2
  } else if(gene_mode %in% c(0,1) && is.na(gene_mode) == FALSE && is.na(genes)) {
    # Case 2: If the given gene mode is in one of the valid selection types, but no genes are given.
    # Only 1 and 0 are selected here as the default is 2 in the program.
    print("No genes were given with gene mode, thus no gene mode can/will be used")
  }
  
  return(gene_mode)
}

.check_pval <- function(pval) {
  # 3. Check that the p-value is valid.
  if(typeof(pval) != "double" && is.na(pval) == FALSE) {
    if(is.na(as.double(pval))){
      print(paste0("p-value threshold could not be converted to a double from the given value: '",
                   pval, "'. No p-value threshold will be used."))
      pval <- NA
    }
  }
  
  return(pval)
}

.check_chrs <-function(chrs) {
  # 4. Check if the chromosome/s are valid and if they are not attempt to fix them.
  if(length(grep("chr[0-9]", x = chrs)) != length(chrs) && is.na(chrs) == FALSE){
    #Separate the bad and good chromosomes
    bad_chrs <- chrs[!(chrs %in% chrs[grep("chr[0-9]", x = chrs)])]
    good_chrs <- chrs[grep("chr[0-9]", x = chrs)]
    
    print(paste0("Chromosomes must be given as a value 'chr3'. Currently given badly as: ",
                 bad_chrs, ". Cleaning up bad entries to fit entry format as possible."))
    
    for(i in 1:length(bad_chrs)){
      if(anyNA(as.double(bad_chrs[i])) == FALSE) {
        bad_chrs[i] <- paste("chr", bad_chrs[i], sep = "")
      } else {
        bad_chrs <- bad_chrs[bad_chrs!= bad_chrs[i]]
      }
    }
    if(length(good_chrs) != 0){
      chrs <- c(good_chrs, bad_chrs)
      out_print <- paste(chrs, collapse = '')
      out_print <- gsub("(\\d+)\\c", "\\1\\ \\c", out_print)
      print(paste0("These are the remaining chromosomes after fixing: ", out_print))
    } else {
      print(paste0("Chromosomes must be given as a value 'chr3'. Currently given badly as: ",
                   bad_chrs, ". All chromosomes to be used by default"))
    }
  }
  
  return(chrs)
}

.check_rsid <-function(rsid, genes, traits, chrs) {
  # 6. If a rsid is given then search just by the rsid and not by any other params even
  # if they have been given. A default title should made for each rsid as well and
  # not just the default; This is done later on.
  if(is.null(rsid)) {
    rsid = NA
  } else if(is.na(rsid) == FALSE) {
    genes = NA
    traits = NA
    chrs = NA
  }
  
  return(list(rsid, genes, traits, chrs))
}

.check_omim <- function(omim, pval_thresh, pval_thresh_show) {
  # 7. Check the OMIM IDs
  if(is.null(omim)) {
    params[[14]] = NA
  } else if(is.null(omim) == FALSE && is.null(pval_thresh == FALSE)) {
    print("No p-value threshold restriction available for OMIM searching.")
    pval_thresh = NA
    pval_thresh_show = 1
  }
  
  return(list(omim, pval_thresh, pval_thresh_show))
}

.check_colors <- function(color, pval_low_col, pval_high_col) {
  # 8. Check that the colors are properly set.
  if(is.null(color) || is.na(color)){
    color = "#da8cde"
  }
  if(is.null(pval_high_col) || is.na(pval_high_col)){
    pval_high_col = "#1254c7"
  }
  if(is.null(pval_low_col) || is.na(pval_low_col)){
    pval_low_col = "#09873e"
  }
  
  return(list(color, pval_low_col, pval_high_col))
}

.check_show_pval_thresh <- function(show_pval_thresh) {
  # 9. Check that the show_pval_thresh is a boolean
  if(show_pval_thresh != 1 && show_pval_thresh != 0) {
    print(paste0("parameter: ", show_pval_thresh, " is not valid. Set to 1"))
    show_pval_thresh = 1
  }
  
  return(show_pval_thresh)
}

.check_trait <- function(trait) {
  # 10. Check that traits not given as null and if so change to be NA.
  if(is.null(trait)) {
    trait = NA
  }
  
  return(trait)
}

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
.kegg_validity_check <- function(params) {
  print("Checking traits to see if they are valid")
  # 1. Check the directory exists for the VarGen data.
  .check_vargen_dir(params[[1]])
  # 2. Check that the gene mode is valid.
  params[[7]] <- .check_gene_mode(params[[7]], params[[6]])
  # 3. Check that the p-value is valid.
  params[[9]] <- .check_pval(params[[9]])
  # 4. Check if the chromosome/s are valid and if they are not attempt to fix them.
  params[[4]] <- .check_chrs(params[[4]])
  # 5. Check if the output directory exists and if it doesn't then make it.
  params[[2]] <- .check_output_dir(params[[2]], "KEGG_images")
  # 6. If a rsid is given then search just by the rsid and not by any other params even
  # if they have been given. A default title should made for each rsid as well and
  # not just the default; This is done later on.
  out <- .check_rsid(params[[8]], params[[3]], params[[6]], params[[4]])
  params[[8]] <- out[[1]]
  params[[3]] <- out[[2]]
  params[[6]] <- out[[3]]
  params[[4]] <- out[[4]]
  # 7. Check the OMIM IDs
  out <- .check_omim(params[[14]], params[[9]], params[[13]])
  params[[14]] <- out[[1]]
  params[[9]] <- out[[2]]
  params[[13]] <- out[[3]]
  # 8. Check that the colors are properly set.
  out <- .check_colors(params[[10]], params[[12]], params[[11]])
  params[[10]] <- out[[1]]
  params[[12]] <- out[[2]]
  params[[11]] <- out[[3]]
  # 9. Check that the show_pval_thresh is a boolean
  params[[13]] <- .check_show_pval_thresh(params[[13]])
  # 10. Check that traits not given as null and if so change to be NA.
  params[[3]] <- .check_trait(params[[3]])
  
  return(params)
}

#' @title Get genes from a kegg pathway.
#' @description Gets all the genes from a given pathway.
#' @param kegg pathway: The kegg pathway.
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
  
  ko_ids <- c()
  # 2. Check the values and get the id split form the actual genes
  if(is.null(kegg_genes) == FALSE){
    kegg_id <- kegg_genes[seq(1, length(kegg_genes), by = 2)]
    kegg_info <- kegg_genes[seq(0, length(kegg_genes), 2)]
    
    kegg_info <- sub(";.*\\[KO:", ",", kegg_info)
    kegg_info <- sub(" \\[EC:.*", "", kegg_info)
    kegg_info <- sub("]", "", kegg_info)
    kegg_info <- strsplit(kegg_info, ",")
    
    # Two column dataset with matching ko ids and genes.
    kegg_info <- t(data.frame(kegg_info))
    kegg_info <- data.frame(kegg_info)
    
    return(list(kegg_info[["X1"]],
                gsub("hsa", "map", kegg_path_id),
                kegg_id, kegg_info[["X2"]]))
  } else {
    return(list(kegg_genes, gsub("hsa", "map", kegg_path_id), kegg_id = NULL, ko_ids))
  }
}

#---- KEGG Functions ----

#' @title Restricts variants by for the kegg_graph function.
#' @description Takes the variant list after it has been checked for
#' validity and restricts the GWAS catelog and restricts by the given
#' parameters from the user.
#' @param params: A vector containing all values possible even if not passed into
#' the function by the user after they have been deemed valid.
#' @return a dataframe for restricted variants from the GWAS catalog
.restrict_snp_gwas <- function(params){
  print("Restricting GWAS catalog by given parameters.")
  # If there is a given rsid or list of then that is the only thing that needs to
  # be looked for. Otherwise, the other restirction paramters can be used.
  if(is.na(params[[8]])) {
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
      
    } else if(is.na(params[[1]]) == FALSE && is.na(params[[3]])
              && is.na(params[[4]]) == FALSE && is.na(params[[6]])) {
      # Case 2: Restrict the variants by chromosome only.
      gwas_cat <- create_gwas(params[[1]])
      variants_traits <- gwascat::subsetByChromosome(x = gwas_cat, ch = params[[4]])
      variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                    gwas_cat))]
    } else if(is.na(params[[1]]) == FALSE && is.na(params[[3]])
              && is.na(params[[4]]) && is.na(params[[6]]) == FALSE) {
      # Case 3: Restrict by gene only.
      gwas_cat <- create_gwas(params[[1]])
      if(params[[7]] == 1) {
        variants_traits <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                                 gwas_cat$"P-VALUE", gwas_cat$"MAPPED_GENE",
                                 gwas_cat$"REPORTED GENE(S)")
        colnames(variants_traits) <- cbind("rsid", "ensembl_gene_id",
                                           "p-value", "mapped_genes",
                                           "reported_genes")
        
        variants_traits <- data.frame(variants_traits)
        variants_single <- subset(variants_traits,
                                  variants_traits[["mapped_genes"]] == params[[6]])
        
        for(i in 1:nrow(variants_traits)) {
          print(paste0("searching ", i, " out of ", nrow(variants_traits), " special cases."))
          inner_genes <- strsplit(variants_traits[["mapped_genes"]], ", ")[[i]]
          if(params[[6]] %in% inner_genes) {
            variants_single <- merge(x = variants_single, y = variants_traits[i, ],
                                     by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                            "mapped_genes", "reported_genes"), all = TRUE)
          }
        }
        variants_traits <- variants_singlef
        
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
          if(params[[6]] %in% inner_genes){
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
          if(params[[6]] %in% inner_genes_reported || params[[6]] %in% inner_genes_mapped){
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
      variants_traits_chrs <- gwascat::subsetByChromosome(x = gwas_cat, ch = params[[4]])
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
          stop(paste0("gwas trait '", trait,
                      "' not found in gwas catalog, stopping now."))
        }
      }
      variants_traits <- gwascat::subsetByTraits(x = gwas_cat, tr = params[[3]])
      variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                    gwas_cat))]
      if(params[[7]] == 1){
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
          if(params[[6]] %in% inner_genes){
            variants_single <- merge(x = variants_single, y = variants_info[i, ],
                                     by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                            "mapped_genes", "reported_genes"), all = TRUE)
          }
        }
        variants_traits <- variants_single
        
      } else if(params[[7]] == 0) {
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
          if(params[[6]] %in% inner_genes){
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
          if(params[[6]] %in% inner_genes_mapped || params[[6]] %in% inner_genes_reported){
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
      variants_traits <- gwascat::subsetByChromosome(x = gwas_cat, ch = params[[4]])
      variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                    gwas_cat))]
      if(params[[7]] == 1){
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
          if(params[[6]] %in% inner_genes){
            variants_single <- merge(x = variants_single, y = variants_info[i, ],
                                     by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                            "mapped_genes", "reported_genes"), all = TRUE)
          }
        }
        variants_traits <- variants_single
        
      }
      else if(params[[7]] == 0) {
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
          if(params[[6]] %in% inner_genes){
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
          if(params[[6]] %in% inner_genes){
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
      variants_traits <- gwascat::subsetByChromosome(x = gwas_cat, ch = params[[4]])
      variants_traits <- variants_traits[which(IRanges::overlapsAny(variants_traits,
                                                                    gwas_cat))]
      if(params[[7]] == 1){
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
          if(params[[6]] %in% inner_genes){
            variants_single <- merge(x = variants_single, y = variants_info[i, ],
                                     by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                            "mapped_genes", "reported_genes"), all = TRUE)
          }
        }
        variants_traits <- variants_single
        
      }
      else if(params[[7]] == 0) {
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
          if(params[[6]] %in% inner_genes){
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
          if(params[[6]] %in% inner_genes){
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
      # Case 7: Restrict by both chromosome and trait.
      gwas_cat <- create_gwas(params[[1]])
      # Holds the indeces for the hits of variants.
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
      variants_traits_chrs <- gwascat::subsetByChromosome(x = variants_traits, ch = params[[4]])
      variants_traits_chrs <- variants_traits_chrs[which(IRanges::overlapsAny(variants_traits_chrs,
                                                                              gwas_cat))]
      
      if(params[[7]] == 1){
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
          if(params[[6]] %in% inner_genes){
            variants_single <- merge(x = variants_single, y = variants_info[i, ],
                                     by = c("rsid", "ensembl_gene_id", "p.value", "reported_genes",
                                            "mapped_genes", "reported_genes"), all = TRUE)
          }
        }
        variants_traits <- variants_single
        
      } else if(params[[7]] == 0) {
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
          if(params[[6]] %in% inner_genes){
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
          if(params[[6]] %in% inner_genes){
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
  } else {
    if(is.na(params[[1]]) == FALSE) {
      gwas_cat <- create_gwas(params[[1]])
      for(rsid in params[[8]]) {
        if(!(rsid %in% gwas_cat$`SNPS`)) {
          stop(paste0("gwas snp '", rsid, "' not found in gwas catalog, stopping now."))
        }
      }
      variant_info <- cbind(gwas_cat$"SNPS", gwas_cat$"SNP_GENE_IDS",
                            gwas_cat$"DISEASE/TRAIT", gwas_cat$"P-VALUE",
                            gwas_cat$"MAPPED_GENE", gwas_cat$"REPORTED GENE(S)")
      colnames(variant_info) <- cbind("rsid", "ensembl_gene_id", "disease", "p_value",
                                      "mapped_genes", "reported_genes")
      variant_info <- data.frame(variant_info)
      variant_info <- subset(variant_info, variant_info[["rsid"]] == params[[8]])
    } else {
      stop(paste("No directory to vargen_data given"))
    }
  }
}

#' @title Restricts variants by for the kegg_graph function.
#' @description Takes the variant list after it has been checked for
#' validity and restricts the OMIM database and restricts by the given
#' parameters from the user.
#' @param params: A vector containing all values possible even if not passed into
#' the function by the user after they have been deemed valid.
#' @param omim_info: An empty vector that will returned with the dataset.
#' @return a dataframe for restricted variants from the OMIM catelog
.restrict_snp_omim <- function(omim_info, params) {
  if(is.na(params[[14]]) == FALSE && is.na(params[[3]]) == FALSE) {
    omim_ids <- list_omim_accessions(connect_to_gene_ensembl(), params[[3]])[["mim_morbid_accession"]]
    if(length(omim_ids) != 0) {
      omim_genes <- get_omim_genes(omim_ids, connect_to_gene_ensembl())
      omim_info <- omim_genes
      #Start the restiction from the other processes.
      if(nrow(omim_info) != 0) {
        if(is.na(params[[8]])){
          if(is.na(params[[3]]) == FALSE && is.na(params[[4]]) && is.na(params[[6]])) {
            return(omim_info)
          } else if(is.na(params[[3]]) == FALSE && is.na(params[[4]]) == FALSE
                    && is.na(params[[6]])) {
            # Restrict by only chromosomes.
            omim_info <- subset(omim_info, omim_info[["chromosome_name"]] == params[[4]])
          } else if(is.na(params[[3]]) == FALSE && is.na(params[[4]])
                    && is.na(params[[6]]) == FALSE) {
            # Restrict by only genes.
            omim_info <- subset(omim_info, omim_info[["hgnc_symbol"]] == params[[6]])
          } else {
            # Restrict by only genes and chromosomes.
            omim_info <- subset(omim_info, omim_info[["chromosome_name"]] == params[[4]])
            omim_info <- subset(omim_info, omim_info[["hgnc_symbol"]] == params[[6]])
          }
        }
      } else {
        stop(paste("No information found in the OMIM database for the traits given."))
      }
    } else {
      stop(print("No OMIM ID found."))
    }
  } else if(is.na(params[[14]]) == FALSE && is.na(params[[3]])) {
    omim_info <- get_omim_genes(params[[14]])
    #Start the restiction from the other processes.
    if(nrow(omim_info) != 0) {
      if(is.na(params[[8]])){
        if(is.na(params[[4]]) && is.na(params[[6]])) {
          return(omim_info)
        } else if(is.na(params[[4]]) == FALSE && is.na(params[[6]])) {
          # Restrict by only chromosomes.
          omim_info <- subset(omim_info, omim_info[["chromosome_name"]] == params[[4]])
        } else if(is.na(params[[4]]) && is.na(params[[6]]) == FALSE) {
          # Restrict by only genes.
          omim_info <- subset(omim_info, omim_info[["hgnc_symbol"]] == params[[6]])
        } else {
          # Restrict by only genes and chromosomes.
          omim_info <- subset(omim_info, omim_info[["chromosome_name"]] == params[[4]])
          omim_info <- subset(omim_info, omim_info[["hgnc_symbol"]] == params[[6]])
        }
      }
    } else {
      stop(paste("No information found in the OMIM database for the traits given."))
    }
  }
}

#' @title Gets the information needed for the KEGG Api.
#' @description Depending on if the user wants to show values falling below
#' the value threshold or if they want the regular mode showing only those
#' found given in the restrictions, then the information is gathered for the
#' KEGG api.
#' @param index: An empty vectory.
#' @param matched: An empty vectory to be returned.
#' @param fg: An empty vectory to be returned.
#' @param bg: An empty vectory to be returned.
#' @param kegg_paths: The kegg pathways.
#' @param kegg_genes: The kegg genes for each pathway ID.
#' @param color: The color if the p-value threshold mode is set to 1.
#' @param params: All of the given user information in a list.
#' @param kegg_paths_low: The dataset contianing all information restricted by
#' p-value that is significant.
#' @param kegg_paths_high: The dataset contianing all information restricted by
#' p-value that is not able to be shown to be significant.
#' @param pval_low_col: Color representing variants by p-value that is not able
#' to be shown to be significant.
#' @param pval_high_col: Color representing variants by p-value are significant.
#' @return A list of containing two values, one with
#' the genes from a given pathway and the other with the
#' kegg path ids, and the kegg id.
.get_final_output_gwas <- function(index, matched, kegg_paths, kegg_genes,
                                   color, params, kegg_paths_low, kegg_paths_high,
                                   pval_low_col, pval_high_col, fg, bg, ko_ids,
                                   ko_matched, kegg_paths_mixed, pval_mix_col){
  if (anyNA(unique(kegg_paths[["mapped_gene"]])) == FALSE) {
    if (params[[7]] == 1 && params[[13]] == 1) {
      # No P-value threshold save for colors and gene mode is 1.
      matched_info <- .mapped_info_helper(kegg_paths, kegg_genes, ko_ids, color)
      matched <- matched_info[[1]]
      fg <- matched_info[[2]]
      bg <- matched_info[[3]]
      ko_matched <- matched_info[[4]]
    } else if(params[[7]] == 1 && params[[13]] == 0){
      # P-value threshold save for colors and gene mode is 1.
      matched_info <- .mapped_low_paths_info_helper(kegg_paths_low, kegg_genes, ko_ids, pval_low_col)
      matched <- matched_info[[1]]
      fg <- matched_info[[2]]
      bg <- matched_info[[3]]
      ko_matched <- matched_info[[4]]
      
      matched_info <- .mapped_high_paths_info_helper(kegg_paths_high, kegg_genes, ko_ids, pval_high_col)
      matched <- c(matched, matched_info[[1]])
      fg <- c(fg, matched_info[[2]])
      bg <- c(bg, matched_info[[3]])
      ko_matched <- c(ko_matched, matched_info[[4]])
      
      if(length(kegg_paths_mixed) != 0) {
        matched_info <- .mapped_mixed_paths_info_helper(kegg_paths_mixed, kegg_genes, ko_ids, pval_mix_col)
        matched <- c(matched, matched_info[[1]])
        fg <- c(fg, matched_info[[2]])
        bg <- c(bg, matched_info[[3]])
        ko_matched <- c(ko_matched, matched_info[[4]])
      }
    } else if (params[[7]] == 0 && params[[13]] == 1) {
      # No P-value threshold save for colors and gene mode is 0.
      matched_info <- .reported_info_helper(kegg_paths, kegg_genes, ko_ids, color)
      matched <- matched_info[[1]]
      fg <- matched_info[[2]]
      bg <- matched_info[[3]]
      ko_matched <- matched_info[[4]]
    } else if (params[[7]] == 0 && params[[13]] == 0) {
      # P-value threshold save for colors and gene mode is 0.
      matched_info <- .reported_low_paths_info_helper(kegg_paths_low, kegg_genes, ko_ids, pval_low_col)
      matched <- matched_info[[1]]
      fg <- matched_info[[2]]
      bg <- matched_info[[3]]
      ko_matched <- matched_info[[4]]
      
      matched_info <- .reported_high_paths_info_helper(kegg_paths_high, kegg_genes, ko_ids, pval_high_col)
      matched <- c(matched, matched_info[[1]])
      fg <- c(fg, matched_info[[2]])
      bg <- c(bg, matched_info[[3]])
      ko_matched <- c(ko_matched, matched_info[[4]])
      
      if(length(kegg_paths_mixed) != 0) {
        matched_info <- .reported_mixed_paths_info_helper(kegg_paths_mixed, kegg_genes, ko_ids, pval_mix_col)
        matched <- c(matched, matched_info[[1]])
        fg <- c(fg, matched_info[[2]])
        bg <- c(bg, matched_info[[3]])
        ko_matched <- c(ko_matched, matched_info[[4]])
      }
    } else if (params[[7]] == 2 && params[[13]] == 1) {
      # No P-value threshold save for colors and gene mode is 2.
      matched_info <- .mapped_info_helper(kegg_paths, kegg_genes, ko_ids, color)
      matched <- matched_info[[1]]
      fg <- matched_info[[2]]
      bg <- matched_info[[3]]
      ko_matched <- matched_info[[4]]
      
      matched_info <- .reported_info_helper(kegg_paths, kegg_genes, ko_ids, color)
      matched <- c(matched, matched_info[[1]])
      fg <- c(fg, matched_info[[2]])
      bg <- c(bg, matched_info[[3]])
      ko_matched <- c(ko_matched, matched_info[[4]])
      
      if(length(kegg_paths_mixed) != 0) {
        matched_info <- .reported_mixed_paths_info_helper(kegg_paths_mixed, kegg_genes, ko_ids, pval_mix_col)
        matched <- c(matched, matched_info[[1]])
        fg <- c(fg, matched_info[[2]])
        bg <- c(bg, matched_info[[3]])
        ko_matched <- c(ko_matched, matched_info[[4]])
      }
    } else if (params[[7]] == 2 && params[[13]] == 0) {
      # P-value threshold save for colors and gene mode is 2.
      matched_info <- .mapped_low_paths_info_helper(kegg_paths_low, kegg_genes, ko_ids, pval_low_col)
      matched <- matched_info[[1]]
      fg <- matched_info[[2]]
      bg <- matched_info[[3]]
      ko_matched <- matched_info[[4]]
      
      matched_info <- .mapped_high_paths_info_helper(kegg_paths_high, kegg_genes, ko_ids, pval_high_col)
      matched <- c(matched, matched_info[[1]])
      fg <- c(fg, matched_info[[2]])
      bg <- c(bg, matched_info[[3]])
      ko_matched <- c(ko_matched, matched_info[[4]])
      
      if(length(kegg_paths_mixed) != 0) {
        matched_info <- .mapped_mixed_paths_info_helper(kegg_paths_mixed, kegg_genes, ko_ids, pval_mix_col)
        matched <- c(matched, matched_info[[1]])
        fg <- c(fg, matched_info[[2]])
        bg <- c(bg, matched_info[[3]])
        ko_matched <- c(ko_matched, matched_info[[4]])
      }
      
      matched_info <- .reported_low_paths_info_helper(kegg_paths_low, kegg_genes, ko_ids, pval_low_col)
      matched <- c(matched, matched_info[[1]])
      fg <- c(fg, matched_info[[2]])
      bg <- c(bg, matched_info[[3]])
      ko_matched <- c(ko_matched, matched_info[[4]])
      
      matched_info <- .reported_high_paths_info_helper(kegg_paths_high, kegg_genes, ko_ids, pval_high_col)
      matched <- c(matched, matched_info[[1]])
      fg <- c(fg, matched_info[[2]])
      bg <- c(bg, matched_info[[3]])
      ko_matched <- c(ko_matched, matched_info[[4]])
      
      if(length(kegg_paths_mixed) != 0) {
        matched_info <- .reported_mixed_paths_info_helper(kegg_paths_mixed, kegg_genes, ko_ids, pval_mix_col)
        matched <- c(matched, matched_info[[1]])
        fg <- c(fg, matched_info[[2]])
        bg <- c(bg, matched_info[[3]])
        ko_matched <- c(ko_matched, matched_info[[4]])
      }
    }
  }
  
  return(list(matched, fg, bg,  ko_matched))
}

#' @title Gets the information needed for the KEGG Api.
#' @description Depending on if the user wants to show values falling below
#' the value threshold or if they want the regular mode showing only those
#' found given in the restrictions, then the information is gathered for the
#' KEGG api.
#' @param index: An empty vectory.
#' @param matched: An empty vectory to be returned.
#' @param fg: An empty vectory to be returned.
#' @param bg: An empty vectory to be returned.
#' @param kegg_paths: The kegg pathways.
#' @param kegg_genes: The kegg genes for each pathway ID.
#' @param color: The color if the p-value threshold mode is set to 1.
#' @param params: All of the given user information in a list.
#' @param kegg_paths_low: The dataset contianing all information restricted by
#' p-value that is significant.
#' @param kegg_paths_high: The dataset contianing all information restricted by
#' p-value that is not able to be shown to be significant.
#' @param pval_low_col: Color representing variants by p-value that is not able
#' to be shown to be significant.
#' @param pval_high_col: Color representing variants by p-value are significant.
#' @return A list of containing two values, one with
#' the genes from a given pathway and the other with the
#' kegg path ids, and the kegg id.
.get_final_output_omim <- function(index, matched, kegg_paths, kegg_genes,
                                   color, params, fg, bg, ko_ids, ko_matched){
  #Get the colors for the matched genes.
  for (i in 1:length(unique(kegg_paths[["hgnc_symbol"]]))) {
    if (anyNA(match(unique(kegg_paths[["hgnc_symbol"]])[i], unlist(kegg_genes))) == FALSE) {
      index <- c(index, match(unique(kegg_paths[["hgnc_symbol"]])[i], unlist(kegg_genes)))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg <- replicate(length(matched), color)
  bg <- replicate(length(matched), "#000000")
  
  return(list(matched, fg, bg, ko_matched))
}

#---- KEGG Helper Functions ----

#' @title Helper for VarGen pipeline dataset.
#' @description Finds the matched values for each of the
#' values in the matched genes column in the VarGen pipeline
#' dataset.
#' @param data: Data set containing resticted values.
#' @param kegg_genes: Genes to be matched by.
#' @param color: Color for foundvariants in pathway.
#' @return The matched values
.info_helper <- function(data, kegg_genes, ko_ids, color){
  matched <- c()
  index <- c()
  ko_matched <- c()
  
  for (i in 1:length(unique(data[["hgnc_symbol"]]))) {
    if (anyNA(match(unique(data[["hgnc_symbol"]])[i], unlist(kegg_genes))) == FALSE) {
      index <- c(index, match(unique(data[["hgnc_symbol"]])[i], unlist(kegg_genes)))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg <- replicate(length(matched), color)
  bg <- replicate(length(matched), "#000000")
  
  return(list(matched, fg, bg, ko_matched))
}

#' @title Helper for mapped genes in p-value restricting
#' @description Helper for mapped genes in p-value restricting that finds the
#' matched values for each of the values in the matched genes column in the
#' kegg_paths.
#' @param kegg_paths: Data set containing resticted values.
#' @param kegg_genes: Genes to be matched by.
#' @param color: Color for foundvariants in pathway.
#' @return The matched values
.mapped_info_helper <- function(kegg_paths, kegg_genes, ko_ids, color){
  matched <- c()
  index <- c()
  ko_matched <- c()
  for (i in 1:length(unique(kegg_paths[["mapped_genes"]]))) {
    if (anyNA(match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes))) == FALSE) {
      index <- c(index, match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes)))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg <- replicate(length(matched), color)
  bg <- replicate(length(matched), "#000000")
  
  return(list(matched, fg, bg, ko_matched))
}


#' @title Helper for reported genes in p-value restricting
#' @description Helper for mapped genes in p-value restricting that finds the
#' matched values for each of the values in the matched genes column in the
#' kegg_paths.
#' @param kegg_paths: Data set containing resticted values.
#' @param kegg_genes: Denes to be matched by.
#' @param color: Color for found variants in pathway.
#' @return The matched values
.reported_info_helper <- function(kegg_paths, kegg_genes, ko_ids, color){
  matched <- c()
  index <- c()
  ko_matched <- c()
  for (i in 1:length(unique(kegg_paths[["reported_genes"]]))) {
    if (is.na(match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes))) == FALSE) {
      # This might cause an error... needs testing.
      append(index, c(index, match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes))), length(index))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg <- replicate(length(matched), color)
  bg <- replicate(length(matched), "#000000")
  
  return(list(matched, fg, bg, ko_matched))
}

#' @title Helper for reported genes in p-value restricting for low threshold
#' coloring
#' @description Helper for reported genes in p-value restricting that finds the
#' matched values for each of the values in the matched genes column in the
#' kegg_paths.
#' @param kegg_paths: data set containing resticted values.
#' @param kegg_genes: genes to be matched by.
#' @param pval_low_col: The color of p-values that are significant given the
#' p-value threshold aka the alpha.
#' @return The matched values
.reported_low_paths_info_helper <- function(kegg_paths, kegg_genes, ko_ids, pval_low_col){
  matched <- c()
  index <- c()
  ko_matched <- c()
  for (i in 1:length(unique(kegg_paths[["reported_genes"]]))) {
    if (anyNA(match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes))) == FALSE) {
      index <- c(index, match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes)))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg_low <- replicate(length(matched), pval_low_col)
  bg_low <- replicate(length(matched), "#000000")
  
  return(list(matched, fg_low, bg_low, ko_matched))
}

#' @title Helper for reported genes in p-value restricting for mixed threshold
#' coloring
#' @description Helper for reported genes in p-value restricting that finds the
#' matched values for each of the values in the matched genes column in the
#' kegg_paths.
#' @param kegg_paths: data set containing resticted values.
#' @param kegg_genes: genes to be matched by.
#' @param pval_low_col: The color of p-values that are significant given the
#' p-value threshold aka the alpha.
#' @return The matched values
.reported_mixed_paths_info_helper <- function(kegg_paths, kegg_genes, ko_ids, pval_mixed_col){
  matched <- c()
  index <- c()
  ko_matched <- c()
  for (i in 1:length(unique(kegg_paths[["reported_genes"]]))) {
    if (anyNA(match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes))) == FALSE) {
      index <- c(index, match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes)))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg_low <- replicate(length(matched), pval_mixed_col)
  bg_low <- replicate(length(matched), "#000000")
  
  return(list(matched, fg_low, bg_low, ko_matched))
}

#' @title Helper for mapped genes in p-value restricting for low threshold
#' coloring
#' @description Helper for reported genes in p-value restricting that finds the
#' matched values for each of the values in the matched genes column in the
#' kegg_paths.
#' @param kegg_paths: data set containing resticted values.
#' @param kegg_genes: genes to be matched by.
#' @param pval_low_col: The color of p-values that are significant given the
#' p-value threshold aka the alpha.
#' @return The matched values
.mapped_low_paths_info_helper <- function(kegg_paths, kegg_genes, ko_ids, pval_low_col){
  matched <- c()
  index <- c()
  ko_matched <- c()
  for (i in 1:length(unique(kegg_paths[["mapped_genes"]]))) {
    if (anyNA(match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes))) == FALSE) {
      index <- c(index, match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes)))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg_low <- replicate(length(matched), pval_low_col)
  bg_low <- replicate(length(matched), "#000000")
  
  return(list(matched, fg_low, bg_low, ko_matched))
}

#' @title Helper for reported genes in p-value restricting for high threshold
#' coloring
#' @description Helper for reported genes in p-value restricting that finds the
#' matched values for each of the values in the matched genes column in the
#' kegg_paths.
#' @param kegg_paths: data set containing resticted values.
#' @param kegg_genes: genes to be matched by.
#' @param pval_nigh_col: The color of p-values that are not able to be dicidely shown
#' to be significant given the p-value threshold aka the alpha.
#' @return The matched values
.mapped_high_paths_info_helper <- function(kegg_paths, kegg_genes, ko_ids, pval_high_col){
  matched <- c()
  index <- c()
  ko_matched <- c()
  for (i in 1:length(unique(kegg_paths[["mapped_genes"]]))) {
    if (anyNA(match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes))) == FALSE) {
      index <- c(index, match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes)))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg_high <- replicate(length(matched), pval_high_col)
  bg_high <- replicate(length(matched), "#000000")
  
  return(list(matched, fg_high, bg_high, ko_matched))
}

#' @title Helper for mapped genes in p-value restricting for mixed threshold
#' coloring
#' @description Helper for reported genes in p-value restricting that finds the
#' matched values for each of the values in the matched genes column in the
#' kegg_paths.
#' @param kegg_paths: data set containing resticted values.
#' @param kegg_genes: genes to be matched by.
#' @param pval_low_col: The color of p-values that are significant given the
#' p-value threshold aka the alpha.
#' @return The matched values
.mapped_mixed_paths_info_helper <- function(kegg_paths, kegg_genes, ko_ids, pval_mixed_col){
  matched <- c()
  index <- c()
  ko_matched <- c()
  for (i in 1:length(unique(kegg_paths[["mapped_genes"]]))) {
    if (anyNA(match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes))) == FALSE) {
      index <- c(index, match(unique(kegg_paths[["mapped_genes"]])[i], unlist(kegg_genes)))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg_low <- replicate(length(matched), pval_mixed_col)
  bg_low <- replicate(length(matched), "#000000")
  
  return(list(matched, fg_low, bg_low, ko_matched))
}

#' @title Helper for reported genes in p-value restricting for low threshold
#' coloring
#' @description Helper for reported genes in p-value restricting that finds the
#' matched values for each of the values in the matched genes column in the
#' kegg_paths.
#' @param kegg_paths: data set containing resticted values.
#' @param kegg_genes: genes to be matched by.
#' @param pval_nigh_col: The color of p-values that are not able to be dicidely shown
#' to be significant given the p-value threshold aka the alpha.
#' @return The matched values
.reported_high_paths_info_helper <- function(kegg_paths, kegg_genes, ko_ids, pval_high_col){
  matched <- c()
  index <- c()
  ko_matched <- c()
  for (i in 1:length(unique(kegg_paths[["reported_genes"]]))) {
    if (anyNA(match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes))) == FALSE) {
      index <- c(index, match(unique(kegg_paths[["reported_genes"]])[i], unlist(kegg_genes)))
      for (j in 1:length(index)) {
        matched <- c(matched, unlist(kegg_genes)[index[j]])
        ko_matched <- c(ko_matched, unlist(ko_ids)[index[j]])
      }
    }
  }
  
  fg_high <- replicate(length(matched), pval_high_col)
  bg_high <- replicate(length(matched), "#000000")
  
  return(list(matched, fg_high, bg_high, ko_matched))
}

#---- LD Functions ----

#' @title Population code finders for 1000 genomes.
#' @description Finds the population code for LD plots and dataframe.
#' @returns Dataframe with codes and description of each
#' # population availble to generate a dataframe with.
get_LD_populations <- function() {
  tryCatch(
    {
      data <- httr::GET(paste("https://rest.ensembl.org",
                              "/info/variation/populations/homo_sapiens?filter=LD",
                              sep = ""),
                        httr::content_type("application/json"))
      data <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(data)))
      data <- data.frame(data, do.call(rbind, stringr::str_split(data[["name"]], ":")))
      
      data <- data[, c(2,6)]
      colnames(data) <- c("Description", "Code")
      
      return(data)
    },
    error = function(UnableToRetrieveDataError) {
      print(paste0("ERROR in get_LD_populations: Either the data was queried properly ",
                   "form the site or couldn't be formatted from url ",
                   "'https://rest.ensembl.org/info/variation/populations/homo_sapiens?filter=LD'.",
                   sep = ""))
      
      print(paste0("R thown error:"))
      message(UnableToRetrieveDataError)
    }
  )
}

#' @title Make LD score dataframe.
#' @description Makes the LD dataframe using the Ensembl API for
#' the 1000 genomes project and a given population as
#' selected by the user.  The population has a code associated
#' with it.
#' @param rsid: rsid to compare to population for LD.
#' @param pop_code: Population code given and outlines by Ensembl and the
#' 1000 genomes project conventions.  This can be found with the
#' get_LD_populations function.
#' @param d_prime_thresh: The D' alpha value that the dataset can be
#' resticted by.
#' @returns Dataframe with LD scores of a specific population.
get_LD_scores <- function(rsid, pop_code, d_prime_thresh = NULL) {
  # 1. Retrieve the data from the 1000 genomes project into a usable form of conversion
  # into a Grange.
  tryCatch(
    {
      get_output <- httr::GET(paste0("https://rest.ensembl.org/ld/human/",
                                     rsid, "/1000GENOMES:phase_3:", pop_code,
                                     "?", sep = ""),
                              httr::content_type("application/json"))
      get_output <- data.frame(jsonlite::fromJSON(jsonlite::toJSON(httr::content(get_output))))
      
      get_output$r2 <- as.numeric(as.character(get_output$r2))
      get_output$r2 <- as.numeric(as.character(get_output$d_prime))
      
      # 2. Restrict the dataset by the D' value.
      if(is.null(d_prime_thresh) == FALSE) {
        get_output <- subset(get_output, get_output$d_prime >= d_prime_thresh)
      }
      
      return(get_output)
    },
    error = function(UnableToRetrieveDataError) {
      print(paste0("ERROR in get_LD_scores: Either the data was queried properly ",
                   "form the site or couldn't be formatted from url ",
                   "'https://rest.ensembl.org/ld/human/rsid/1000GENOMES:phase_3:pop_code'.",
                   sep = ""))
      print(paste0("R thown error:"))
      message(UnableToRetrieveDataError)
    }
  )
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
#' @param traits: A vector containing the gwas traits that can be seached for.
#' Can be given as either a vector or as a single string.
#' @param chrs: A vector containing the chromosomes that can be seached for.
#' Can be given as either a vector or as a single string. Ex: "chr1"
#' @param title: The title that each of the graphics will be given with a
#' counter added.
#' @param pval_thresh: The cut off threshold for the p-values of the variants.
#' @param omim_ids: The omim id or ids to be searched only in a KEGG pathway.
#' @param rsids: A single snp rsid identifier or a vector of identifiers which
#' will return a list of found information from the gwas catalog in a dataframe.
#' @param genes: The mapped and reported genes that are found in the gwas file.
#' @param gene_mode: The mode for selecting if the reported genes and the
#' mapped genes should both be search for the given genes or if only one or the
#' other should be.  0 is for mapped genes only, 1 is for only reported, 2 is
#' for both. The defualt is 2.
#' @param color: The color for the variants given no p-value threshold has been given.
#' @param pval_low_col: The color for the genes given both a low p-value is found.
#' @param pval_high_col: The color for the genes given both a high p-value is found.
#' @param pval_mix_col: The color for the genes given both a low and high p-value is found.
#' @param pval_thresh_show: Show the p-value or not. 1 if not showing values and no
#' separation of colors if a p-value has been given and 0 if showing separation of p-value
#' by colors
#' @return Nothing; The KEGG pathways figures in a given or default directory.
kegg_search_graph <- function(vargen_dir = NA, output_dir = NA, traits = NA, chrs = NA, title = NA,
                              genes = NA, gene_mode = 2, rsids = NA, omim_ids = NA,
                              pval_thresh = NA, color = "#da8cde",
                              pval_low_col = "#14e06a", pval_high_col = "#14c5e0",
                              pval_mix_col = "#14e0c1", pval_thresh_show = 1) {
  # 1. Check if the parameters are valid.
  params = .kegg_validity_check(list(vargen_dir, output_dir, traits, chrs, title,
                                     genes, gene_mode, rsids, pval_thresh, color,
                                     pval_low_col, pval_high_col, pval_thresh_show,
                                     omim_ids))
  omim_info <- c()
  variants_traits <- c()
  # 2. Restrict the gwas catelog by the given paramaters.
  if(is.na(params[[14]])) {
    variants_traits <- .restrict_snp_gwas(params)
    if(length(variants_traits) == 0) {
      stop(paste("No information found in the GWAS Catalog for the restictions given."))
    } else {
      print(paste("Number of variants found with current restictions: ", length(variants_traits)))
    }
  } else {
    # 3. Get the OMIM variants.
    variants_traits <- .restrict_snp_omim(omim_info, params)
    if(length(variants_traits) == 0) {
      stop(paste("No information found in the GWAS Catalog for the restictions given."))
    } else {
      print(paste("Number of variants found with current restictions: ", length(variants_traits)))
    }
  }
  
  # 4. Get KEGG pathways for the filtered variants above.
  if(length(variants_traits) != 0  &&
     (length(!(variants_traits$"SNP_GENE_IDS" %in% "")) != 0 ||
      length(!(variants_traits$"ensembl_gene_id" %in% "")) != 0)) {
    if(is.na(params[[14]])) {
      if(is.na(omim_ids) == FALSE) {
        kegg_paths <- biomaRt::getBM(
          attributes = c("kegg_enzyme", "ensembl_gene_id"),
          filters = c("ensembl_gene_id"),
          values = variants_traits$"ensembl_gene_id",
          mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
      } else {
        kegg_paths <- biomaRt::getBM(
          attributes = c("kegg_enzyme", "ensembl_gene_id"),
          filters = c("ensembl_gene_id"),
          values = variants_traits$"SNP_GENE_IDS",
          mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
      }
    } else {
      if(is.na(omim_ids) == FALSE) {
        kegg_paths <- biomaRt::getBM(
          attributes = c("kegg_enzyme", "ensembl_gene_id"),
          filters = c("ensembl_gene_id"),
          values = variants_traits$"ensembl_gene_id",
          mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
      } else {
        kegg_paths <- biomaRt::getBM(
          attributes = c("kegg_enzyme", "ensembl_gene_id"),
          filters = c("ensembl_gene_id"),
          values = variants_traits$"SNP_GENE_IDS",
          mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
      }
    }
  } else {
    stop(print("No variants found for given restrictions."))
  }
  
  # 5. Make the directory and title as fit given the parameter inputs.
  if(nrow(kegg_paths) == 0 || is.na(unique(kegg_paths[["kegg_enzyme"]])[1])){
    stop(paste0("No KEGG pathways found.  Try broading the search."))
  } else {
    # Make output directory.
    output_dir <- .make_dir(params[[2]], "KEGG_images")
    # Check if the title is given.
    titleKey <- FALSE
    if(is.na(params[[5]])) {
      titleKey <- TRUE
    }
    
    # Get the information out of the grange.
    if(is.na(params[[14]])) {
      variant_info <- cbind(variants_traits$"SNPS", variants_traits$"SNP_GENE_IDS",
                            variants_traits$"P-VALUE", variants_traits$"MAPPED_GENE",
                            variants_traits$"REPORTED GENE(S)")
      colnames(variant_info) <- cbind("rsid", "ensembl_gene_id", "p_value",
                                      "mapped_genes", "reported_genes")
      variant_info <- merge(x = variant_info, y = kegg_paths,
                            by = "ensembl_gene_id", all.x = TRUE)
    } else {
      variant_info <- merge(x = variants_traits, y = kegg_paths,
                            by = "ensembl_gene_id", all.x = TRUE)
    }
    
    # 6. Subset and restrict the variants by the pathways and the p-values.
    # In this section the genes are also gone through for each mode (either
    # reported, mapped or both) as well.
    cnt <- 1
    filecnt <- 1
    # 6.1 Remove values without a KEGG pathway.
    kegg_paths <- subset(variant_info, is.na(variant_info[["kegg_enzyme"]]) == FALSE)
    kegg_paths <- subset(kegg_paths, kegg_paths[["kegg_enzyme"]] != "")
    print(paste("Number of KEGG pathways found for variants via the kegg api: ", nrow(kegg_paths)))
    
    # 6.2 Get the values that are missing KEGG pathways.
    kegg_paths_na <- subset(variant_info, is.na(variant_info[["kegg_enzyme"]]) == TRUE)
    kegg_paths_na <- subset(kegg_paths_na, kegg_paths_na[["kegg_enzyme"]] != "")
    
    if(is.na(params[[9]]) == FALSE){
      kegg_paths_na <- subset(kegg_paths_na, kegg_paths_na[["p_value"]] >= is.na(params[[9]]))
    }
    # 6.3 Write output to text file in same directory of values not found with KEGG pathways
    # in the next for loop.
    file_out <- file(paste(output_dir, "variants_without_kegg.txt"))
    
    kegg_paths_low <- kegg_paths
    kegg_paths_high <- kegg_paths
    kegg_paths_mixed <- c()
    # 6.4 Restrict it by p-value threshold
    if(is.na(pval_thresh) == FALSE){
      # Column must be changed to be numeric else it's not calculated properly.
      kegg_paths[["p_value"]] <- as.numeric(as.character(kegg_paths[["p_value"]]))
      if(params[[13]] == 0){
        # If the user wants to show the p-value threshold colors that are below
        # and above the threshold. For any p-values that are found to be both above and
        # below the threshold they will be colored according to the mixed color.
        # 6.4.1 Find the values in the high and low categories.
        kegg_paths_high <- kegg_paths[which(kegg_paths[["p_value"]] >= params[[9]]), ]
        kegg_paths_low <- kegg_paths[which(kegg_paths[["p_value"]] < params[[9]]), ]
        
        # 6.4.2 Find if there were any intersection between them and if there were remove and
        # add them to the mixed paths.
        intersect <- intersect(kegg_paths_high[["rsid"]], kegg_paths_low[["rsid"]])
        if(length(intersect) != 0) {
          kegg_paths_mixed.1 <- subset(kegg_paths_high, kegg_paths_high[["rsid"]] %in% intersect[["rsid"]])
          kegg_paths_mixed.2 <- subset(kegg_paths_low, kegg_paths_low[["rsid"]] %in% intersect[["rsid"]])
          kegg_paths_mixed <- merge(kegg_paths_mixed.1, kegg_paths_mixed.2, by = "rsid")
          
          kegg_paths_high[!kegg_paths_high[["rsid"]] %in% intersect[["rsid"]],]
          kegg_paths_low[!kegg_paths_low[["rsid"]] %in% intersect[["rsid"]],]
        }
        
      } else {
        kegg_paths <- kegg_paths[which(kegg_paths[["p_value"]] < params[[9]]), ]
      }
    }
    
    if(nrow(kegg_paths) == 0){
      writeLines(c(variant_info[["rsid"]][cnt]), file_out)
      print("No KEGG pathways available for this search request. Try broading search if possible.")
      print(paste0("File written in"), getwd())
    }
    
    # 7. Output KEGG pathways for each pathway in the output directory
    # assuming it doesn't already exist. Also write the file containing
    # the rsids that were aquired by filtered out due to the p-value
    # threshold.
    index <- c()
    for (kegg_pathway in kegg_paths[["kegg_enzyme"]]) {
      # 7.1 Get the genes from a given pathway.
      kegg_info <- .make_genes(kegg_pathway)
      kegg_genes <- kegg_info[1]
      kegg_path_id <- kegg_info[2]
      ko_ids <- kegg_info[4]
      
      fg <- c()
      bg <- c()
      matched <- c()
      ko_matched <- c()
      # 7.2 Find the unique values that match in the mapped or the reported genes
      # depending on the mode.
      if(is.na(params[[14]])) {
        final_info <- .get_final_output_gwas(index, matched, kegg_paths, kegg_genes,
                                             color, params, kegg_paths_low, kegg_paths_high,
                                             pval_low_col, pval_high_col, fg, bg, ko_ids,
                                             ko_matched, kegg_paths_mixed, pval_mix_col)
        
        matched <- final_info[[1]]
        fg <- final_info[[2]]
        bg <- final_info[[3]]
        ko_matched <- final_info[[4]]
        
      } else {
        final_info <- .get_final_output_omim(index, matched, kegg_paths, kegg_genes,
                                             color, params, fg, bg, ko_ids,
                                             ko_matched)
        matched <- final_info[[1]]
        fg <- final_info[[2]]
        bg <- final_info[[3]]
        ko_matched <- final_info[[4]]
      }
      
      # 7.3 Make sure that some value is found so that a pathway can be made.
      if (is.null(unlist(kegg_genes)) == FALSE
          && identical(unlist(kegg_genes), character(0)) == FALSE
          && length(unlist(kegg_genes, use.names=FALSE)) != 0
          && length(matched) != 0) {
        
        matched <- unlist(matched)
        fg <- unlist(fg)
        bg <- unlist(bg)
        
        df <- data.frame(unlist(kegg_genes), unlist(kegg_info[[3]]))
        df <- df[is.element(df$"unlist.kegg_genes.", matched),]
        colnames(df) <- cbind("matched", "id")
        
        # Get unique fg and bg colors for the macthed information as duplicated can occur if reported
        # and mapped have same information.
        temp_df <- data.frame(matched, fg, bg, ko_matched)
        df <- merge(x = df, y = temp_df[!duplicated(temp_df[["matched"]]), ], by = "matched", all.x = TRUE)
        
        temp <- df[["id"]]
        for(i in 1:length(df[["id"]])) {
          temp[i] <- paste("hsa:", df[["id"]][i], sep = "")
        }
        kegg_id <- temp
        
        kegg_genes <- df[["matched"]]
        fg <- df[["fg"]]
        bg <- df[["bg"]]
        ko_matched <- df[["ko_matched"]]
        
        print(paste0("Number of unique variants found for pathway ", kegg_pathway, " : ", length(kegg_genes)))
        
        if(length(ko_matched) != 0) {
          url <- KEGGREST::color.pathway.by.objects(paste("path:", kegg_path_id, sep = ""),
                                                    unlist(ko_matched),
                                                    fg.color.list = unlist(fg),
                                                    bg.color.list = unlist(bg))
          
          # 6.4 Make the correct titles for each of the pathways.
          if (titleKey == TRUE) {
            title <- paste0("hsa_", kegg_path_id, "_kegg.png")
          } else {
            if (file.exists(paste0(output_dir, filecnt, "_", title)) == FALSE) {
              title <- paste0(cnt, "_", title)
            }
          }
          
          if (file.exists(paste(output_dir, title, sep = ""))) {
            file_cnt <- sapply(output_dir, function(output_dir) {
              length(list.files(output_dir, pattern = title))
            })
            title <- paste0(file_cnt, "_", title, sep = "")
          }
          # 6.6 Download the files from the website.
          download.file(url, paste(output_dir, title, sep = ""), mode = "wb")
          filecnt <- filecnt + 1
        } else {
          print("Added variant to file containing those not found in directory.")
          writeLines(c(kegg_paths[["rsid"]][cnt]), file_out)
        }
      } else {
        # 6.5 Write the output files for the paths not found.
        if (is.na(params[[14]] == FALSE)) {
          print("Added variant to file containing those not found in directory.")
          writeLines(c(kegg_paths[["rsid"]][cnt]), file_out)
        } else {
          print(kegg_paths[["mim_morbid_accession"]][cnt])
          write.table(x = kegg_paths[["mim_morbid_accession"]], file = file_out)
        }
      }
      cnt <- cnt + 1
    }
  }
}

#' @title Makes the KEGG graphs from VarGen pipeline
#' @description Allows for the KEGG graphs to be output into a directory,
#' from a dataset generated from the VarGen pipeline. The function takes
#' data directly from the VarGen pipeline and vargen_custom. please refer to:
#' \code{\link{vargen_pipeline}}
#' @param data: The output from the VarGen pipeline.
#' download to or is being stored in.
#' @param color: The color for when a gene is found found in a pathway.
#' @param output_dir: The directory where the files will be put as desired
#' @return Nothing; The KEGG pathways figures in a given or default directory.
kegg_vargen_graph <- function(data, output_dir = NULL, color = "#da8cde") {
  # 1. Check the dataset and the output directory.  Make the directory if necessary.
  .check_dataframe(data)
  output_dir <- .check_output_dir(output_dir, "KEGG_images")
  
  kegg_paths <- biomaRt::getBM(
    attributes = c("kegg_enzyme", "ensembl_gene_id", "entrezgene_id"),
    filters = c("ensembl_gene_id"),
    values = data[["ensembl_gene_id"]],
    mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
  
  variants_info <- c()
  
  if(length(unique(kegg_paths[["kegg_enzyme"]])) != 1) {
    variants_info <- merge(data[, c(1,4:5)], kegg_paths, by = c("ensembl_gene_id"),
                           all.x = TRUE)
    # 3.1 Remove any newly added variants that were not found by vargen.
    variants_found <- subset(variants_info, variants_info[["ensembl_gene_id"]] != "")
    # 3.2 Remove those without an Entrez ID as they can't be used.
    variants_Entrez <- subset(variants_found, is.na(variants_found[["entrezgene_id"]]) == FALSE)
    print(paste0("Number of variants found with Entrez ID: ", nrow(variants_Entrez), " out of ",
                 nrow(variants_found)))
    # 3.3 Remove those without a Kegg pathway as they also can't be used.
    variants_info_kegg <- subset(variants_Entrez, is.na(variants_Entrez[["kegg_enzyme"]]) == FALSE)
    variants_info_kegg <- subset(variants_info_kegg, variants_info_kegg[["kegg_enzyme"]] != "")
    variants_info_kegg <- subset(variants_info_kegg,
                                 !duplicated(subset(variants_info_kegg,
                                                    select=c("kegg_enzyme", "hgnc_symbol"))))
    print(paste0("Number of variants found with KEGG ID: ", nrow(variants_info_kegg), " out of ",
                 nrow(variants_Entrez)))
    variants_info <- variants_info_kegg
  } else if(na(unique(kegg_paths[["kegg_enzyme"]])[1]) && length(unique(kegg_paths[["kegg_enzyme"]])) == 1) {
    stop(paste("No KEGG pathways are available for the variants found in the KEGG database"))
  }
  
  for (kegg_pathway in variants_info[["kegg_enzyme"]]) {
    # 4. Get the matching genes.
    kegg_info <- .make_genes(kegg_pathway)
    kegg_genes <- kegg_info[1]
    kegg_path_id <- kegg_info[2]
    ko_ids <- kegg_info[4]
    
    fg <- c()
    bg <- c()
    matched <- c()
    ko_matched <- c()
    # 5. Find the unique values that match in the mapped or the reported genes
    # depending on the mode.
    matched_info <- .info_helper(variants_info, kegg_genes, ko_ids, color)
    matched <- matched_info[[1]]
    fg <- matched_info[[2]]
    bg <- matched_info[[3]]
    ko_matched <- matched_info[[4]]
    
    if (is.null(unlist(kegg_genes)) == FALSE
        && identical(unlist(kegg_genes), character(0)) == FALSE
        && length(unlist(kegg_genes, use.names=FALSE)) != 0
        && length(matched) != 0) {
      
      matched <- unlist(matched)
      fg <- unlist(fg)
      bg <- unlist(bg)
      
      df <- data.frame(unlist(kegg_genes), unlist(kegg_info[[3]]))
      df <- df[is.element(df$"unlist.kegg_genes.", matched),]
      colnames(df) <- cbind("matched", "id")
      
      # Get unique fg and bg colors for the macthed information as duplicated can occur if reported
      # and mapped have same information.
      temp_df <- data.frame(matched, fg, bg, ko_matched)
      df <- merge(x = df, y = temp_df[!duplicated(temp_df[["matched"]]), ], by = "matched", all.x = TRUE)
      
      temp <- df[["id"]]
      for(i in 1:length(df[["id"]])) {
        temp[i] <- paste("hsa:", df[["id"]][i], sep = "")
      }
      kegg_id <- temp
      
      kegg_genes <- df[["matched"]]
      fg <- df[["fg"]]
      bg <- df[["bg"]]
      ko_matched <- df[["ko_matched"]]
      
      print(paste0("Number of unique variants found for pathway ", kegg_pathway, " : ", length(kegg_genes)))
      
      if(length(ko_matched) != 0) {
        url <- KEGGREST::color.pathway.by.objects(paste("path:", kegg_path_id, sep = ""),
                                                  unlist(ko_matched),
                                                  fg.color.list = unlist(fg),
                                                  bg.color.list = unlist(bg))
        print(paste(output_dir, kegg_pathway, ".png", sep = ""))
        download.file(url, paste(output_dir, kegg_pathway, ".png", sep = ""), mode = "wb")
      }
    }
  }
}

#' @title Visualize the online OMIM Graphs.
#' @description Takes the user to the OMIM website's graphic for a given OMIM ID. Users
#' can choose either radial or linear as a mode - with the default being linear. For
#' more information on OMIM's PheneGene graphics, please see the offical OMIM website.
#' @param type: the type of graph to go to either as 'linear' or 'radial'.
#' The type can also be given as 'l', 'L' or 'r', 'R'.
#' @param omim_id: A single OMIM ID given as a string.
#' @return Nothing; Takes user to url of omim id selected.
omim_graph_online <- function(type = "linear", omim_id){
  # 1. Check the type to see if it is valid then check if the omim id is not null or na.
  type <- tolower(type)
  if(!(type %in% c("linear", "l", "r", "radial"))) {
    stop(print(paste0("There is no option '", type,
                      "' for parameter type. Options are 'linear', 'radial', ",
                      "'l', 'L' or 'r', 'R'.")))
  }
  if(is.null(omim_id) || is.na(omim_id)) {
    stop(print("There is no OMIM ID given."))
  }
  
  # 2. Get the url from the OMIM database depending on the omim id.
  if(type %in% c("linear", "l")) {
    url <- paste0("https://omim.org/graph/linear/", omim_id, sep = "")
    browseURL(url)
  } else {
    url <- paste0("https://omim.org/graph/radial/", omim_id, sep = "")
    browseURL(url)
  }
}

#' @title Makes a treemap using data from the VarGen pipeline.
#' @description Makes a treemap from the VarGen pipeline using the treemap package.
#' The number of variants are counted per gene which is then displayed are a treemap.
#' Each block is colored accoding to the raw variant counts and the size is also dependent
#' on the raw variant counts. Normalized variant counts are available to be used if the
#' mode is set to be 'norm'. The variants with the highest counts can also be
#' retrieved given at an n amount using 'top_thresh' parameter. In Ontology Mode two graphs
#' will be displayed. One shows all information while the other allows for the user to
#' iteractivly click through each catagory/matched gene to see the traits. The graph with all
#' information displayed may be a bit hard to read, soit should be used as a guide for exploring
#' the data if large subsets are origionally graphed. This method only works with the annotated
#' VarGen pipeline as it requires the OMIM ID's in the final column of the dataset with the column
#' name of 'traits'.
#' @param data: The output from the VarGen pipeline.
#' @param title: The title of the figure.
#' @param top_thresh: The number of the top most frequent hits of counts in the dataset
#' relating to the variants found on that gene.
#' @param mode: Sets the counts to be either the raw or normalized counts.
#' #' (# variants / kbp). The modes are 'raw' for the raw counts, 'norm' for
#' normalized counts using # variants/kbp.
#' @param ontology: Sets the ontology mode which allows the traits to be grouped on their
#' corresponding gene.
#' @param ontol_genes: The genes that are to be included in the ontology mode's graph.
#' This may be useful if a certain gene is obscuring a large portion of the other genes.
#' @return Nothing; A treemap graph.
treemap_graph <- function(data, title = "VarGen Treemap Graph", top_thresh = NULL,
                          mode = "raw", ontology = FALSE, ontol_genes = NULL) {
  # 1. Check that each of the parameters to see if they are valid.
  # 1.1 Check if the dataset is null, na, or contains no data.
  .check_dataframe(data)
  
  # 1.2 Check if the threshold is a number or can be converted to one if
  # given as a string.
  if(typeof(top_thresh) != "double" && is.na(top_thresh) == FALSE
     && is.null(top_thresh) == FALSE) {
    if(is.na(as.double(top_thresh))) {
      print(paste0("The top threshold '", top_thresh,
                   "' could not be converted to a numeric value, ",
                   "so threshold will be used."))
      top_thresh <- NULL
    }
  }
  # 1.3 Check that the mode is either 'raw' or 'norm' and nothing else.
  if(!(mode %in% c("raw", "norm"))) {
    print(paste0("Setting mode to raw as mode '", mode, "' is not a valid mode."))
    mode = "raw"
  }
  # 1.4 If the title is set to be null or na, then make the title empty.
  if(is.null(title) || is.na(title)) {
    title = ""
  }
  # 1.5 Check if the ontology is either true or false
  if(ontology != TRUE && ontology != FALSE) {
    ontology == FALSE
  }
  # 1.6 Give a warning to the user if restriction genes are given and
  # the ontology mode is not set.
  if(ontology != TRUE && is.null(ontol_genes) == FALSE) {
    print("Genes can only be restricted when using ontology mode.")
  }
  
  
  # Option: ontology is being counted or just the number of variants are being counted.
  if(ontology == FALSE) {
    # 2.a Get the count data for the genes.
    hgnc_cnt <- data.frame(table(data[["hgnc_symbol"]]))
  } else {
    # 2.b Get the count data for the ontology
    hgnc_cnt <- data.frame(table(data[,c("trait", "hgnc_symbol")]))
    hgnc_cnt <- subset(hgnc_cnt, hgnc_cnt[["Freq"]] != 0)
    # Option subset by a certain gene only for traits
    if(is.na(ontol_genes) == FALSE && is.null(ontol_genes) == FALSE) {
      attempt <- hgnc_cnt[hgnc_cnt[["hgnc_symbol"]] %in% ontol_genes,]
      if(nrow(attempt) != 0) {
        hgnc_cnt <- attempt
      } else {
        print(paste0("No traits found for gene '", ontol_genes,
                     "', so all genes were used."))
      }
    }
    colnames(hgnc_cnt) <- c("Var1", "Var2", "Freq")
  }
  
  # Option: normilization of the variant counts.
  if(mode == "norm") {
    # 2.1.a Get the gene length data for normalization.  All data that doesn't have a start
    # and stop position can't be normalized and is removed from the dataset. The user is
    # alerted to this point being removed.
    length_info <- get_omim_genes(omim_ids = data[["trait"]], connect_to_gene_ensembl())
    data <- merge(data[3:7], length_info, by = c("ensembl_gene_id", "hgnc_symbol"),
                  all.x = TRUE)
    data <- data[!is.na(data[["end_position"]]), ]
    data <- unique(data)
    
    if(length(hgnc_cnt) == 2) {
      # 2.2.a. Normalize the count data.
      colnames(hgnc_cnt) <- c("hgnc_symbol", "Freq")
      norm_hgnc <- merge(hgnc_cnt, data[, c(2,7,8)], by = c("hgnc_symbol"),
                         all.x = TRUE)
      norm_hgnc <- unique(norm_hgnc)
      norm_hgnc$"len" <- (norm_hgnc[["end_position"]] - norm_hgnc[["start_position"]])
      norm_hgnc$"Freq" <- (norm_hgnc[["Freq"]] / norm_hgnc[["len"]])
      hgnc_cnt <- norm_hgnc
      # rounding off the counts.
      hgnc_cnt[["Freq"]] <- round(hgnc_cnt[["Freq"]], digit = 3)
      colnames(hgnc_cnt) <- c("Var1", "Freq")
    } else {
      # 2.2.b. Normalize the count data for ontology.
      colnames(hgnc_cnt) <- c("trait", "hgnc_symbol", "Freq")
      norm_hgnc <- merge(hgnc_cnt, data[, c(2,7,8)], by = c("hgnc_symbol"),
                         all.x = TRUE)
      norm_hgnc <- unique(norm_hgnc)
      norm_hgnc$"len" <- (norm_hgnc[["end_position"]] - norm_hgnc[["start_position"]])
      norm_hgnc$"Freq" <- (norm_hgnc[["Freq"]] / norm_hgnc[["len"]])
      hgnc_cnt <- norm_hgnc
      # rounding off the counts.
      hgnc_cnt[["Freq"]] <- round(hgnc_cnt[["Freq"]], digit = 3)
      colnames(hgnc_cnt) <- c("Var1", "Var2", "Freq")
      
      #Check if no positions are found and report it to the user.
      if(is.na(unique(hgnc_cnt[["Freq"]])[1]) && nrow(hgnc_cnt) == 1) {
        stop(print("None of the genes found can be normalized due to not having positional data."))
      }
    }
  }
  
  # Option: restrict by n number of the most hit counts as given by the user.
  if(is.null(top_thresh) == FALSE) {
    hgnc_cnt <- hgnc_cnt[order(-hgnc_cnt[["Freq"]]),]
    hgnc_cnt <- hgnc_cnt[1:top_thresh,]
  }
  
  # 2.3. Alert the user to the genes that were removed due to not having positions
  # if any have been removed during the 'norm' normalization process.
  if(mode == "norm") {
    removed <- intersect(hgnc_cnt[["Var1"]], data[["hgnc_symbol"]])
    removed <- data.frame(removed)
    if(nrow(removed) != nrow(hgnc_cnt) && mode == "norm" && nrow(removed) != 0) {
      removed <- setdiff(hgnc_cnt[["Var1"]], removed[["V1"]])
      removed <- data.frame(removed)
      print(paste0("There were some genes without positions and couldn't be normalized",
                   " thus were removed: ", removed[["V1"]]))
    }
    
    if(length(na.omit(removed[["V1"]])) == length(na.omit(unique(hgnc_cnt[["Var1"]])))) {
      stop("No genes found with proper elements for this analysis with given data.")
    }
  }
  
  # 3. Add in the labels that include the counts and genes.
  if((length(hgnc_cnt) == 2 && mode == "raw") ||
     (length(hgnc_cnt) == 5 && mode == "norm")) {
    hgnc_cnt$"labels" <- paste(hgnc_cnt[["Var1"]], hgnc_cnt[["Freq"]], sep = ": ")
    
    # 4. Grpah the information gathered for the treemap package.
    if(nrow(hgnc_cnt) != 0) {
      x11()
      treemap::treemap(dtf = hgnc_cnt, index = "labels", vSize = "Freq",
                       type = "value", vColor = "Freq", title = title,
                       border.lwds=c(5,2), border.col=c("black","white"))
    }
  } else {
    hgnc_cnt$"labels" <- paste(hgnc_cnt[["Var1"]], hgnc_cnt[["Freq"]], sep = ": ")
    
    # 4. Graph the information gathered for the treemap package.
    if(nrow(hgnc_cnt) != 0) {
      x11()
      treemap::treemap(dtf = hgnc_cnt, index = c("Var2", "labels"), vSize = "Freq",
                       type = "value", vColor = "Freq", title = title,
                       overlap.labels = 0.5, align.labels = list(
                         c("left", "top"),
                         c("right", "bottom")
                       ), fontsize.labels = c(14, 11),
                       border.col = c("black","white"), border.lwds = c(5,2),
                       bg.labels = "#abe1f5", fontfamily.labels = c("serif", "sans"))
      
      d3treeR::d3tree(treemap::treemap(dtf = hgnc_cnt, index = c("Var2", "labels"), vSize = "Freq",
                                       type = "value", vColor = "Freq", title = title,
                                       overlap.labels = 0.5, align.labels = list(
                                         c("left", "top"),
                                         c("right", "bottom")
                                       ), fontsize.labels = c(14, 11),
                                       border.col = c("black","white"), border.lwds = c(5,2),
                                       bg.labels = "#abe1f5", fontfamily.labels = c("serif", "sans")))
      
    }
  }
}

#' @title Make pathview figures for a VarGen variants.
#' @description Uses the pathview package to make figures for the variants found via the
#' VarGen pipeline.
#' @param data: The output from the vargen pipeline.
#' @param output_dir: The directory where the files will be put as desired
#' by the user.
#' @param title: The title that each of the graphics will be given with a
#' counter added.
#' @param key_pos: The position of the key for the pathview figure given as either
#' "topleft", "bottomleft", "topright", or "bottomright".
#' @param top_thresh: The number of the top most found values in the dataset.
#' @return Nothing; The KEGG pathways figures in a given or made directory.
pathview_vargen <- function(data, output_dir = NULL, title = NULL,
                            top_thresh = NULL, key_pos = "topright") {
  # 1. Check input parameters to see if they are valid.
  # 1.1 Check if the dataset is null, na, or contains no data.
  .check_dataframe(data)
  # 1.2 Check if the output directory exists and if it doesn't then make it.
  output_dir <- .check_output_dir(output_dir, "pathview")
  # 1.3 If the title is set to be na, then make the title null so it will be the generic title.
  if(is.null(title) || is.na(title)) {
    title = NULL
  }
  # 1.4 Check if the threshold is a number or can be converted to one if
  # given as a string.
  if(typeof(top_thresh) != "double" && is.na(top_thresh) == FALSE
     && is.null(top_thresh) == FALSE) {
    if(is.na(as.double(top_thresh))) {
      print(paste0("The top threshold '", top_thresh,
                   "' could not be converted to a numeric value, ",
                   "so threshold will be used."))
      top_thresh <- NULL
    }
  }
  # 1.5 Check that the given key is a proper position and if not sets a default.
  key_pos <- tolower(key_pos)
  if(!(key_pos %in% c("topleft", "bottomleft", "topright", "bottomright"))) {
    print(paste0("There is no option '", key_pos,
                 "' for parameter key_pos. Options are 'topleft', 'bottomleft',
                      'topright', or 'bottomright'. The default 'topright' will be used."))
    key_pos <- "topright"
  }
  
  # 2. Get the Entrez IDs and KEGG pathway IDs for each of the variants.
  kegg_paths <- biomaRt::getBM(
    attributes = c("kegg_enzyme", "ensembl_gene_id", "entrezgene_id"),
    filters = c("ensembl_gene_id"),
    values = data[["ensembl_gene_id"]],
    mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
  
  # 3. Merge the datasets together if there are KEGG pathways found.
  variants_info <- c()
  if(is.na(unique(kegg_paths[["kegg_enzyme"]]))) {
    variants_info <- merge(data[, c(1,4:5)], kegg_paths, by = c("ensembl_gene_id"),
                           all.x = TRUE)
    # 3.1 Remove any newly added variants that were not found by vargen.
    variants_found <- subset(variants_info, variants_info[["ensembl_gene_id"]] != "")
    # 3.2 Remove those without an Entrez ID as they can't be used.
    variants_Entrez <- subset(variants_found, is.na(variants_found[["entrezgene_id"]]) == FALSE)
    print(paste0("Number of variants found with Entrez ID: ", nrow(variants_Entrez), " out of ",
                 nrow(variants_found)))
    # 3.3 Remove those without a Kegg pathway as they also can't be used.
    variants_info_kegg <- subset(variants_Entrez, is.na(variants_Entrez[["kegg_enzyme"]]) == FALSE)
    variants_info_kegg <- subset(variants_info_kegg, variants_info_kegg[["kegg_enzyme"]] != "")
    variants_info_kegg <- subset(variants_info_kegg,
                                 !duplicated(subset(variants_info_kegg,
                                                    select=c("kegg_enzyme", "hgnc_symbol"))))
    print(paste0("Number of variants found with KEGG ID: ", nrow(variants_info_kegg), " out of ",
                 nrow(variants_Entrez)))
    variants_info <- variants_info_kegg
  } else {
    stop(paste("No KEGG pathways are available for the variants found in the KEGG database"))
  }
  hgnc_cnt <- data.frame(table(variants_info[["hgnc_symbol"]]))
  colnames(hgnc_cnt) <- c("hgnc_symbol", "cnts")
  
  # 4. Restrict by the top threshold if it is not null.
  if(is.null(top_thresh) == FALSE) {
    hgnc_cnt <- hgnc_cnt[order(-hgnc_cnt[["cnts"]]),]
    hgnc_cnt <- hgnc_cnt[1:top_thresh,]
  }
  
  # 5. Get the sum of the total genes found, not just those with a pathway.
  # This makes the counts into frequencies.
  sum <- sum( hgnc_cnt[["cnts"]])
  for(i in 1:nrow(hgnc_cnt)) {
    hgnc_cnt[["cnts"]][[i]]  <- hgnc_cnt[["cnts"]][[i]] / sum
  }
  # 6. Get the restricted information since pathways are the limiting factor.
  hgnc_cnt <- merge(hgnc_cnt, variants_info, by = "hgnc_symbol", all = TRUE)
  
  # 7. Make the proper input for each of the pathways for the pathview and then create the Pathview images.
  # The gene data is the frequency of the counts of variants found by VarGen with entrez genes IDs as the
  # requested ID.  The pathway id is the kedd_id.
  old_dir = getwd()
  setwd(output_dir)
  for(kegg_pathway in hgnc_cnt[["kegg_enzyme"]]) {
    kegg_id <- unlist(strsplit(kegg_pathway, "\\+"))[1]
    sel_info <- unique(hgnc_cnt[, c(6,2)])
    sel_genes <- as.matrix(sel_info[, 2])
    rownames(sel_genes) <- c(sel_info[, 1])
    
    pathview::pathview(gene.data = sel_genes, pathway.id = kegg_id,
                       gene.idtype = "entrez", species = "hsa",
                       kegg.native = FALSE, key_pos = key_pos,
                       out.suffix = title)
  }
  setwd(old_dir)
}

#' @title Create a Lolliplot.
#' @description Makes the LD dataframe using the Ensembl API for
#' the 1000 genomes project and a given population as
#' selected by the user.  The population has a code associated
#' with it.
#' @param rsid: rsid to compare to population for LD.
#' @param pop_code: Population code given and outlines by Ensembl and the
#' 1000 genomes project conventions.  This can be found with the
#' get_LD_populations function.
#' @param d_prime_thresh: The D' alpha value that the dataset can be
#' resticted by.
#' @param xaxis_show: Allows for the positions to be turned off if they are
#' overlapping with the Genes.
#' @param tableOff: Allows for the table of variants that are linked but not
#' able to be graphed to be shown on the figure on not. Sometimes the figure
#' holds too many variants and overlaps so it may be best to turn it off.
#' @returns Dataframe with LD scores of a specific population.
lolliplot_LD <- function(get_output, title = NULL, pdf_out = FALSE,
                         xaxis_show = FALSE, tableOff = TRUE) {
  # 1. Check to make sure the dataframe input isn't empty.
  .check_dataframe(get_output)
  # 1.1 Check if the title is null or na.
  if(is.null(title) || is.na(title)) {
    title = "Vargen Sample Plot Title"
  }
  # 1.2 Check if the pdf_out value is a boolean
  if(is.null(pdf_out) || is.na(pdf_out)) {
    pdf_out <- FALSE
  }
  # 1.3 This needs to be removed as it creates unnessary
  # categories if it is included.
  get_output$d_prime <- NULL
  
  # 2. Make a Grange object.
  # 2.1 get the rsids with only their numeric parts.
  get_output$"num_variation2" <- as.numeric(gsub("rs([0-9]+).*",
                                                 "\\1",
                                                 get_output[["variation2"]]))
  get_output <- get_output[order(get_output[["num_variation2"]]),]
  
  # 3. Get the variables needed for the GRanges.
  # 3A. The chromosome information and genes names.
  snp_mart <- biomaRt::useEnsembl(biomart = "snp",
                                  host = "www.ensembl.org",
                                  version = 100,
                                  dataset = "hsapiens_snp")
  
  rsids_info <- biomaRt::getBM(
    attributes = c("refsnp_id", "chr_name", "ensembl_gene_stable_id",
                   "chrom_start", "chrom_end"),
    filters = "snp_filter",
    values = unlist(get_output[["variation2"]]),
    mart = snp_mart, uniqueRows = TRUE)
  
  # 3B.1 Removing variants without any stable gene id and adding them into a separate
  # dataframe to be added to the figure after.
  rsids_unknown <- subset(rsids_info, rsids_info[["ensembl_gene_stable_id"]] == "")
  rsids_info <- subset(rsids_info, rsids_info[["ensembl_gene_stable_id"]] != "")
  
  chrom <- unique(rsids_info[["chr_name"]])
  # 3C. Get the start and stop position for the genes as well as their symbol.
  gene_info <- biomaRt::getBM(
    attributes = c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"),
    filters = "ensembl_gene_id",
    values = rsids_info[["ensembl_gene_stable_id"]],
    mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
  
  cnt <- 1
  rep_col <- c()
  rep_index <- which(gene_info[["hgnc_symbol"]] %in% "")
  for(item in rep_index) {
    gene_info[["hgnc_symbol"]][[item]] <- paste("Uknown", cnt, sep = "")
    cnt <- cnt + 1
  }
  
  colnames(rsids_info) <- c("refsnp_id", "chr_name", "ensembl_gene_id",
                            "start_position_rsid","end_position_rsid")
  all_info <- merge(rsids_info, gene_info, by = "ensembl_gene_id", all.x = TRUE)
  all_info$"dup" <- duplicated(all_info[["refsnp_id"]])
  all_info <- subset(all_info, all_info[["dup"]] == FALSE)
  
  res_all <- unique(all_info[,c(6,7,8)])
  
  genes <- unlist(res_all[["hgnc_symbol"]])
  genes <- genes[!is.na(genes)]
  
  # 3D. The variant IDs and names.
  rsid_pos <- unique(unlist(all_info[["start_position_rsid"]]))
  rsid <- unique(unlist(all_info[["refsnp_id"]]))
  
  # 3E. The width of each gene.
  gene_width <- unlist(unique(res_all[["start_position"]]))
  
  # 4. Get all the information into the same dataset restricted by the origional data
  get_output <- as.data.frame(lapply(get_output, unlist))
  colnames(get_output) <- gsub("variation2", "refsnp_id", colnames(get_output))
  all_info <- merge(get_output, all_info, by = "refsnp_id", all.x = TRUE)
  
  unknown <- subset(all_info, is.na(all_info[["chr_name"]]))
  unknown <- merge(rsids_unknown, unknown, by = "refsnp_id", all = TRUE)
  
  # This will contain all information that has a gene and positional information
  # needed for graphing.
  all_info <- subset(all_info, is.na(all_info[["chr_name"]]) == FALSE)
  
  rsid_unknown <- unlist(unknown[["refsnp_id"]])
  
  rsid_pos_unknown <- c()
  multi <- 10
  for(i in 1:length(rsid_unknown)) {
    rsid_pos_unknown <- c(rsid_pos_unknown, multi)
    multi <- multi + 10
  }
  
  # 5. Get the total width of the gene in a column with the start and end right after the other.
  all_info <- all_info[order(all_info[["start_position"]]),]
  
  height_score <- all_info[["r2"]] * 100
  gene_length_height <- replicate(length(genes), .018)
  
  height_score_unknown <- unknown[["r2"]] * 100
  gene_length_height_unknown <- .018
  
  start <- unlist(unique(all_info[["start_position"]]))
  stop <- unlist(unique(all_info[["end_position"]]))
  diff <- stop - start
  
  start <- start[!is.na(start)]
  stop <- stop[!is.na(stop)]
  diff <- diff[!is.na(diff)]
  
  # 6. Make the scores much like a heatmap for the gradiant.
  gradient <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
                "#4292C6", "#2171B5", "#08519C", "#08306B")
  all_info$"color_score" <- "#F7FBFF"
  r2 <- unlist(all_info[["r2"]])
  rep_col <- c()
  for(value in r2) {
    if(value*100 <= 100/length(gradient)) {
      rep_col <- c(rep_col, "#F7FBFF")
    }
    if(value*100 > 100/length(gradient) && value*100 <= 100/length(gradient)*2) {
      rep_col <- c(rep_col, "#DEEBF7")
    }
    if(value*100 > 100/length(gradient)*2 && value*100 <= 100/length(gradient)*3) {
      rep_col <- c(rep_col, "#C6DBEF")
    }
    if(value*100 > 100/length(gradient)*3 && value*100 <= 100/length(gradient)*4) {
      rep_col <- c(rep_col, "#9ECAE1")
    }
    if(value*100 > 100/length(gradient)*4 && value*100 <= 100/length(gradient)*5) {
      rep_col <- c(rep_col, "#6BAED6")
    }
    if(value*100 > 100/length(gradient)*5 && value*100 <= 100/length(gradient)*6) {
      rep_col <- c(rep_col, "#4292C6")
    }
    if(value*100 > 100/length(gradient)*6 && value*100 <= 100/length(gradient)*7) {
      rep_col <- c(rep_col, "#2171B5")
    }
    if(value*100 > 100/length(gradient)*7 && value*100 <= 100/length(gradient)*8) {
      rep_col <- c(rep_col, "#08519C")
    }
    if(value*100 > 100/length(gradient)*8 && value*100 <= 100/length(gradient)*9) {
      rep_col <- c(rep_col, "#08306B")
    }
  }
  all_info[["color_score"]] <- rep_col
  
  all_info <- all_info[order(all_info[["start_position"]]),]
  genes <- unique(unlist(all_info[["hgnc_symbol"]]))
  genes <- genes[!is.na(genes)]
  
  # 7. Make the GRanges and plot the track.
  lolliplot_vargen <- GenomicRanges::GRanges(chrom, IRanges::IRanges(rsid_pos, width = 1, names = rsid))
  lolliplot_vargen_unknown <- GenomicRanges::GRanges(chrom, IRanges::IRanges(rsid_pos_unknown, width = 1, names = rsid_unknown))
  
  features <- GenomicRanges::GRanges(chrom, IRanges::IRanges(c(start),
                                                             width = c(diff),
                                                             names = genes))
  
  features_unknown <- GenomicRanges::GRanges(chrom, IRanges::IRanges(c(1, max(rsid_pos_unknown)+10)))
  # 8. Get the colors for the sections, ends, and the caps.
  colors <- colorspace::qualitative_hcl((length(genes)), palette = "Set 3")
  lolliplot_vargen$dashline.col <- lolliplot_vargen$color
  
  lolliplot_vargen$color <- all_info[["color_score"]]
  lolliplot_vargen_unknown$color <- unknown[["color_score"]]
  
  temp <- cbind(genes, colors)
  colnames(temp) <- c("hgnc_symbol", "colors")
  all_info <- merge(all_info, temp, by = "hgnc_symbol")
  all_info <- all_info[order(all_info[["start_position"]]),]
  
  label.parameter.gp.col <- grid::gpar(col = all_info[["colors"]])
  lolliplot_vargen$label.parameter.gp <- label.parameter.gp.col
  
  # 9. Change the height.
  lolliplot_vargen$score <- height_score
  features$height <- gene_length_height
  lolliplot_vargen_unknown$score <- height_score_unknown
  features_unknown$height <- gene_length_height_unknown
  
  # 10. Get the overlap of the genes and make separate layers for each one by finding
  # if the start and stop lengths have an overlap.  If they do, then they are
  # separated into different layers.
  pos <- c() # Holds the positions that are overlaped.
  ranges <- split(IRanges::IRanges(all_info[["start_position"]],
                                   all_info[["end_position"]]),
                  all_info[["hgnc_symbol"]])
  
  for(gr1 in ranges) {
    for(gr2 in ranges) {
      overlap <- IRanges::findOverlaps(gr1, gr2)
      if(length(overlap) != 0 && gr1 != gr2) {
        temp <- data.frame(gr2)
        pos <- c(pos, temp[1,1])
        temp <- data.frame(gr1)
        pos <- c(pos, temp[1,1])
      }
    }
  }
  pos <- unique(pos)
  # Holds the unique genes that are overlaped.
  gene_overlap <- unique(subset(all_info,
                                all_info[["start_position"]] == pos)[,c("hgnc_symbol")])
  
  # 11. Make the GRanges with each of the genes that overlap not being on the same layers.
  features <- GenomicRanges::GRanges(c(seqnames = NULL, ranges = NULL, strand = NULL)) # Holds the feature layers
  genes <- genes[!(genes %in% gene_overlap)]
  colors <- c()
  layerIDs <- c()
  for(gene in gene_overlap) {
    # Add gene to the genes vector
    genes <- c(genes, gene)
    # Subset the information
    all_info_ss <- all_info[all_info[["hgnc_symbol"]] %in% genes, ]
    start <- unlist(unique(all_info_ss[["start_position"]]))
    stop <- unlist(unique(all_info_ss[["end_position"]]))
    diff <- stop - start
    
    genes <- unique(all_info_ss[["hgnc_symbol"]])
    colors <- c(colors, unique(all_info_ss[["colors"]]))
    
    feature <- GenomicRanges::GRanges(chrom, IRanges::IRanges(c(start),
                                                              width = c(diff),
                                                              names = genes))
    # Add features to the vector
    features <- c(features, feature)
    layerIDs <- c(layerIDs, unique(all_info_ss[["hgnc_symbol"]]))
  }
  
  features$fill <- c(colors)
  
  # 12. Make the table for the Unknown variants that are still linking. These
  # Variants do not have positional data usable in a manner on this plot and as
  # such should not be represented in the positional graph. They are instead
  # put into a table with a separate warning for the user and displayed below.
  unknown <- cbind(unknown[1], unknown[2], unknown[4], unknown[5], unknown[8])
  colnames(unknown) <- c("SNP", "Chomosome", "Gene_Start", "Gene_Stop", "r2")
  unknown$"Reason" <- "*** No position was found in data for SNP on the gene ***"
  # Make a table Grob with the gridExtra package.
  table <- gridExtra::tableGrob(unknown)
  
  print(start)
  print(stop)
  # 12. Change the axis to fit the data.
  xaxis <- c(start, (start[length(start)] + diff[length(diff)]))
  if(is.null(layerIDs) == FALSE) {
    features$featureLayerID <- paste("tx", layerIDs, sep = "_")
    # 13. Make the file for the output of the pdf if the boolean is True. The
    # output which the file will be saved is the current working directory.
    if(tableOff == FALSE){
      if(pdf_out == TRUE) {
        pdf(file = paste0(getwd(), "lolliplot_vargen", sep = ""), paper = "a4")
        
        if(xaxis_show == TRUE) {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = xaxis)
        } else {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = FALSE)
        }
        grid::grid.text(title, x = .5, y = .98, just = "top",
                        gp = grid::gpar(cex = 1.5, fontface = "bold"))
        lolli <- grid::grid.grab()
        cowplot::plot_grid(lolli, table)
        dev.off()
      } else {
        print("hit3")
        if(xaxis_show == TRUE) {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = xaxis)
        } else {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = FALSE)
        }
        grid::grid.text(title, x = .5, y = .98, just = "top",
                        gp = grid::gpar(cex = 1.5, fontface = "bold"))
        
        lolli <- grid::grid.grab()
        cowplot::plot_grid(lolli, table, ncol = 1)
      }
    } else {
      if(pdf_out == TRUE) {
        pdf(file = paste0(getwd(), "lolliplot_vargen", sep = ""), paper = "a4")
        
        trackViewer::lolliplot(lolliplot_vargen, features,
                               legend = legend, xaxis = xaxis)
        dev.off()
      } else {
        if(xaxis_show == TRUE) {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = xaxis)
        } else {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = FALSE)
        }
      }
    }
  } else {
    start <- unlist(unique(all_info[["start_position"]]))
    stop <- unlist(unique(all_info[["end_position"]]))
    genes <- unlist(res_all[["hgnc_symbol"]])
    genes <- genes[!is.na(genes)]
    
    features <- GenomicRanges::GRanges(chrom, IRanges::IRanges(c(start),
                                                               width = c(diff),
                                                               names = genes))
    colors <- colorspace::qualitative_hcl((length(genes)), palette = "Set 3")
    features$fill <- c(colors)
    lolliplot_vargen$dashline.col <- lolliplot_vargen$color
    
    if(tableOff == FALSE){
      if(pdf_out == TRUE) {
        pdf(file = paste0(getwd(), "lolliplot_vargen", sep = ""), paper = "a4")
        
        if(xaxis_show == TRUE) {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = xaxis)
        } else {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = FALSE)
        }
        grid::grid.text(title, x = .5, y = .98, just = "top",
                        gp = grid::gpar(cex = 1.5, fontface = "bold"))
        lolli <- grid::grid.grab()
        cowplot::plot_grid(lolli, table)
        dev.off()
      } else {
        if(xaxis_show == TRUE) {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = xaxis)
        } else {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = FALSE)
        }
        grid::grid.text(title, x = .5, y = .98, just = "top",
                        gp = grid::gpar(cex = 1.5, fontface = "bold"))
        
        lolli <- grid::grid.grab()
        cowplot::plot_grid(lolli, table, ncol = 1)
      }
    } else {
      if(pdf_out == TRUE) {
        pdf(file = paste0(getwd(), "lolliplot_vargen", sep = ""), paper = "a4")
        
        trackViewer::lolliplot(lolliplot_vargen, features,
                               legend = legend, xaxis = xaxis)
        dev.off()
      } else {
        if(xaxis_show == TRUE) {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = xaxis)
        } else {
          trackViewer::lolliplot(lolliplot_vargen, features,
                                 legend = legend, xaxis = FALSE)
        }
      }
    }
  }
}