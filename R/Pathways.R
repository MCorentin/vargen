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
.check_kegg_validity <- function(params) {
  print("Checking traits to see if they are valid")
  # 1. Check the directory exists for the VarGen data.
  if(dir.exists(params[[1]]) == FALSE) {
    stop(paste0("The directory used for vargen_dir: ", params[[1]],
                " was not found in the current directory ", getwd()))
  }

  # 2. Check that the gene mode is valid.
  if((!(params[[7]] %in% c(0,1,2)) && is.na(params[[7]]) == FALSE && is.na(params[[6]])) ||
     (!(params[[7]] %in% c(0,1,2)) && is.na(params[[7]]) == FALSE && is.na(params[[6]]) == FALSE)) {
    # Case 1: If the given gene mode is not in one of the valid selection types.
    print(paste0("Gene_mode: ", params[[7]],
                 " does not exist so the default mode of 2 is being used instead."))
    params[[7]] <- 2
  } else if(params[[7]] %in% c(0,1) && is.na(params[[7]]) == FALSE && is.na(params[[6]])) {
    # Case 2: If the given gene mode is in one of the valid selection types, but no genes are given.
    # Only 1 and 0 are selected here as the default is 2 in the program.
    print("No genes were given with gene mode, thus no gene mode can/will be used")
  }

  # 3. Check that the p-value is valid.
  if(typeof(params[[9]]) != "double" && is.na(params[[9]]) == FALSE) {
    if(is.na(as.double(params[[9]]))){
      print(paste0("p-value threshold could not be converted to a double from the given value: '",
                   params[[9]], "'. No pvalue threshold will be used."))
      params[[9]] <- NA
    }
  }

  # 4. Check if the chromosome/s are valid and if they are not attempt to fix them.
  if(length(grep("chr[0-9]", x = params[[4]])) != length(params[[4]]) && is.na(params[[4]]) == FALSE){
    #Separate the bad and good chromosomes
    bad_chrs <- params[[4]][!(params[[4]] %in% params[[4]][grep("chr[0-9]", x = params[[4]])])]
    good_chrs <- params[[4]][grep("chr[0-9]", x = params[[4]])]

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
      params[[4]] <- c(good_chrs, bad_chrs)
      out_print <- paste(params[[4]], collapse = '')
      out_print <- gsub("(\\d+)\\c", "\\1\\ \\c", out_print)
      print(paste0("These are the remaining chromosomes after fixing: ", out_print))
    } else {
      print(paste0("Chromosomes must be given as a value 'chr3'. Currently given badly as: ",
                   bad_chrs, ". All chromosomes to be used by default"))
    }
  }

  # 5. Check if the output directory exists and if it doesn't then make it.
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
        print(paste0("Using output directory ", params[[2]]))
      }
    } else {
      print(paste0("Using output directory ", params[[2]]))
    }
  }

  # 6. If a rsid is given then search just by the rsid and not by any other params even
  # if they have been given. A default title should made for each rsid as well and
  # not just the default; This is done later on.
  if(is.null(params[[8]])) {
    params[[8]] = NA
  } else if(is.na(params[[8]]) == FALSE) {
    params[[3]] = NA
    params[[6]] = NA
    params[[4]] = NA
  }

  # 7 .Check the OMIM IDs
  if(is.null(params[[14]])) {
    params[[14]] = NA
  } else if(is.null(params[[14]]) == FALSE && is.null(params[[9]] == FALSE)) {
    print("No p-value threshold restiction available for OMIM searching.")
    params[[9]] = NA
    params[[13]] = 1
  }

  # 8. Check that the colors are properly set.
  if(is.null(params[[10]]) || is.na(params[[10]])){
    params[[10]] = "#da8cde"
  }
  if(is.null(params[[12]]) || is.na(params[[12]])){
    params[[12]] = "#1254c7"
  }
  if(is.null(params[[11]]) || is.na(params[[11]])){
    params[[11]] = "#09873e"
  }

  # 9. Check that the show_pval_thresh is a boolean
  if(params[[13]] != 1 && params[[13]] != 0) {
    print(paste0("parameter: ", params[[13]], " is not valid. Set to 1"))
    params[[13]] = 1
  }

  # 10. Check that traits not given as null and if so change to be NA.
  if(is.null(params[[3]])) {
    params[[3]] = NA
  }

  return(params)
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

  return(omim_info)
}

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

#' @title Makes output directory for Windows, Linux, or Unix system.
#' @description Makes a directory if none is given or returns the
#' directory of the computer system so that figures can be placed
#' in a given directory.
#' @param dir: The path name that should be made.
#' @param name: The name for the final directory.
#' @return a the directory or make a directory is none is given.
.make_dir <- function(dir, name = ".") {
  if(is.na(dir[[1]])  || dir[[1]] == " " || is.double(dir[[1]])) {
    # Make the proper output directory if one isn't given.
    if(.Platform$OS.type == "windows") {
      dir.create(name)
      output_dir <- paste0(getwd(), "\\", name, "\\")
    } else {
      dir.create(name)
      output_dir <- paste0(getwd(), "/", name, "/")
    }
  } else {
    if(.Platform$OS.type == "windows") {
      output_dir <- paste0(dir[[1]], "\\", name, "\\")
    } else {
      output_dir <- paste0(dir[[1]], "/", name, "/")
    }
  }
  return(output_dir)
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
                                   ko_matched){
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

#---- Helpers ----

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
#' @param pval_thresh: The cut off threshold for the p-values of the variants.
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
#' @param omim_ids: The omim id or ids to be searched only in a KEGG pathway.
#' @param rsids: A single snp rsid identifier or a vector of identifiers which
#' will return a list of found information from the gwas catalog in a dataframe.
#' @param genes: The mapped and reported genes that are found in the file.
#' @param gene_mode: The mode for selecting if the reported genes and the
#' mapped genes should both be search for the given genes or if only one or the
#' other should be.  0 is for mapped genes only, 1 is for only reported, 2 is
#' for both. The defualt is 2.
#' @param pval_thresh: The cut off threshold for the p-values of the variants.
#' @return Nothing; The KEGG pathways figures in a given or default directory.
kegg_graph <- function(vargen_dir, output_dir = NA, traits = NA,
                       chrs = NA, title = NA, genes = NA, gene_mode = 2,
                       rsids = NA, omim_ids = NA, pval_thresh = NA,
                       color = "#da8cde", pval_low_col = "#09873e",
                       pval_high_col = "#1254c7", pval_thresh_show = 1) {

  # 1. Check if the parameters are valid.
  params = .check_kegg_validity(list(vargen_dir, output_dir, traits, chrs, title,
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
      kegg_paths <- biomaRt::getBM(
        attributes = c("kegg_enzyme", "ensembl_gene_id"),
        filters = c("ensembl_gene_id"),
        values = variants_traits$"SNP_GENE_IDS",
        mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
    } else {
      kegg_paths <- biomaRt::getBM(
        attributes = c("kegg_enzyme", "ensembl_gene_id"),
        filters = c("ensembl_gene_id"),
        values = variants_traits$"SNP_GENE_IDS",
        mart = connect_to_gene_ensembl(), uniqueRows = TRUE)
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
    # 6.4 Restrict it by p-value threshold
    if(is.na(pval_thresh) == FALSE){
      # Column must be changed to be numeric else it's not calculated properly.
      kegg_paths[["p_value"]] <- as.numeric(as.character(kegg_paths[["p_value"]]))
      if(params[[13]] == 0){
          # If the user wants to show the p-value threshold colors that are below
          # and above the threshold.
          kegg_paths_high <- kegg_paths[which(kegg_paths[["p_value"]] >= params[[9]]), ]
          kegg_paths_low <- kegg_paths[which(kegg_paths[["p_value"]] < params[[9]]), ]
      } else {
          kegg_paths <- kegg_paths[which(kegg_paths[["p_value"]] < params[[9]]), ]
      }
    }

    if(nrow(kegg_paths) == 0){
      writeLines(c(variant_info[["rsid"]][cnt]), file_out)
      print("No KEGG pathways available for this search request. Try broading search if possible.")
      print(paste0("File written in"), getwd())
    }

    # 7. Output KEGG pathways for each pathway in the output firectory
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
                                             ko_matched)

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
#' @param threshold: The number of the top most found values in the dataset.
#' @return Nothing; The KEGG pathways figures in a given or made directory.
pathview_vargen <- function(data, output_dir = NA, title = NULL, top_thresh = NULL,
                            key_pos = "topright") {
  # 1. Check input parameters to see if they are valid.
  if(is.na(data) || is.null(data)) {
    stop(paste0("No data found inside the given dataset from the VarGen pipeline.",
                " The data frame is either NA or null."), sep = "")
  } else if(nrow(data) == 0) {
    stop(paste("No data found inside the given dataset from the VarGen pipeline.",
               " The data frame has 0 rows."), sep = "")
  }
  # 1.2 Check if the output directory exists and if it doesn't then make it.
  if(is.na(output_dir) == FALSE) {
    if(dir.exists(output_dir) == FALSE) {
      # Make the proper output directory if one isn't given.
      if(.Platform$OS.type == "windows"){
        if(dir.exists(paste0(getwd(), "\\", "pathview\\")) == FALSE){
          dir.create("pathview")
        }
        output_dir <- paste0(getwd(), "\\", "pathview\\")
      } else {
        if(dir.exists(paste0(getwd(), "\\", "pathview\\")) == FALSE){
          dir.create("pathview")
        }
        output_dir <- paste0(getwd(), "\\", "pathview\\")
        print(paste0("Using output directory ", output_dir))
      }
    } else {
      print(paste0("Using output directory ", output_dir))
    }
  }
  output_dir <- .make_dir(output_dir, "pathview")

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
    variants_info <- subset(variants_info, variants_info[["ensembl_gene_id"]] != "")
    # 3.2 Remove those without an Entrez ID as they can't be used.
    variants_info <- subset(variants_info, is.na(variants_info[["entrezgene_id"]]) == FALSE)
    # 3.3 Remove those without a Kegg pathway as they also can't be used.
    variants_info <- subset(variants_info, is.na(variants_info[["kegg_enzyme"]]) == FALSE)
    variants_info <- subset(variants_info, variants_info[["kegg_enzyme"]] != "")
    variants_info <- subset(variants_info,
                            !duplicated(subset(variants_info,
                                               select=c("kegg_enzyme", "hgnc_symbol"))))
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
  sum <- sum(data.frame(table(data[["hgnc_symbol"]]))[["Freq"]])
  for(i in 1:nrow(hgnc_cnt)) {
    hgnc_cnt[["cnts"]][[i]]  <- hgnc_cnt[["cnts"]][[i]] / sum
  }

  # 6. Get the restricted information since pathways are the limiting factor.
  hgnc_cnt <- merge(hgnc_cnt, variants_info, by = "hgnc_symbol", all = TRUE)

  # 6. Make the proper input for each of the pathways for the pathview and then create the Pathview images.
  # The gene data is the frequency of the counts of variants found by VarGen with entrez genes IDs as the
  # requested ID.  The pathway id is the kedd_id.
  old_dir = getwd()
  for(kegg_pathway in hgnc_cnt[["kegg_enzyme"]]) {
    kegg_id <- unlist(strsplit(kegg_pathway, "\\+"))[1]
    sel_info <- unique(hgnc_cnt[, c(6,2)])
    sel_genes <- as.matrix(sel_info[, 2])
    rownames(sel_genes) <- c(sel_info[, 1])

    setwd(output_dir)
    pv.out <- pathview::pathview(gene.data = sel_genes, pathway.id = kegg_id,
                                 gene.idtype = "entrez",
                                 kegg.native = FALSE, key_pos = key_pos,
                                 out.suffix = title)
  }
  setwd(old_dir)
}

#' @title Visualize the online OMIM Graphs.
#' @description Takes the user to the OMIM website's graphic for a given OMIM
#' ID. Users can choose either radial or linear as a mode - with the default
#' being linear.
#' @param type: the type of graph to go to either as "linear" or "radial". Can
#' also be given as 'l' or 'r'.
#' @param omim_id: A single OMIM ID given as a string.
#' @return Nothing; Takes user to url of omim id selected.
omim_graph_online <- function(type = "linear", omim_id){
  # 1. Check the type to see if it is valid then check if the omim id is not null or na.
  if(!(type %in% c("linear", "l", "Linear", "r", "radial", "Radial"))) {
    stop(print(paste0("There is no option '", type,
                      "' for parameter type. Type options are 'linear' or 'radial'.")))
  }
  if(is.null(omim_id) || is.na(omim_id)) {
    stop(print("There is no OMIM ID given."))
  }

  # 2. Get the url from the OMIM database depending on the omim id.
  if(type == "linear") {
    url <- paste0("https://omim.org/graph/linear/", omim_id, sep = "")
    browseURL(url)
  } else {
    url <- paste0("https://omim.org/graph/radial/", omim_id, sep = "")
    browseURL(url)
  }
}

#' @title Makes a treemap from the VarGen pipeline.
#' @description Makes a treemap from the vargen pipeline.
#' @param data: The output from the vargen pipeline.
#' @param title: The title of the figure.
#' @param threshold: The number of the top most found values in the dataset.
#' @param mode: Sets the counts to be either the raw or normalized counts
#' (# variants / kbp). The modes are 'raw' for the raw counts, and 'norm' for
#' normalized counts.
#' @return Nothing; A treemap graph.
treemap_graph <- function(data, title = "VarGen Treemap Graph", top_thresh = NULL, mode = "raw") {
  # 1. Get the count data.
  hgnc_cnt <- data.frame(table(data[["hgnc_symbol"]]))
  colnames(hgnc_cnt) <- c("hgnc_symbol", "cnts")

  if(mode == "norm") {
    # 2. Get the gene length data for normalization.
    length_info <- get_omim_genes(omim_ids = data[["trait"]], connect_to_gene_ensembl())
    data <- merge(data[3:7], length_info, by = c("ensembl_gene_id", "hgnc_symbol"),
                         all.x = TRUE)
    data <- data[!is.na(data[["end_position"]]), ]
    data <- unique(data)

    # 3. Normalize the count data.
    norm_hgnc <- merge(hgnc_cnt, data[, c(2,7,8)], by = c("hgnc_symbol"),
                       all.x = TRUE)
    norm_hgnc <- unique(norm_hgnc)
    norm_hgnc$"len" <- (norm_hgnc[["end_position"]] - norm_hgnc[["start_position"]])
    norm_hgnc$"cnts" <- (norm_hgnc[["cnts"]] / norm_hgnc[["len"]])
    hgnc_cnt <- norm_hgnc
  }

  if(is.null(top_thresh) == FALSE) {
    hgnc_cnt <- hgnc_cnt[order(-hgnc_cnt[["cnts"]]),]
    hgnc_cnt <- hgnc_cnt[1:top_thresh,]
  }

  hgnc_cnt$"labels" <- paste(hgnc_cnt[["hgnc_symbol"]], hgnc_cnt[["cnts"]], sep = ": ")

  x11()
  treemap(dtf = hgnc_cnt, index = "labels", vSize = "cnts",
           type = "value", vColor = "cnts", title = title)
}
