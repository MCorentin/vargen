# VarGen

VarGen is an R package designed to get a list of variants related to a disease. It just need an OMIM morbid ID as input 
and optionally a list of tissues / gwas traits of interest to complete the results. You can use your own customised list 
of genes instead of an OMIM ID. VarGen is also capable of annotating the variants to help you rank and identify the most 
impactful ones.

All the coordinates are based on the hg38 version of the human genome.

## Table of Contents

- [VarGen](#vargen)
- [Table of Contents](#table-of-contents)
- [Workflows](#workflows)
    - [VarGen workflow](#vargen-workflow)
    - [VarPhen workflow](#varphen-workflow)
- [Installation](#installation)
    - [Dependencies](#dependencies)
    - [Install VarGen with devtools](#install-vargen-with-devtools)
    - [Install VarGen from source](#install-vargen-from-source)
    - [Database files](#database-files) 
- [How to use VarGen](#how-to-use-vargen)
    - [Examples](#examples)
- [How to use VarPhen](#how-to-use-varphen)
    - [Example](#example)
- [Tips](#Tips)
    - [How to get the OMIM morbid ID](#how-to-get-the-omim-morbid-id)
    - [How to use a local gwas catalog file](#how-to-use-a-local-gwas-catalog-file)
    - [How to list the available tissues in GTEx](#how-to-list-the-available-tissues-in-gtex)
    - [How to list the available GWAS traits](#how-to-list-the-available-gwas-traits)
    - [How to use a custom list of genes](#how-to-use-a-custom-list-of-genes)
    - [How to annotate the variants](#how-to-annotate-the-variants)
    - [How to plot the gwas variants](#how-to-plot-the-gwas-variants)
    - [How to plot the omim variants](#how-to-plot-the-omim-variants)
    - [How to filter VarGen output](#how-to-filter-vargen-output)

## Workflows

### VarGen workflow

This pipeline is centred on the genes linked to the disease of interest in OMIM (subsequently called the "OMIM genes").

VarGen outputs variants from the following sources:
- **OMIM:** Variants located directly on the "OMIM genes".
- **FANTOM5:** Variants located on the enhancers / promoters of the "OMIM genes".
- **GTEx:** Variants associated with a change in expression for the "OMIM genes", in certain tissues. Currently GTEx v7 
and v8 are supported.
- **GWAS catalog:** Variants associated with the phenotype of interest.

![VarGen workflow](./images/VarGen_workflow.png?raw=true)

The variants are then annotated with [myvariant](http://bioconductor.org/packages/release/bioc/html/myvariant.html "myvariant bioconductor package"). 

This pipeline is designed as a discovery analysis, to identify potential new variants, **you should not expect every variants from the pipeline to have an effect on the phenotype**. 
The annotation will help you defining which variants to keep or discard. The annotation contains the [CADD Phred score](https://cadd.gs.washington.edu/ "CADD main page"),
annotation type (eg: "Intergenic"), consequence (eg: "DOWNSTREAM"), [clinvar clinical significance](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/ "Representation of clinical significance in ClinVar and other variation resources at NCBI") and [snpEff impact](http://snpeff.sourceforge.net/SnpEff_manual.html "snpEff Manual").

### VarPhen workflow

An alternative pipeline is available as part of this package, called "VarPhen", it outputs a smaller list of variants, but 
directly related to the disease of interest. It relies on biomaRt to link variants to phenotypes.

![VarPhen workflow](./images/VarPhen_workflow.png?raw=true)

## Installation

VarGen is available on [GitHub](https://github.com/MCorentin/VarGen "GitHub VarGen page")

### Dependencies

VarGen needs the following:
- **R** (tested on version 3.6)
- **An internet connection**
- **The following R libraries:** (The number is the version tested during development)
```` 
    Bioconductor (3.9)      biomaRt (2.40.3)        gtools (3.8.1)         
    gwascat (2.16.0)        jsonlite (1.6)          GenomeInfoDb (1.20.0)
    IRanges (2.18.1)        httr (1.4.1)            BiocGenerics (0.30.0)
    stringr (1.4.0)         utils (3.6.0)           splitstackshape (1.4.8)
    ggplot2 (3.2.1)         rtracklayer (1.44.2)    BiocManager (1.30.4)
    R.utils (2.9.0)         myvariant (1.14.0)      GenomicRanges (1.36.0)
````

- Optional R libraries (used for the visualisation functions: "vargen_visualisation" and "plot_manhattan_gwas"):
````
    Gviz (1.28.1)           ggbio (>= 1.32.0)       grDevices (>= 3.6.0)
````

To install the dependencies you can use the following command in R (it might take a while depending on your connection):
````
# If not already installed
install.packages("BiocManager")

BiocManager::install(c("biomaRt", "gtools", "GenomicRanges", "gwascat", "jsonlite",
                       "GenomeInfoDb", "IRanges", "httr", "BiocGenerics", "stringr",
                       "utils", "splitstackshape", "ggplot2", "rtracklayer",
                       "R.utils", "myvariant"), dependencies = TRUE)
````
**note:** "R.methodsS3" and "R.oo" will be installed as dependencies of "R.utils"

### Install VarGen with devtools

The easiest way to get VarGen is to install it directly from R using “devtools”:
To install devtools:
````
install.packages("devtools")
library(devtools)
install_github("MCorentin/VarGen")
library(VarGen)
````

### Install VarGen from source

Alternatively you can clone the GitHub repository:
````
git clone https://github.com/MCorentin/VarGen
````
Then open R and install the VarGen.R script from source:
````
library(utils)
install.packages("./VarGen/", repos = NULL, type = "source")
````

### Database files

VarGen is fetching data from public databases, it needs the following files (the files should all be in the same folder). 
The easiest way to get them is to use the *vargen_install* function from this package.

````
vargen_install(install_dir = "./vargen_data", gtex_version = "v8", verbose = T)
````

Alternatively, they can be installed manually:

 - __enhancer_tss_associations.bed__, this is the enhancer to transcript start site association file from FANTOM5. 
 (available at:  http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed)
 
 - __hg19ToHg38.over.chain__, some databases are still using information from the human reference genome "hg19". VarGen 
 will use this file to liftover the information to "hg38" (available at: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)

- __GTEx_Analysis_v8_eQTL__, a folder containing the significant variant gene pairs from the Genotype-Tissues Expression 
database (GTEx) (available at: https://gtexportal.org/home/datasets). v7 and v8 are supported by VarGen.
    
    Direct link for v7: https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz

    Direct link for v8: https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar 

## How to use VarGen

The main entry point is *vargen_pipeline*, this function will get a list of variants possibly related to the phenotype of 
interest. The only mandatory input is an OMIM morbid ID (or alternatively a list of gene IDs).

You can also add a list of tissues, if you are interested in variants that are affecting the expression of the genes (from 
GTEx) and a list of GWAS traits to get variants from the [GWAS catalog](https://www.ebi.ac.uk/gwas/ "gwas catalog main page").

The pipeline outputs variants with their positions, IDs and related gene (if applicable). It is then possible to annotate 
the results with the function *annotate_variants*.

Output after annotation:
![Output after annotation](./images/Example_output_VarGen.PNG?raw=true)

### Examples

In the following example we are assuming that the files described in "Installation" have been installed in the folder 
"./vargen_data", either manually or by running:
````
 vargen_install("./vargen_data")
````

 - Example 1: "Simple query for Type 1 Diabetes Mellitus"
````
T1DM_variants <- vargen_pipeline(vargen_dir = "./vargen_data",
                                 omim_morbid = "222100")
````

 - Example 2: "Type 1 Diabetes Mellitus with GTEx"
````
# First select the tissues of interest (here pancreas)
# it is possible to select more than one tissue (eg: c("pancreas", "adipose"))
pancreas_gtex <- select_gtex_tissues(gtex_dir = "./vargen_data/GTEx_Analysis_v8_eQTL/",
                                     tissues_query = "pancreas")

T1DM_variants <- vargen_pipeline(vargen_dir = "./vargen_data", 
                                 omim_morbid = "222100", 
                                 gtex_tissues = pancreas_gtex)
````

 - Example 3: "Obesity with GTEx, GWAS and annotation"
````
# First select the tissues of interest (here adipose) 
# it is possible to select more than one tissue (eg: c("pancreas", "adipose"))
adipose_gtex <- select_gtex_tissues(gtex_dir = "./vargen_data/GTEx_Analysis_v8_eQTL/", 
                                    tissues_query = "adipose")

# List the available gwas traits:
list_gwas_traits(keywords = "obesity")

# Select gwas traits of interest
gwas_obesity <- c("Obesity (extreme)", "Obesity-related traits", "Obesity", "Obesity (early onset extreme)")

Obesity_variants <- vargen_pipeline(vargen_dir = "./vargen_data",
                                    omim_morbid = "601665", 
                                    gtex_tissues = adipose_gtex,
                                    gwas_traits = gwas_obesity,
                                    fantom_corr = 0.25, 
                                    verbose = TRUE)

# Annotation of the variants obtained with the previous command:
Obesity_annotation <- annotate_variants(rsid = Obesity_variants$rsid)

# Merging the original output with the annotation:
Obesity_variants_annotated <- merge(Obesity_variants, Obesity_annotation)
View(Obesity_variants_annotated)
````

## How to use VarPhen

VarPhen is the alternative pipeline of this package. It is designed to be much more specific than "VarGen", it should give 
less variants but all associated to the phenotypes. It relies on biomaRt to link variants to phenotypes (cf Workflow).

### Example

Example with obesity:
````
# First connect to snp ensembl
snp_mart <- connect_to_snp_ensembl()

obesity_phens <- get_phenotype_terms(keyword = "obesity", 
                                     snp_mart = snp_mart)

obesity_varphen_snps <- get_variants_from_phenotypes(phenotypes = obesity_phens, 
                                                     snp_mart = snp_mart)

obesity_varphen_snps_annotated <- annotate_variants(obesity_varphen_snps$refsnp_id)
````

Output example:
![VarPhen output](./images/Example_output_VarPhen.PNG?raw=true)


## Tips 

### How to get the OMIM morbid ID

You can search on the Online Mendelian Inheritance in Man website (https://www.omim.org/) or use the function *list_omim_accessions*.
````
gene_mart <-  connect_to_gene_ensembl()
list_omim_accessions(gene_mart)
````

### How to use a local gwas catalog file

Depending on your connection, creating the lastest gwas catalog using the *makeCurrentGwascat* function can take a long time. 
You can instead download a gwas catalog file from http://www.ebi.ac.uk/gwas/api/search/downloads/alternative and give it as 
input for the VarGen pipeline. VarGen uses the name of the file to get the extract date, so it needs to be in the format: 

\[filename\]_r**YYYYY**-**MM**-**DD**.tsv

For example: "gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv" will be parsed to get the extract date "2019-07-30".

You can then use the **gwascat_file** option in the *vargen_pipeline* function.
````
# Getting the gwas variants will be faster if the "gwascat_file" option is used
T1DM_variants <- vargen_pipeline(vargen_dir = "./vargen_data/", 
                                 omim_morbid = "222100", 
                                 gwas_traits = "Type 1 diabetes", 
                                 gwascat_file = "./gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv")
````

### How to list the available tissues in GTEx

The function __select_gtex_tissues__ will list the files in the "./vargen_data/GTEx_Analysis_v8_eQTL" folder corresponding to the keyword 
entered as input.

Example:
````
select_gtex_tissues(gtex_dir = "./vargen_data/GTEx_Analysis_v8_eQTL/", tissues_query = c("pancreas", "adipose"))

[1] "./vargen_data/GTEx_Analysis_v8_eQTL/Pancreas.v8.signif_variant_gene_pairs.txt.gz"
[2] "./vargen_data/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.signif_variant_gene_pairs.txt.gz"    
[3] "./vargen_data/GTEx_Analysis_v8_eQTL/Adipose_Visceral_Omentum.v8.signif_variant_gene_pairs.txt.gz"
````

### How to list the available GWAS traits

The function __list_gwas_traits__ will list the available gwas traits. It builds a gwas object from "makeCurrentGwascat" 
and list the unique list of `DISEASE/TRAIT`.

Example:
````
list_gwas_traits(keyword = "Obesity")

 [1] "Obesity (extreme)"
 [2] "Obesity-related traits"
 [3] "Obesity"
 [4] "Obesity (early onset extreme)"
 [5] "Obesity in adult survivors of childhood cancer exposed to cranial radiation"
 [6] "Obesity in adult survivors of childhood cancer not exposed to cranial radiation"
 [7] "Bilirubin levels in extreme obesity"
 [8] "Type 2 diabetes (young onset) and obesity"
 [9] "Obesity and osteoporosis"
[10] "Hepatic lipid content in extreme obesity"
[11] "Hyperinsulinemia in obesity"
````
**note** The function also accept a gwas catalog file (cf "How to use a local gwas catalog file" above):
````
list_gwas_traits(keyword = "Obesity", 
                 gwascat_file = "./gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv")
````

### How to use a custom list of genes

The *vargen_custom* pipeline has been designed to accept a list of gene IDs instead of a OMIM morbid id. 

The pipeline will return the same information as with the omim query (ie: variants on the genes, on the enhancers, changing 
expression in certain tissues and from gwas analysis).

Example:
````
# First select the tissues of interest (here adipose) 
# it is possible to select more than one tissue (eg: c("pancreas", "adipose"))
pancreas_tissues <- select_gtex_tissues(gtex_dir = "./vargen_data/GTEx_Analysis_v8_eQTL/", 
                                        tissues_query = "pancreas")

custom_variants <- vargen_custom(vargen_dir = "./vargen_data/", 
                                 gene_ids = c("ENSG00000134242", "ENSG00000049768"),
                                 outdir = "./custom_test/", 
                                 gtex_tissues = pancreas_tissues,
                                 gwas_traits = c("Type 1 diabetes"))
````

### How to annotate the variants

As the VarGen pipeline outputs a lot of variants, it is necessary to prioritise them. The *annotate_variants* function has 
been designed to associate each variant to a consequence, snpEff impact, cadd phred score and clinical significance.

/!\ Due to the different transcripts for a same gene, some variants will appear more than once with a different annotation, 
this is expected. You can check the variant position on the different isoforms with the *vargen_visualisation* function.

````
# Get the variants:
obesity_vargen <- vargen_pipeline(vargen_dir = "./vargen_data/", 
                                  omim_morbid = "601665",
                                  outdir = "./vargen_obesity/")

# Annotate the variants:
obesity_ann <- annotate_variants(obesity_vargen$rsid, verbose = T)

# Then you will just have to merge the variants (they will merge automatically by rsids)
obesity_vargen_ann <- merge(obesity_vargen, obesity_ann)

# We advise you to write the variants in a file, so you will not have to run the pipeline again.
write.table(x = obesity_vargen_ann,
            file = "./vargen_obesity/vargen_variants_annotated.tsv", sep = "\t")

obesity_vargen_ann <- read.table(file = "./vargen_obesity/vargen_variants_annotated.tsv",
                                 header = T, sep = "\t", stringsAsFactors = F)
````

### How to plot the gwas variants

If you want to visualise the variants in a manhattan plot, you can use the *plot_manhattan_gwas* function:
````
gwas_cat <- create_gwas(gwascat_file = "./gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv")
 
plot_manhattan_gwas(gwas_cat = gwas_cat, traits = c("Type 1 diabetes", "Type 2 diabetes"))

# Optional: if you want to save the plot as a pdf
grDevices::dev.print(pdf, "./manhanttan_diabetes.pdf")
````

![Example of manhanttan plot for diabetes](./images/manhattan_diabetes.png?raw=true)

The two thresholds, suggestive and significant, correspond to the definition given by Lander and Kruglyak:
````
Lander E, Kruglyak L. Genetic dissection of complex traits: guidelines for interpreting and reporting linkage results. Nat Genet. 1995;11(3):241-7.
````

### How to plot the omim variants

You can use the *vargen_visualisation* function to have an overview of the variants located on the omim genes. Note that 
one variant can be represented multiple times, as it consequence and phred score will be different according to each transcript.

````
# cf: 'How to annotate the variants' to create "vargen_variants_annotated.tsv":
obesity_vargen_ann <- read.table("vargen_obesity/vargen_variants_annotated.tsv")
 
gene_mart <- connect_to_gene_ensembl()
vargen_visualisation(annotated_snps = obesity_vargen_ann,
                     outdir = "./obesity_gviz/", 
                     device = "png", 
                     gene_mart = gene_mart)

````
The plot contains 4 tracks:
 - The chromosome, with a red marker on the gene location
 - The ensembl transcripts
 - The variant consequences, grouped by type (eg: INTRONIC, STOP LOST etc...), each green bar represent a variant
 - The cadd phred score, each dot represent a variant. 

![vargen visualisation](./images/SIM1_ENSG00000112246_GVIZ.png?raw=true)

The **rsid_highlight** parameter allows you to highlight some variants (by rsid) in red:
````
obesity_vargen_ann <- read.table("vargen_obesity/vargen_variants_annotated.tsv")

gene_mart <- connect_to_gene_ensembl()
vargen_visualisation(annotated_snps = obesity_vargen_ann, verbose = T,
                     outdir = "./obesity_gviz/", device = "png", 
                     gene_mart = gene_mart,
                     rsid_highlight = unique(obesity_vargen_ann[obesity_vargen_ann$cadd_phred > 20, "rsid"]))

````
![vargen visualisation with highlighted variants](./images/SIM1_ENSG00000112246_GVIZ_highlighted.png?raw=true)

### How to filter VarGen output

The filtering strategy is dependent on your study. VarGen outputs an R data.frame containing information about each variant. 
You can focus on clinically significant variants ("pathogenic", "likely pathogenic" etc...), on variants with a high phred 
score etc...

You can even combine different filtering, for example we found that keeping all the variants from the gwas catalog and with 
clinical significance while removing the variants with a cadd phred score lower than 10 was removing a lot of potentially
false positive results:

````
# If you have not installed the files already:
# vargen_install(install_dir = "./vargen_data/", verbose = T)

adipose_tissues <- select_gtex_tissues(gtex_dir = "./vargen_data/GTEx_Analysis_v8_eQTL/",
                                        tissues_query = "adipose")

obesity_vargen <- vargen_pipeline(vargen_dir = "./vargen_data/", omim_morbid = "601665",
                                  fantom_corr = 0.20, outdir = "./vargen_obesity/",
                                  gtex_tissues = adipose_tissues,
                                  gwascat_file = "./gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv",
                                  gwas_traits = "Obesity", verbose = T)

obesity_vargen_ann <- annotate_variants(obesity_vargen$rsid, verbose = T)
obesity_vargen_ann <- merge(obesity_vargen, obesity_vargen_ann)

# Filtering the vargen variants:
    vargen_phred_10 <- obesity_vargen_ann[obesity_vargen_ann$cadd_phred > 10,]
    vargen_phred_10 <- vargen_phred_10[!is.na(vargen_phred_10$cadd_phred),]
    vargen_clinVar <- obesity_vargen_ann[obesity_vargen_ann$clinical_significance != "",]
    vargen_gwas <- obesity_vargen_ann[obesity_vargen_ann$source == "gwas",]

    filtered_vargen <- rbind(vargen_phred_10, vargen_clinVar, vargen_gwas)
````
