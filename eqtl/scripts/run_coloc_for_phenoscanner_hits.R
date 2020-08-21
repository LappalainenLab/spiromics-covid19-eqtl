#!/usr/bin/env Rscript
#---------------------------------------------------------
# Colocalization analysis for COVID-19 relevant eGenes
# Author: Silva Kasela
#---------------------------------------------------------

# Run coloc for the best trait based on GWAS P-value for each of the different `efo_category` from every `efo_category_parent`

# If `coloc_method = mask` then using masking in the GWAS dataset
##    LD from the corresponding 1000G population (for example, CEU if EUR was the discovery population in GWAS)

# Load libaries
library(here)
library(data.table)
library(coloc)

# Source functions
source(here("scripts", "functions_fetch_api_data.R"))
source(here("scripts", "functions_run_coloc.R"))
source(here("scripts", "functions_plot_locus.R"))
source(here("scripts", "functions_get_genomatrix.R"))
## Some little fixes (using the development version of condmask)
source(here("scripts", "functions_coloc_signals_adj.R"))

# Arguments from command line
args <- commandArgs(trailingOnly = TRUE)
cat(args, fill = T)
gene_name <- args[1]
coloc_method <- args[2] # standard or mask

QTL_PREFIX <- sub("_", " ", "cis_eQTL")
QTL_TRAIT <- "cis"
QTL_N <- 144
QTL_SDYFILE <- "coloc/input/qtl/spiromics_sdY.txt"
QTL_COLOC <- 'beta,single'
GWAS_COLOC <- ifelse(coloc_method %in% "standard", "pvalues,single",
                     ifelse(coloc_method %in% "mask", "pvalues,mask", NA))
COLOC_MODE <- "iterative"
PRIORS <- '1e-4,1e-4,5e-6'
LD_SHOW <- FALSE
LD_VARIANT <- NA
LOCUSCOMPARE <- TRUE
LOCUSZOOM <- TRUE
SENSITIVITY <- TRUE
OUTPATH <- ifelse(coloc_method %in% "standard", "coloc/result/eqtl/single",
                  ifelse(coloc_method %in% "mask", "coloc/result/eqtl/mask", NA))
CIS_WINDOW <- 250000 # 500-kb region centered on lead eQTL (+/-250 kb from the lead variant)

# Data ----------

# 1. PhenoScanner results
version <- "v11"
phenoscanner <- read.table(here("eqtl", "cis", "phenoscanner", paste0("phewas_phenoscanner.covid19_related_egenes_", version, ".txt")), header = T, sep = "\t", stringsAsFactors = F)
phenoscanner <- phenoscanner[phenoscanner$gene_name %in% gene_name,]
phenoscanner$chr <- sapply(phenoscanner$hg38_coordinates, function(x) sub("chr", "", unlist(strsplit(x, ":"))[1]))
phenoscanner$position <- as.numeric(sapply(phenoscanner$hg38_coordinates, function(x) unlist(strsplit(x, ":"))[2]))
addmargins(table(phenoscanner$efo_category_parent))
gene_id <- unique(phenoscanner$gene_id)
cat(gene_id, fill = T)

## Double-check - does this eGene have multiple eVariants?
indep_variants <- read.table(here("eqtl", "cis", "independent", "spiromics.cis_independent_qtl.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)
indep_variants <- indep_variants[indep_variants$phenotype_id %in% gene_id,]
if (nrow(indep_variants) != 1) stop("Provide conditional summary stats!")

# 2. UKBB v3 manifest, https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=178908679
ukbb_manifest <- read.table(here("coloc", "UKBB_GWAS_Imputed_v3_File_Manifest_Release_20180731_Manifest_201807.tsv"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")

# 3. Studies with summary statistics from GWAS catalog, https://www.ebi.ac.uk/gwas/downloads/summary-statistics
gwas_catalog <- read.table(here("coloc", "list_gwas_summary_statistics_13_Apr_2020.csv"), header = T, sep = ",", stringsAsFactors = F)

# 4. eQTL allpairs data
qtl <- fread(here("coloc", "input", "qtl", paste0("spiromics.cis_eqtl.covid19_related_egenes.allpairs.", unique(phenoscanner$chr), ".txt.gz")), header = T, sep = "\t", stringsAsFactors = F, data.table = F)
stopifnot(gene_id %in% qtl$phenotype_id)
qtl <- qtl[qtl$phenotype_id %in% gene_id, ]
colnames(qtl)[c(4, 7:9)] <- c("qtl_maf", "qtl_pval", "qtl_beta", "qtl_se")
qtl$id <- paste0(qtl$chr, "_", qtl$pos) #great id to link with gwas file
qtl$chr <- sub("chr", "", qtl$chr)
# Add study type and sample size to merged dataset
qtl$qtl_n <- QTL_N
# Exclude probes that fall into the MHC region on chr6
if (paste0("chr", unique(phenoscanner$chr)) == "chr6") {
  qtl$pos_pheno <- qtl$pos - qtl$tss_distance - 1
  qtl <- qtl[!(qtl$pos_pheno > 28510120 & qtl$pos_pheno < 33480577), ]
  message(paste0("No. of significant QTLs not falling into MHC region: ", length(unique(qtl$phenotype_id))))
  qtl$pos_pheno <- NULL
}
# Add sdY
if (!is.na(QTL_SDYFILE)) {
  qtl_sdy <- read.table(here(QTL_SDYFILE), header = F, sep = "\t", stringsAsFactors = F)
  stopifnot(gene_id %in% qtl_sdy$V1)
  qtl$qtl_sdY <- qtl_sdy[match(gene_id, qtl_sdy$V1), "V2"]
  cat("sdY is ", qtl_sdy[match(gene_id, qtl_sdy$V1), "V2"], fill = TRUE)
}
# Replace Ensembl IDs with HGNC IDs
qtl$phenotype_id <- gene_name

# 5. rsID lookup table
lookup <- fread(here("genotype", "freeze9.pass_only.spiromics_144samples.maf05.biallelic.lookup_table.txt.gz"), header = T, sep = "\t", stringsAsFactors = F, data.table = F)
stopifnot(qtl$variant_id %in% lookup$variant_id)
lookup <- lookup[lookup$variant_id %in% qtl$variant_id,]
qtl$qtl_rsid <- lookup[match(qtl$variant_id, lookup$variant_id), "rs_id_dbSNP151_GRCh38p7"]

# Process ----

# Replace traits with the exact UKBB trait names
if ("Forced vital capacity, best measure" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Forced vital capacity, best measure" & phenoscanner$pmid == "UKBB", "trait"] <- "Forced vital capacity (FVC), Best measure"
}
if ("Peak expiratory flow" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Peak expiratory flow" & phenoscanner$pmid == "UKBB", "trait"] <- "Peak expiratory flow (PEF)"
}
if ("Forced expiratory volume in 1-second, predicted" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Forced expiratory volume in 1-second, predicted" & phenoscanner$pmid == "UKBB", "trait"] <- "Forced expiratory volume in 1-second (FEV1), predicted"
}
if ("Height" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Height" & phenoscanner$pmid == "UKBB", "trait"] <- "Standing height"
}
if ("Body mass index" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Body mass index" & phenoscanner$pmid == "UKBB", "trait"] <- "Body mass index (BMI)"
}
if ("Arm fat mass right" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Arm fat mass right" & phenoscanner$pmid == "UKBB", "trait"] <- "Arm fat mass (right)"
}
if ("Heel bone mineral density" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Heel bone mineral density" & phenoscanner$pmid == "UKBB", "trait"] <- "Heel bone mineral density (BMD)"
}
if ("Diastolic blood pressure" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Diastolic blood pressure" & phenoscanner$pmid == "UKBB", "trait"] <- "Diastolic blood pressure, automated reading"
}
if ("Treatment with metformin" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Treatment with metformin" & phenoscanner$pmid == "UKBB", "trait"] <- "Treatment/medication code: metformin"
}
if ("Hayfever, allergic rhinitis or eczema" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Hayfever, allergic rhinitis or eczema" & phenoscanner$pmid == "UKBB", "trait"] <- "Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor: Hayfever, allergic rhinitis or eczema"
}
if ("Hand grip strength left" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Hand grip strength left" & phenoscanner$pmid == "UKBB", "trait"] <- "Hand grip strength (left)"
}
if ("Pain type experienced in last month: knee pain" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Pain type experienced in last month: knee pain" & phenoscanner$pmid == "UKBB", "trait"] <- "Pain type(s) experienced in last month: Knee pain"
}
if ("Pain type experienced in last month: none of the above" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Pain type experienced in last month: none of the above" & phenoscanner$pmid == "UKBB", "trait"] <- "Pain type(s) experienced in last month: None of the above"
}
if ("Qualifications: college or university degree" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Qualifications: college or university degree" & phenoscanner$pmid == "UKBB", "trait"] <- "Qualifications: College or University degree"
}

# GWAS catalog trait names
if ("Age-related macular degeneration" %in% phenoscanner$trait) {
  phenoscanner[phenoscanner$trait %in% "Age-related macular degeneration" & phenoscanner$pmid != "UKBB", "trait"] <- "Advanced age-related macular degeneration"
}

addmargins(table(phenoscanner$efo_category_parent, phenoscanner$efo_category))

# Setup coloc -------------------

# Specify priors
PRIORS <- unlist(strsplit(PRIORS, ","))
p1used <- as.numeric(PRIORS[1])
p2used <- as.numeric(PRIORS[2])
p12used <- as.numeric(PRIORS[3])
cat("Using p1: ", p1used, fill = T)
cat("Using p2: ", p2used, fill = T)
cat("Using p12: ", p12used, fill = T)

# Coloc mode
qtl_coloc_input <- unlist(strsplit(QTL_COLOC, ","))[1]
qtl_coloc_mode <- unlist(strsplit(QTL_COLOC, ","))[2]
stopifnot(!is.na(GWAS_COLOC))
gwas_coloc_input <- unlist(strsplit(GWAS_COLOC, ","))[1]
gwas_coloc_mode <- unlist(strsplit(GWAS_COLOC, ","))[2]

# Process -----------------------

j <- 1
coloc_result <- list()
cat("-----------------------------", fill = TRUE)
for (efo_category_parent in unique(phenoscanner$efo_category_parent)) {
  cat(efo_category_parent, fill = TRUE)
  sel <- phenoscanner[phenoscanner$efo_category_parent %in% efo_category_parent, ]
  for (efo_category in unique(sel$efo_category)) {
    cat(efo_category, fill = TRUE)
    sel_1 <- sel[sel$efo_category %in% efo_category, ]

    # Sort by GWAS p-values
    sel_1 <- sel_1[order(sel_1$p), ]

    # Find one trait that we can get summary stats for
    for (k in 1:nrow(sel_1)) {
      print(k)
      gwas_data <- NULL
      # Select region centered on lead eQTL
      bp_lower <- sel_1$position[k] - CIS_WINDOW
      bp_upper <- sel_1$position[k] + CIS_WINDOW

      # Get GWAS data
      if (sel_1$pmid[k] %in% "UKBB") {
        u <- ukbb_manifest[ukbb_manifest$Phenotype.Description %in% sel_1$trait[k], ]
        if (nrow(u) == 0) {
          cat("Warning: no match from UKBB manifest", fill = TRUE)
          next
        }
        if (sum(grepl("irnt", u$Phenotype.Code)) > 0) u <- u[grepl("irnt", u$Phenotype.Code),]
        if (sum(grepl("both_sexes", u$Sex)) > 0) u <- u[grepl("both_sexes", u$Sex),]
        # BMI has two codes (Anthropometry and Impedance measure), select one
        if (sum(grepl("21001", u$Phenotype.Code)) > 0) u <- u[grepl("21001", u$Phenotype.Code),]
        # Weight has two codes (Anthropometry and Impedance measure), select one
        if (sum(grepl("21002", u$Phenotype.Code)) > 0) u <- u[grepl("21002", u$Phenotype.Code),]
        if (nrow(u) != 1) {
          cat("Warning: Check match from UKBB manifest", fill = TRUE)
          next
        }

        # Summary stats from UKBB
        command <- unlist(strsplit(u$wget.command, " -O "))
        if (!command[2] %in% list.files(path = here("coloc", "input", "gwas"))) {
          # Get summary stats if not downloaded already
          system(paste0(command[1], " -O ", here("coloc/input/gwas/"), command[2]))
        }
        gwas_data_lookup <- fread(here("coloc", "input", "gwas", "ukbb_variants.grch38.txt.gz"), header = T, sep = "\t", stringsAsFactors = F, data.table = F)
        gwas_data_lookup <- gwas_data_lookup[gwas_data_lookup$GRCh38_chr %in% paste0("chr", sel_1$chr[k]), ]
        gwas_data_lookup <- gwas_data_lookup[gwas_data_lookup$GRCh38_pos > bp_lower & gwas_data_lookup$GRCh38_pos < bp_upper, ]

        gwas_data <- fread(cmd = paste0("gunzip -c ", here("coloc", "input", "gwas", command[2])), header = T, sep = "\t", stringsAsFactors = F, data.table = F)
        gwas_data <- gwas_data[gwas_data$variant %in% gwas_data_lookup$variant, ]
        cat("Downloaded ", nrow(gwas_data), " associations from UKBB (", command[2], ")", fill = TRUE)
        gwas_data_file <- command[2]

        # case-control or quant trait
        if ("expected_case_minor_AC" %in% colnames(gwas_data)) {
          cat("Case-control trait", fill = TRUE)
          pheno_split <- unlist(strsplit(command[2], "[.]"))
          ukbb_phenotypes <- read.table(here("coloc", "input", "gwas", paste0("phenotypes.", pheno_split[4], ".tsv.gz")), header = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")
          ukbb_phenotypes <- ukbb_phenotypes[ukbb_phenotypes$phenotype %in% pheno_split[1], ]
          stopifnot(nrow(ukbb_phenotypes) == 1)
          gwas_data$gwas_type <- "cc"
          gwas_data$gwas_s <- ukbb_phenotypes$n_cases/ukbb_phenotypes$n_non_missing
          colnames(gwas_data)[c(3, 6, 12)]  <- c("gwas_maf", "gwas_n", "gwas_pval")
        } else if ("expected_min_category_minor_AC" %in% colnames(gwas_data)) {
          colnames(gwas_data)[c(3, 6, 12)]  <- c("gwas_maf", "gwas_n", "gwas_pval")
          gwas_data$gwas_type <- "quant"
        } else {
          cat("Quantitative trait", fill = TRUE)
          colnames(gwas_data)[c(3, 5, 11)]  <- c("gwas_maf", "gwas_n", "gwas_pval")
          gwas_data$gwas_type <- "quant"
        }
        gwas_data$gwas_rsid <- gwas_data_lookup[match(gwas_data$variant, gwas_data_lookup$variant), "rsid"]
        gwas_data$id <- gwas_data_lookup[match(gwas_data$variant, gwas_data_lookup$variant), "GRCh38_id"]

        # GWAS_TRAIT for figures
        GWAS_TRAIT <- unlist(strsplit(sel_1$trait[k], " "))
        GWAS_TRAIT <- gsub("[^[:alnum:][:blank:]+?&/\\-]", "", GWAS_TRAIT)
        GWAS_TRAIT <- sub("/", "_", GWAS_TRAIT)
        GWAS_TRAIT <- paste(GWAS_TRAIT, collapse = "_")
        cat(GWAS_TRAIT, fill = TRUE)

      } else {
        g <- gwas_catalog[gwas_catalog$PubMed.ID %in% sel_1$pmid[k] &
                            (gwas_catalog$Reported.trait %in% sel_1$trait[k] | gwas_catalog$Trait.s. %in% tolower(sel_1$trait[k]) | gwas_catalog$Trait.s. %in% sel_1$efo_category[k]), ]
        if (nrow(g) > 1) {
          g <- gwas_catalog[gwas_catalog$PubMed.ID %in% sel_1$pmid[k] & gwas_catalog$Reported.trait %in% sel_1$trait[k],]
          if (nrow(g) > 1) stop("More than one match from GWAS catalog")
        }
        if (nrow(g) == 0) {
          cat("Warning: no match from GWAS catalog", fill = TRUE)
          next
        }
        if (!grepl("API", g$Data.access)) {
          cat("Warning: FTP download only for GWAS catalog", fill = TRUE)
          next
        }

        # Download summary stats using API
        gwas_query_str <- paste0("https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/", sel_1$chr[k],
                                 "/associations?study_accession=", g$Study.accession,
                                 "&bp_lower=", bp_lower, "&bp_upper=", bp_upper, "&size=1000")
        gwas_data <- fetch_from_catalog_API(link = gwas_query_str, is_gwas = TRUE)
        cat("Downloaded ", nrow(gwas_data), " associations from GWAS Catalogue (study accession code", g$Study.accession, ")", fill = TRUE)
        gwas_data_file <- g$Study.accession
        colnames(gwas_data)[c(1, 5)] <- c("gwas_rsid", "gwas_pval")
        gwas_data$id <- paste0("chr", gwas_data$chromosome, "_", gwas_data$base_pair_location)
        gwas_data$gwas_maf <- sapply(gwas_data$effect_allele_frequency, function(x){min(x, 1 - x)})
        dupl <- gwas_data$id[duplicated(gwas_data$id)]
        cat("Remove duplicated positions in GWAS file: ", length(dupl), fill = TRUE)
        gwas_data <- gwas_data[!gwas_data$id %in% dupl,]
        gwas_data$gwas_n <- sel_1$n[k]
        if (sel_1$n_cases[k] != 0) {
          gwas_data$gwas_type <- "cc"
          gwas_data$gwas_s <- as.numeric(sel_1$n_cases[k]) / as.numeric(sel_1$n[k])
        } else {
          gwas_data$gwas_type <- "quant"
        }
        if (gwas_coloc_input == "beta") {
          # Get se based on p-value and beta
          gwas_data$gwas_beta <- gwas_data$beta
          gwas_data$zscore <- ifelse(gwas_data$gwas_beta < 0,
                                     qnorm(p = gwas_data$gwas_pval/2, lower.tail = F) * (-1),
                                     qnorm(p = gwas_data$gwas_pval/2, lower.tail = F))
          gwas_data$gwas_se <- gwas_data$gwas_beta / gwas_data$zscore
        }
        GWAS_TRAIT <- paste(unlist(strsplit(tolower(g$Reported.trait), " ")), collapse = "_")
        cat(GWAS_TRAIT, fill = TRUE)
      }
      # End for loop, if we have found one example GWAS trait
      break
    }

    # if summary data not available
    null_results_standard <- data.frame("nsnps" = NA,	"PP.H0.abf" = NA,	"PP.H1.abf" = NA, "PP.H2.abf" = NA,	"PP.H3.abf" = NA,	"PP.H4.abf" = NA,
                                        "qtl" = NA,	"phenotype_id" = gene_name,	"qtl_lead_snv" = NA, "qtl_lead_snv_pval" = NA, "qtl_lead_snv_pos" = NA,
                                        "gwas" = sel_1$trait[1], "gwas_lead_snv" = sel_1$rsid[1],	"gwas_lead_snv_pval" = sel_1$p[1], "gwas_lead_snv_pos" = sel_1$position[1],
                                        "gwas_hit_nearest_log10p_5" = NA,	"chr" = NA,	"pos_phenotype" = NA,	"region_start" = NA, "region_end" = NA,
                                        "gwas_file" = NA, "efo_category" = sel_1$efo_category[1], "efo_category_parent" = sel_1$efo_category_parent[1],
                                        stringsAsFactors = F)
    null_results_mask <- data.frame("hit2" = NA, "hit1" = NA, "nsnps" = NA, "PP.H0.abf" = NA,	"PP.H1.abf" = NA, "PP.H2.abf" = NA,	"PP.H3.abf" = NA,	"PP.H4.abf" = NA,
                                    "best1" = NA, "best2" = NA, "best4" = NA, "hit1.margz" = NA, "hit2.margz" = NA,
                                    "qtl" = NA,	"phenotype_id" = gene_name,	"qtl_lead_snv" = NA, "qtl_lead_snv_pval" = NA, "qtl_lead_snv_pos" = NA,
                                    "gwas" = sel_1$trait[1], "gwas_lead_snv" = sel_1$rsid[1],	"gwas_lead_snv_pval" = sel_1$p[1], "gwas_lead_snv_pos" = sel_1$position[1],
                                    "gwas_hit_nearest_log10p_5" = NA,	"chr" = NA,	"pos_phenotype" = NA,	"region_start" = NA, "region_end" = NA,
                                    "gwas_file" = NA, "efo_category" = sel_1$efo_category[1], "efo_category_parent" = sel_1$efo_category_parent[1],
                                    stringsAsFactors = F)

    if (is.null(gwas_data)) {
      cat("No summary data for this trait", fill = TRUE)
      cat("-----------------------------", fill = TRUE)
      if (coloc_method %in% "standard") {
        coloc_result[[j]] <- null_results_standard
      } else {
        coloc_result[[j]] <- null_results_mask
      }
      j <- j + 1
      next
    }

    # Clean gwas data
    gwas_data <- gwas_data[!is.na(gwas_data$gwas_maf), ]
    gwas_data <- gwas_data[!is.na(gwas_data$gwas_pval), ]

    if (nrow(gwas_data) == 0) {
      cat("No summary data for this trait after filtering", fill = TRUE)
      cat("-----------------------------", fill = TRUE)
      if (coloc_method %in% "standard") {
        coloc_result[[j]] <- null_results_standard
      } else {
        coloc_result[[j]] <- null_results_mask
      }
      j <- j + 1
      next
    }

    if (sum(gwas_data$gwas_pval < 1e-5) == 0) {
      cat("No GWAS hits with P < 1e-5", fill = TRUE)
      cat("-----------------------------", fill = TRUE)
      if (coloc_method %in% "standard") {
        coloc_result[[j]] <- null_results_standard
      } else {
        coloc_result[[j]] <- null_results_mask
      }
      j <- j + 1
      next
    }

    # Merge eqtl and gwas data -----------------------------
    dat <- merge(qtl, gwas_data, by = "id", all = F, sort = F)
    message("No. of variants in the merged file: ", nrow(dat))
    ## in case there is no overlap between qtl and gwas, set all.x=T and have gwas pvalue =1 and frequency copied from qtl file
    if (nrow(dat) == 0) {
      stop("ERROR: no overlap between molQTL and GWAS")
      ## . GWAS p-values will be NA"
      #dat.mm <- merge(dat, gwas, by = "id", all.x = T, sort = F)
      #dat.mm$pvalue = 1
      #dat.mm$frequency = dat.mm$maf
    }

    # Order by phenotype_id
    ord <- order(dat$pos)
    dat <- dat[ord,]

    # Run coloc ------------------

    # Specify 1000G population for LD if needed
    if (coloc_method %in% "standard") {
      INDIV <- NULL
      GENOFILE <- NA
    } else {
      # There is one study that is using individuals from Japanese population (the IDE gene varaint and association with BMI, pmid 28892062),
      # other studies using individuals of European ancestry
      pop_name <- ifelse(sel_1$ancestry[k] %in% "European", "CEU",
                         ifelse(sel_1$ancestry[k] %in% "East Asian", "JPT", NA))
      stopifnot(!is.na(pop_name))
      pop_1kg <- read.table("~/data/1kg/gazal_et_al_2019.table_S4.filtered_unrelated_outbred.txt", header = T, sep = "\t", stringsAsFactors = F)
      cat("For LD using", pop_name, "population from 1000G", fill = T)
      INDIV <- pop_1kg[pop_1kg$POP %in% pop_name, "IID"]
      GENOFILE <- paste0("~/data/1kg/phase3_GRCh38/ALL.chr", sel_1$chr[k], ".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
    }

    coloc_result[[j]] <- run_coloc(phenotype = gene_name, data = dat,
                  qtl_coloc_input, qtl_coloc_mode,
                  gwas_coloc_input, gwas_coloc_mode,
                  genofile = GENOFILE, indiv = INDIV,
                  locuscompare = LOCUSCOMPARE,
                  locuscompare_fig = paste0(OUTPATH, "/fig_locuscompare/locuscompare.", QTL_TRAIT, "_", tolower(unlist(strsplit(QTL_PREFIX, " "))[2]), "_", gene_name,
                                            ".GWAS_", unlist(strsplit(sel_1$study[k], " "))[1], "_", GWAS_TRAIT),
                  locuszoom = LOCUSZOOM,
                  locuszoom_fig = paste0(OUTPATH, "/fig_locuszoom/locuszoom.", QTL_TRAIT, "_", tolower(unlist(strsplit(QTL_PREFIX, " "))[2]), "_", gene_name,
                                         ".GWAS_", unlist(strsplit(sel_1$study[k], " "))[1], "_", GWAS_TRAIT),
                  plot_main = paste0(gsub("_", " ", GWAS_TRAIT), " GWAS and ", QTL_PREFIX),
                  ld_show = LD_SHOW, ld_variant = LD_VARIANT,
                  p1 = p1used, p2 = p2used, p12 = p12used,
                  coloc_mode = COLOC_MODE,
                  sensitivity = SENSITIVITY,
                  sensitivity_fig = paste0(OUTPATH, "/fig_sensitivity/sensitivity.", QTL_TRAIT, "_", tolower(unlist(strsplit(QTL_PREFIX, " "))[2]), "_", gene_name,
                                           ".GWAS_", unlist(strsplit(sel_1$study[k], " "))[1], "_", GWAS_TRAIT),
                  qtl_trait = paste0(QTL_TRAIT, "_", tolower(unlist(strsplit(QTL_PREFIX, " "))[2]), "_bronch"),
                  gwas_trait = paste0(unlist(strsplit(sel_1$study[k], " "))[1], "_", GWAS_TRAIT))
    coloc_result[[j]]$gwas_file <- gwas_data_file
    coloc_result[[j]]$efo_category <- sel_1$efo_category[k]
    coloc_result[[j]]$efo_category_parent <- sel_1$efo_category_parent[k]
    j <- j + 1
    cat("-----------------------------", fill = TRUE)
  }
  cat("-----------------------------", fill = TRUE)
}

coloc_result <- do.call(rbind, coloc_result)

# Write out results -------------------
write.table(coloc_result, file = here(OUTPATH, paste0("coloc_results.", gene_name, ".txt")), col.names = T, row.names = F, sep = "\t", quote = F)

cat("Done!", fill = TRUE)
