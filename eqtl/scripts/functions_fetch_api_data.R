#-----------------------------------------------------------------
# Function: fetch data from GWAS catalog API or eQTL catalog API
# Author: Silva Kasela
#-----------------------------------------------------------------

# Modified from: http://htmlpreview.github.io/?https://github.com/kauralasoo/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.html

library(jsonlite)
library(dplyr)

# Function to easily fetch data from both GWAS and eQTL APIs ----

fetch_from_catalog_API <- function(link, is_gwas = FALSE){
  nullToNA <- function(x) {
    x[sapply(x, is.null)] <- NA
    return(x)
  }

  if (is_gwas) {
    cols_to_nest <- c("variant_id", "chromosome", "base_pair_location", "trait",
                      "p_value", "ci_lower", "ci_upper", "beta",
                      "effect_allele", "other_allele", "effect_allele_frequency",
                      "odds_ratio", "study_accession", "code")
  } else {
    cols_to_nest <- c("study_id", "qtl_group", "rsid",
                      "chromosome", "position", "pvalue", "condition_label",
                      "tissue_label", "molecular_trait_id", "gene_id", "ac",
                      "ref", "beta",  "variant", "an", "median_tpm",  "condition",
                      "r2", "alt", "type", "maf",  "tissue")
  }

  is_paginated <- !stringr::str_detect(link, "paginate=False")
  # message("isPagined:", is_paginated)
  page <- 1
  merged_summaries <- data.frame()
  while (!is.null(link)) {
    # print(paste0("Fetching page #",page))
    api_raw_data <- fromJSON(link, simplifyDataFrame = TRUE, flatten = TRUE)
    link <- api_raw_data$`_links`$`next`$href
    if (purrr::is_empty(api_raw_data$`_embedded`$associations)) {
      return(merged_summaries)
    }
    raw_list_data <- do.call(rbind, lapply(api_raw_data$`_embedded`$associations, rbind))
    raw_data <- nullToNA(raw_list_data) %>% as.matrix() %>% dplyr::as_tibble()
    if (is_paginated & is_gwas) { raw_data <- dplyr::select(raw_data, -c("_links")) }
    raw_data <- tidyr::unnest(raw_data, cols = cols_to_nest)
    if (!is.null(link)) {
      page <- page + 1
    }
    merged_summaries <- merged_summaries %>% rbind(raw_data) %>% as.data.frame()
  }
  return(merged_summaries[cols_to_nest])
}
