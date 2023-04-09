###############
# BindCompare: geneont.R
# Author: Pranav Mahableshwarkar
# Date: 04/12/2023
# Note: Performs GO analysis using gprofiler2.
# Output: GeneOntology HTML, Table Results
###############

library(optparse)
library(gprofiler2)

option_list <- list(
  make_option(c("-o", "--outdir"),
    type = "character", default = NULL,
    help = "directory for the gene ontology results to be saved"
  ),
  make_option(c("-c", "--expcondition"),
    type = "character", default = NULL,
    help = "the name of the experimental condition being examined"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

run_go <- function(output_dir, condition_string) {
  if (substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") {
    output_dir <- paste(output_dir, "/", sep = "")
  }
  significant <- read.csv("/tmp/gene_list.csv")
  print(significant)
  go_one_all <- as.character(significant[, c('Gene.ID')])
  print(go_one_all)

  gostres_one <- gost(
    query = go_one_all, organism = "dmelanogaster",
    ordered_query = FALSE, multi_query = FALSE,
    significant = TRUE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.05, correction_method = "g_SCS",
    domain_scope = "annotated", custom_bg = NULL,
    numeric_ns = "", sources = NULL, as_short_link = FALSE
  )
  p <- gostplot(gostres_one, capped = FALSE, interactive = FALSE)

  outpath <- paste(output_dir, condition_string, "_GeneOnt.png", sep = "")
  publish_gostplot(p,
    highlight_terms = gostres_one$result[c(1:3), ],
    width = NA, height = NA, filename = outpath
  )

  cols <- c("source", "term_name", "term_size", "intersection_size")
  outpath <- paste(output_dir, condition_string, "_GeneOntTable.pdf", sep = "")
  publish_gosttable(gostres_one,
    highlight_terms = gostres_one$result[c(1:150), ],
    use_colors = TRUE, show_columns = cols, filename = outpath
  )
}

run_go(opt$outdir, opt$expcondition)