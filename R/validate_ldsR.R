# ldsR compatibility gate for an exported directory.
# Usage: Rscript validate_ldsR.R <ldsR-dir> [bundled | <mask.parquet>]
# "bundled" (default) checks against ldsR's own extdata/common_snps.parquet, which is what
# partition_h2 uses; a mask path is only for test bundles that are not ldsR-compatible.
# Prints one JSON line: {"status": "pass"|"fail", ...}. Exit status 0 when the check ran.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ldsR)
})

args <- commandArgs(trailingOnly = TRUE)
dir <- args[[1]]
mask_arg <- if (length(args) >= 2) args[[2]] else "bundled"

result <- tryCatch({
  d <- ldsR::parse_parquet_dir(dir, read_ref = TRUE)
  mask_path <- if (mask_arg == "bundled") system.file("extdata/common_snps.parquet", package = "ldsR") else mask_arg
  mask <- arrow::read_parquet(mask_path, col_select = "common")
  annots <- d$annot$annot

  # partition_h2 indexes annot_ref rows with the bundled mask, so lengths must agree exactly
  stopifnot(nrow(d$annot_ref) == nrow(mask))
  m50_from_mask <- annots |>
    map_dbl(\(a) sum(d$annot_ref[[a]][mask$common]))
  m_from_ref <- annots |>
    map_dbl(\(a) sum(d$annot_ref[[a]]))

  checks <- tibble(
    annot = annots,
    m = d$annot$m,
    m50 = d$annot$m50,
    m_from_ref = m_from_ref,
    m50_from_mask = m50_from_mask
  ) |>
    mutate(ok = m == m_from_ref & m50 == m50_from_mask & m50 <= m)

  list(
    status = if (all(checks$ok) && !anyDuplicated(d$ld$SNP) && all(colnames(d$ld)[-1] == annots)) "pass" else "fail",
    ldsR_version = as.character(packageVersion("ldsR")),
    mask = mask_arg,
    n_ld_rows = nrow(d$ld),
    n_annot_ref_rows = nrow(d$annot_ref),
    checks = checks
  )
}, error = \(e) list(status = "fail", error = conditionMessage(e)))

cat(jsonlite::toJSON(result, auto_unbox = TRUE, dataframe = "rows"), "\n", sep = "")
