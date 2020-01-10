## ----options, include=FALSE---------------------------------------------------
library("errRt")
basedir <- system.file("extdata", package="errRt")
knitr::opts_knit$set(width=120,
                     root.dir=basedir,
                     progress=TRUE,
                     verbose=TRUE,
                     echo=TRUE)
knitr::opts_chunk$set(error=TRUE,
                      dpi=96)
old_options <- options(digits=4,
                       stringsAsFactors=FALSE,
                       knitr.duplicate.label="allow")
ggplot2::theme_set(ggplot2::theme_bw(base_size=10))
rundate <- format(Sys.Date(), format="%Y%m%d")
##tmp <- sm(loadme(filename=paste0(gsub(pattern="\\.Rmd", replace="", x=previous_file), "-v", ver, ".rda.xz")))
##rmd_file <- "03_expression_infection_20180822.Rmd"

## ----installation, eval=FALSE-------------------------------------------------
#  devtools::install_github("abelew/errRt")

## ----getting_started----------------------------------------------------------
library(errRt)
basedir <- system.file("extdata", package="errRt")
metadata_file <- file.path(basedir, "all_samples.xlsx")
metadata_file
setwd(basedir)
meta <- hpgltools::extract_metadata(metadata_file)

## ----metadata, results='asis'-------------------------------------------------
knitr::kable(meta)

## ----create_matrices----------------------------------------------------------
error_data <- create_matrices(sample_sheet=metadata_file, ident_column="identtable",
                              mut_column="mutationtable", min_reads=3, min_indexes=3,
                              min_sequencer=10, min_position=24, prune_n=TRUE,
                              verbose=TRUE)

## ----mismatches_by_position---------------------------------------------------
gplots::heatmap.2(error_data[["normalized_by_counts"]][["miss_indexes_by_position"]],
                  trace="none",
                  Rowv=FALSE)
gplots::heatmap.2(error_data[["normalized_by_counts"]][["miss_sequencer_by_position"]],
                  trace="none",
                  Rowv=FALSE)

## ----reference_nucleotides, results='asis'------------------------------------
knitr::kable(error_data[["normalized_by_counts"]][["miss_indexes_by_ref_nt"]])
knitr::kable(error_data[["normalized_by_counts"]][["miss_sequencer_by_ref_nt"]])

## ----mutated_nucleotides, results='asis'--------------------------------------
knitr::kable(error_data[["normalized_by_counts"]][["miss_indexes_by_hit_nt"]])
knitr::kable(error_data[["normalized_by_counts"]][["miss_sequencer_by_hit_nt"]])

## ----type_nucleotides, results='asis'-----------------------------------------
knitr::kable(error_data[["normalized_by_counts"]][["miss_indexes_by_type"]])
knitr::kable(error_data[["normalized_by_counts"]][["miss_sequencer_by_type"]])

## ----trans_nucleotides, results='asis'----------------------------------------
knitr::kable(error_data[["normalized_by_counts"]][["miss_indexes_by_trans"]])
knitr::kable(error_data[["normalized_by_counts"]][["miss_sequencer_by_trans"]])

## ----strength_nucleotides, results='asis'-------------------------------------
knitr::kable(error_data[["normalized_by_counts"]][["miss_indexes_by_strength"]])
knitr::kable(error_data[["normalized_by_counts"]][["miss_sequencer_by_strength"]])

## ----insert_position, results='asis'------------------------------------------
knitr::kable(error_data[["normalized_by_counts"]][["insert_indexes_by_position"]])
knitr::kable(error_data[["normalized_by_counts"]][["insert_sequencer_by_position"]])

## ----insert_nt, results='asis'------------------------------------------------
knitr::kable(error_data[["normalized_by_counts"]][["insert_indexes_by_nt"]])
knitr::kable(error_data[["normalized_by_counts"]][["insert_sequencer_by_nt"]])

