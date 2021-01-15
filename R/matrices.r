#' Given a samples sheet with some metadata, create a big pile of matrices
#' describing the data.
#'
#' @param sample_sheet xlsx/csv/whatever file of metadata.
#' @param ident_column Column containing the files of reads identical to the
#'  template.
#' @param mut_column Column containing the files of reads not identical to the
#'  template.
#' @param min_reads Filter for the minimum number of reads / index.
#' @param min_indexes Filter for the minimum numer of indexes / mutation.
#' @param min_sequencer Filter defining the minimum number of reads when looking
#'  for sequencer-based mutations.
#' @param min_position Filter against weird data at the beginning of the
#'  template.
#' @param max_position Filter against weird data at the end.
#' @param max_mutations_per_read Drop reads with more than this number of mutations.
#' @param prune_n Remove mutations which are something to N.
#' @param excel Save data to this file.
#' @param verbose Print information about the running of this function?
#' @param plot_order When plotting, force the order to this.
#' @param savefile Save the timeconsuming portion to this file.
#' @param overwrite If the savefile exists, overwrite it?
#' @return List with lots of matrices of the resulting data.
#' @export
create_matrices <- function(sample_sheet="sample_sheets/all_samples.xlsx",
                            ident_column="identtable",
                            mut_column="mutationtable",
                            min_reads=NULL, min_indexes=NULL,
                            min_sequencer=10, min_position=NULL,
                            max_position=NULL, max_mutations_per_read=NULL,
                            prune_n=TRUE, excel=NULL, verbose=TRUE,
                            plot_order=c("dna_control", "dna_low", "dna_high",
                                         "rna_control", "rna_low", "rna_high"),
                            savefile="sample_tables.rda", overwrite=FALSE) {
  meta <- hpgltools::extract_metadata(metadata=sample_sheet)
  if (is.null(meta[["sampletype"]])) {
    meta[["sampletype"]] <- paste0(meta[["type"]], "_", meta[["condition"]])
  }
  meta[["sampletype"]] <- factor(meta[["sampletype"]], levels=plot_order)
  meta <- with(meta, meta[order(sampletype, sampleid), ])
  meta[["sampleid"]] <- factor(meta[["sampleid"]], levels=meta[["sampleid"]])
  samples <- list()
  filtered_matrix <- NULL
  reads_remaining <- NULL
  indexes_remaining <- NULL
  sample_names <- meta[["sampleid"]]

  ## 1. Mutant indexes before filtering.
  pre_chng_index_density_dt <- data.table()
  ## 2. Mutant indexes after filtering.
  post_chng_index_density_dt <- data.table()
  ## 3. Table identical indexes before filtering.
  pre_ident_index_density_dt <- data.table()
  ## 4. identical indexes after filtering.
  post_ident_index_density_dt <- data.table()
  ## 5. Table of both, before filtering.
  pre_index_density_dt <- data.table()
  ## 6. Table of both, after filtering.
  post_index_density_dt <- data.table()

  for (s in 1:nrow(meta)) {
    identical <- meta[s, ident_column]
    changed <- meta[s, mut_column]
    sname <- as.character(meta[s, "sampleid"])
    if (verbose) {
      message("Starting sample: ", sname, ".")
    }
    samples[[s]] <- quantify_parsed(changed=changed, identical=identical, min_reads=min_reads,
                                    min_indexes=min_indexes, min_sequencer=min_sequencer,
                                    min_position=min_position, max_position=max_position,
                                    max_mutations_per_read=max_mutations_per_read,
                                    prune_n=prune_n, verbose=verbose)
    filtered_matrix <- cbind(filtered_matrix, samples[[s]][["reads_filtered"]])
    reads_remaining <- cbind(reads_remaining, samples[[s]][["reads_remaining"]])
    indexes_remaining <- cbind(indexes_remaining, samples[[s]][["indexes_remaining"]])

    ## Get index density before filters, mutant indexes
    ## suck it r cmd check
    index <- NULL
    index_chng_density <- data.table(
      "index" = samples[[s]][["starting_chng_indexes"]],
      "sample" = sname)
    index_chng_density <- index_chng_density %>%
      group_by(index, sample) %>%

  summarise(index_count = n(), .groups = "keep")
    ## Get index density after filters, mutant indexes
    changed_dt <- data.table(
      "index" = samples[[s]][["changed_table"]][["index"]],
      "sample" = sname)
    changed_dt <- changed_dt %>%
      group_by(index, sample) %>%
      summarise(index_count = n(), .groups = "drop")


    index_chng_density[["index"]] <- NULL
    ## Table 1 above, e.g. a table with rows:reads/index, columns:samples, cells:counts.
    pre_chng_index_density_dt <- make_index_table(pre_chng_index_density_dt,
                                                  index_chng_density, sname)


    changed_dt[["index"]] <- NULL
    ## Table 2 above.
    post_chng_index_density_dt <- make_index_table(post_chng_index_density_dt,
                                                   changed_dt, sname)

    ## Get index density before filters, indentical indexes
    index_ident_density <- data.table(
      "index" = samples[[s]][["starting_ident_indexes"]],
      "sample" = sname)
    index_ident_density <- index_ident_density %>%
      group_by(index, sample) %>%
      summarise(index_count = n(), .groups = "drop")
    index_ident_density[["index"]] <- NULL
    ## Table 3 above.
    pre_ident_index_density_dt <- make_index_table(pre_ident_index_density_dt,
                                                   index_ident_density, sname)

    ## Get index density after filter, identical indexes
    ident_indexes <- data.table(
      "index" = samples[[s]][["ident_table"]][["index"]],
      "sample" = sname)
    ident_indexes <- ident_indexes %>%
      group_by(index, sample) %>%
      summarise(index_count = n(), .groups = "drop")
    ident_indexes[["index"]] <- NULL
    ## Table 4 above.
    post_ident_index_density_dt <- make_index_table(post_ident_index_density_dt,
                                                    ident_indexes, sname)

    ## Table 5 above, the sum of pre_ident and pre_chng density tables.
    pre_index_density_dt <- add_index_densities(pre_index_density_dt,
                                                pre_ident_index_density_dt,
                                                pre_chng_index_density_dt, sname)

    ## Table 6 above, sum of post densities
    post_index_density_dt <- add_index_densities(post_index_density_dt,
                                                 post_ident_index_density_dt,
                                                 post_chng_index_density_dt, sname)
  } ## End iterating over every sample.

  names(samples) <- meta[["sampleid"]]
  colnames(filtered_matrix) <- meta[["sampleid"]]
  colnames(reads_remaining) <- meta[["sampleid"]]
  colnames(indexes_remaining) <- meta[["sampleid"]]

  ## I want to make it possible to reprocess raw data via savefiles,
  ## but at this time I am not sure where the best place is to add that logic.
  ## I think it must be earlier than this, and indeed would require the for()
  ## loop above to be split so that this logic could be added immediately after
  ## reading the raw data.
  write_rda <- TRUE
  if (is.null(savefile)) {
    write_rda <- FALSE
  }
  if (isFALSE(overwrite)) {
    write_rda <- FALSE
  }

  ## Remove summary rows which were not used
  used_rows <- rowSums(filtered_matrix) > 0
  filtered_matrix <- filtered_matrix[used_rows, ]
  used_rows <- rowSums(reads_remaining) > 0
  reads_remaining <- reads_remaining[used_rows, ]
  used_rows <- rowSums(indexes_remaining) > 0
  indexes_remaining <- indexes_remaining[used_rows, ]

  total_ident_reads <- reads_remaining["ident_index_pruned", ]
  total_mut_reads <- reads_remaining["mut_index_pruned", ]
  total_reads <- total_ident_reads + total_mut_reads
  reads_remaining <- rbind(reads_remaining, total_reads)

  tables <- c("miss_reads_by_position", "miss_indexes_by_position", "miss_sequencer_by_position",
              "miss_reads_by_string", "miss_indexes_by_string", "miss_sequencer_by_string",
              "miss_reads_by_refnt", "miss_indexes_by_refnt", "miss_sequencer_by_refnt",
              "miss_reads_by_hitnt", "miss_indexes_by_hitnt", "miss_sequencer_by_hitnt",
              "miss_reads_by_type", "miss_indexes_by_type", "miss_sequencer_by_type",
              "miss_reads_by_trans", "miss_indexes_by_trans", "miss_sequencer_by_trans",
              "miss_reads_by_strength", "miss_indexes_by_strength", "miss_sequencer_by_strength",
              "insert_reads_by_position", "insert_indexes_by_position", "insert_sequencer_by_position",
              "insert_reads_by_nt", "insert_indexes_by_nt", "insert_sequencer_by_nt",
              "delete_reads_by_position", "delete_indexes_by_position", "delete_sequencer_by_position",
              "delete_reads_by_nt", "delete_indexes_by_nt", "delete_sequencer_by_nt")

  ## Now build up the matrices of data where each matrix should have rownames which make
  ## sense for its name and the column names should be by sample.  Then the cells should
  ## be filled in with the data for each sample.
  pre_normalization_data <- matrices_from_tables(tables, samples, verbose=verbose)
  matrices <- pre_normalization_data[["matrices"]]
  indexes_per_sample <- pre_normalization_data[["indexes_per_sample"]]
  reads_per_sample <- pre_normalization_data[["reads_per_sample"]]

  ## Perform some normalizations of the data.
  normalized <- normalize_matrices(matrices, reads_per_sample, indexes_per_sample)

  retlist <- list(
    "metadata" = meta,
    "samples" = samples,
    "filtered" = filtered_matrix,
    "reads_remaining" = reads_remaining,
    "indexes_remaining" = indexes_remaining,
    "reads_per_sample" = reads_per_sample,
    "indexes_per_sample" = indexes_per_sample,
    "matrices" = matrices,
    "matrices_cpm" = normalized[["mat_cpm"]],
    "matrices_cpmlength" = normalized[["mat_cpmlength"]],
    "matrices_counts" = normalized[["mat_over_counts"]],
    "matrices_countslength" = normalized[["mat_over_countslength"]],
    "pre_chng_index_density_dt" = pre_chng_index_density_dt,
    "post_chng_index_density_dt" = post_chng_index_density_dt,
    "pre_ident_index_density_dt" = pre_ident_index_density_dt,
    "post_ident_index_density_dt" = post_ident_index_density_dt,
    "pre_index_density_dt" = pre_index_density_dt,
    "post_index_density_dt" = post_index_density_dt
    )
  retlist[["plots"]] <- barplot_matrices(retlist)
  if (!is.null(excel)) {
    retlist[["excel"]] <- write_matrices(retlist, excel, verbose=verbose)
  }
  return(retlist)
}

#' Given the many tables provided by quantify, merge them into matrices of error
#' types.
#'
#' @param tables The big pile of error tables.
#' @param samples The set of samples in the data.
#' @param verbose Print some information while running.
#' @return Matrices with error types as rows, samples as columns.
#' @export
matrices_from_tables <- function(tables, samples, verbose=FALSE) {
  matrices <- list()
  indexes_per_sample <- c()
  reads_per_sample <- c()
  num_tables <- length(tables)
  for (t in 1:num_tables) {
    table_name <- tables[t]
    col_names <- c()
    num_samples <- length(samples)
    for (s in 1:num_samples) {
      sample_name <- names(samples)[s]
      col_names <- c(col_names, sample_name)
      table <- samples[[sample_name]][[table_name]]
      last_col <- ncol(table)
      column <- table[[last_col]]
      if (is.null(matrices[[table_name]])) {
        matrices[[table_name]] <- data.table::data.table("names"=table[[1]], "s"=column)
        colnames(matrices[[table_name]]) <- c("names", sample_name)
      } else {
        tmpdt <- data.table::data.table(keep.rownames=TRUE, "names"=table[[1]], "s"=column)
        matrices[[table_name]] <- merge(matrices[[table_name]], tmpdt, by="names", all=TRUE)
        colnames(matrices[[table_name]])[s + 1] <- sample_name
      }
      indexes_per_sample[sample_name] <- samples[[sample_name]][["filtered_index_count"]]
      reads_per_sample[sample_name] <- sum(samples[[sample_name]][["changed_table"]][["all_reads"]]) +
        sum(samples[[sample_name]][["ident_table"]][["all_reads"]])
    } ## End iterating over every sample name.
    if (nrow(matrices[[table_name]]) > 0) {
      na_idx <- is.na(matrices[[t]])
      matrices[[t]][na_idx] <- 0
      numeric_rownames <- stringr::str_detect(string=rownames(matrices[[t]]),
                                              pattern="[[:digit:]]")
      if (sum(numeric_rownames) == length(numeric_rownames)) {
        numeric_order <- suppressWarnings(order(as.numeric(rownames(matrices[[t]]))))
        matrices[[t]] <- matrices[[t]][numeric_order, ]
      }
    }
  }
  retlist <- list(
    "matrices" = matrices,
    "indexes_per_sample" = indexes_per_sample,
    "reads_per_sample" = reads_per_sample)
  return(retlist)
}

#' Invoke some normalization methods on the various matrices from
#' matrices_from_tables().
#'
#' @param matrices Tables provided by matrices_from_tables().
#' @param reads_per_sample Numeric vector defining how many reads survived
#'  filtering for each sample.
#' @param indexes_per_sample Numeric vector defining how many indexes survived
#'  filtering for each sample.
#' @return Normalized matrices
#' @export
normalize_matrices <- function(matrices, reads_per_sample, indexes_per_sample) {
  ## I still need to figure out how I want to normalize these numbers.
  mat_cpm <- matrices
  mat_cpmlength <- matrices
  mat_over_counts <- matrices
  mat_over_countslength <- matrices
  normalized_by_counts <- matrices
  matrices_by_counts <- matrices
  sequence_length <- nrow(matrices[["miss_indexes_by_position"]])
  mnames <- names(matrices)

  for (t in 1:length(mnames)) {
    table_name <- mnames[t]
    if (nrow(matrices[[table_name]]) == 0) {
      mat_cpm[[table_name]] <- NULL
      mat_cpmlength[[table_name]] <- NULL
      mat_over_counts[[table_name]] <- NULL
      mat_over_countslength[[table_name]] <- NULL
      normalized_by_counts[[table_name]] <- NULL
      matrices_by_counts[[table_name]] <- NULL
      next
    }

    ## str_detect checks the table name for whether it is reads vs. indexes vs. sequencer.
    read_table <- stringr::str_detect(string=table_name, pattern="reads")

    ## Recast the data table as a matrix for doing cpm and such, also pull the names
    ## back out to their original purpose, rownames.
    data_columns <- 2:ncol(mat_cpm[[table_name]])
    mtrx <- as.matrix(mat_cpm[[table_name]][, 2:ncol(mat_cpm[[table_name]])])
    rownames(mtrx) <- mat_cpm[[table_name]][["names"]]

    ## For a few tables, the control sample is all 0s, so let us take that into account.
    usable_columns <- colSums(mtrx) > 0
    mtrx <- mtrx[, usable_columns]
    ## Now copy in place the new data
    mat_cpm[[table_name]] <- mtrx
    mat_cpmlength[[table_name]] <- mtrx
    mat_over_counts[[table_name]] <- mtrx
    mat_over_countslength[[table_name]] <- mtrx
    normalized_by_counts[[table_name]] <- mtrx
    matrices_by_counts[[table_name]] <- mtrx

    ## Now perform some normalizations
    if (isTRUE(read_table)) {
      mat_cpm[[table_name]] <- try(edgeR::cpm(y=mat_cpm[[table_name]],
                                              lib.sizes=reads_per_sample), silent=TRUE)
      used_columns <- colnames(mat_over_counts[[table_name]])
      mat_over_counts[[table_name]] <- mat_over_counts[[table_name]] / reads_per_sample[used_columns]
    } else {
      mat_cpm[[table_name]] <- try(edgeR::cpm(y=mat_cpm[[table_name]],
                                              lib.sizes=indexes_per_sample), silent=TRUE)
      used_columns <- colnames(mat_over_counts[[table_name]])
      mat_over_counts[[table_name]] <- mat_over_counts[[table_name]] / indexes_per_sample[used_columns]
    }
    mat_over_countslength[[table_name]] <- mat_over_counts[[table_name]] / sequence_length
    if (class(mat_cpm[[table_name]])[1] == "try-error") {
      message("Normalization failed for table: ", table_name, ", setting it to null.")
      mat_cpm[[table_name]] <- NULL
      mat_cpmlength[[table_name]] <- NULL
    } else {
      mat_cpmlength[[table_name]] <- mat_cpmlength[[table_name]] / sequence_length
    }
  }

  retlist <- list(
    "mat_cpm" = mat_cpm,
    "mat_cpmlength" = mat_cpmlength,
    "mat_over_counts" = mat_over_counts,
    "mat_over_countslength" = mat_over_countslength)
  return(retlist)
}


add_density_tables <- function(first, second) {
  first_df <- as.data.frame(first)
  colnames(first_df) <- c("num", "first")
  second_df <- as.data.frame(second)
  colnames(second_df) <- c("num", "second")
  both_df <- merge(first_df, second_df, by="num", all=TRUE)
  rownames(both_df) <- both_df[["num"]]
  both_df[["num"]] <- NULL
  sum <- rowSums(both_df, na.rm=TRUE)
  idx <- order(as.numeric(names(sum)), decreasing=FALSE)
  sum <- sum[idx]
  return(sum)
}

make_frequency_matrix <- function(mtrx) {
  new <- mtrx
  csums <- colSums(mtrx)
  for (col in colnames(mtrx)) {
    new[, col] <- new[, col] / csums[col]
  }
  multiplier <- 1 / max(new)
  new <- new * multiplier
  return(new)
}

#' Given the large table of index counts for one sample, bring them together into a small matrix.
#'
#' @param composite The matrix to build.
#' @param piece The piece to add.
#' @param sname Sample name.
make_index_table <- function(composite, piece, sname) {
  piece_table <- table(piece[["index_count"]])
  composite <- merge_index_table(composite, piece_table, sname)
  return(composite)
}

merge_index_table <- function(composite, piece_table, sname) {
  sname <- as.character(sname)
  piece_dt <- data.table("num" = names(piece_table), sname=as.numeric(piece_table))
  colnames(piece_dt) <- c("num", sname)
  if (ncol(composite) == 0) {
    composite <- piece_dt
  } else {
    composite <- merge(composite, piece_dt, by="num", all=TRUE)
  }
  return(composite)
}

#' Extract the numbers from the growing density matrix for this specific sample.
#'
#' Then add them and make a new matrix with these numbers.
#'
#' @param composite Table to build up.
#' @param first First piece to add to the composite.
#' @param second Second piece to add to the composite.
#' @param sname Sample name.
add_index_densities <- function(composite, first, second, sname) {
  sname <- as.character(sname)
  first_piece <- first[, c("num", sname), with=FALSE]
  colnames(first_piece) <- c("num", "first")
  second_piece <- second[, c("num", sname), with=FALSE]
  colnames(second_piece) <- c("num", "second")
  both <- merge(first_piece, second_piece, by="num", all=TRUE)
  num_names <- as.numeric(both[["num"]])
  both[["num"]] <- NULL
  sum_table <- rowSums(both, na.rm=TRUE)
  names(sum_table) <- num_names
  composite <- merge_index_table(composite, sum_table, sname)
  return(composite)
}
