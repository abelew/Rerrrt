#' Given a samples sheet with some metadata, create a big pile of matrices
#' describing the data.
#'
#' @param sample_sheet xlsx/csv/whatever file of metadata.
#' @param ident_column Column containing the files of reads identical to the
#'   template.
#' @param mut_column Column containing the files of reads not identical to the
#'   template.
#' @param min_reads Filter for the minimum number of reads / index.
#' @param min_indexes Filter for the minimum numer of indexes / mutation.
#' @param min_sequencer Filter defining the minimum number of reads when looking
#'   for sequencer-based mutations.
#' @param min_position Filter against weird data at the beginning of the
#'   template.
#' @param prune_n Remove mutations which are something to N.
#' @param verbose Print information about the running of this function?
#' @return List with lots of matrices of the resulting data.
#' @export
create_matrices <- function(sample_sheet="sample_sheets/all_samples.xlsx",
                            ident_column="identtable",
                            mut_column="mutationtable",
                            min_reads=NULL, min_indexes=NULL,
                            min_sequencer=10, min_position=NULL,
                            max_position=NULL, max_mutations_per_read=NULL,
                            prune_n=TRUE, excel=NULL, verbose=TRUE) {
  meta <- hpgltools::extract_metadata(metadata=sample_sheet)
  samples <- list()
  filtered_matrix <- NULL
  reads_remaining <- NULL
  indexes_remaining <- NULL
  pre_index_density_df <- NULL
  post_index_density_df <- NULL
  for (s in 1:nrow(meta)) {
    if (verbose) {
      message("Starting sample: ", s, ".")
    }
    identical <- meta[s, ident_column]
    changed <- meta[s, mut_column]
    sname <- meta[s, "sampleid"]
    samples[[s]] <- quantify_parsed(changed=changed, identical=identical, min_reads=min_reads,
                                    min_indexes=min_indexes, min_sequencer=min_sequencer,
                                    min_position=min_position, max_position=max_position,
                                    max_mutations_per_read=max_mutations_per_read,
                                    prune_n=prune_n, verbose=verbose)
    filtered_matrix <- cbind(filtered_matrix, samples[[s]][["reads_filtered"]])
    reads_remaining <- cbind(reads_remaining, samples[[s]][["reads_remaining"]])
    indexes_remaining <- cbind(indexes_remaining, samples[[s]][["indexes_remaining"]])
    index_chng_density <- as.data.frame(samples[[s]][["starting_chng_indexes"]])
    index_chng_density[["sample"]] <- paste0(sname, "_changed")
    colnames(index_chng_density) <- c("index", "sample")
    pre_index_density_df <- rbind(pre_index_density_df, index_chng_density)
    index_chng_density <- samples[[s]][["changed_table"]][, "index"]
    index_chng_density[["sample"]] <- paste0(sname, "_changed")
    colnames(index_chng_density) <- c("index", "sample")
    post_index_density_df <- rbind(post_index_density_df, index_chng_density)
    index_ident_density <- as.data.frame(samples[[s]][["starting_ident_indexes"]])
    index_ident_density[["sample"]] <- paste0(sname, "_identical")
    colnames(index_ident_density) <- c("index", "sample")
    pre_index_density_df <- rbind(pre_index_density_df, index_ident_density)
    index_chng_density <- samples[[s]][["ident_table"]][, "index"]
    index_chng_density[["sample"]] <- paste0(sname, "_identical")
    colnames(index_chng_density) <- c("index", "sample")
    post_index_density_df <- rbind(post_index_density_df, index_chng_density)
  }
  names(samples) <- meta[["sampleid"]]
  colnames(filtered_matrix) <- meta[["sampleid"]]
  colnames(reads_remaining) <- meta[["sampleid"]]
  colnames(indexes_remaining) <- meta[["sampleid"]]

  if (verbose) {
    message("Plotting index densities.")
  }
  pre_index_density_plot <- plot_index_density(pre_index_density_df)
  post_index_density_plot <- plot_index_density(post_index_density_df)

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
    "pre_index_density_plot" = pre_index_density_plot,
    "post_index_density_plot" = post_index_density_plot)
  retlist[["plots"]] <- barplot_matrices(retlist)
  if (!is.null(excel)) {
    retlist[["excel"]] <- write_matrices(retlist, excel)
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
  for (t in tables) {
    inc <- 0
    col_names <- c()
    for (s in names(samples)) {
      inc <- inc + 1
      col_names <- c(col_names, s)
      table <- samples[[s]][[t]]
      last_col <- ncol(table)
      column <- table[[last_col]]
      if (is.null(matrices[[t]])) {
        matrices[[t]] <- data.frame(row.names=table[[1]], s=column)
      } else {
        tmpdf <- data.frame(row.names=table[[1]], s=column)
        matrices[[t]] <- merge(matrices[[t]], tmpdf, by="row.names", all=TRUE)
        rownames(matrices[[t]]) <- matrices[[t]][["Row.names"]]
        matrices[[t]][["Row.names"]] <- NULL
      }
      indexes_per_sample[s] <- samples[[s]][["filtered_index_count"]]
      reads_per_sample[s] <- sum(samples[[s]][["changed_table"]][["all_reads"]]) +
        sum(samples[[s]][["ident_table"]][["all_reads"]])
    } ## End iterating over every sample name.
    if (nrow(matrices[[t]]) > 0) {
      colnames(matrices[[t]]) <- col_names
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
#'   filtering for each sample.
#' @param indexes_per_sample Numeric vector defining how many indexes survived
#'   filtering for each sample.
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
      message("Skipping table: ", table_name)
      next
    }
    read_table <- stringr::str_detect(string=table_name, pattern="reads")

    ## For a few tables, the control sample is all 0s, so let us take that into account.
    usable_columns <- colSums(mat_cpm[[table_name]]) > 0
    mat_cpm[[table_name]] <- mat_cpm[[table_name]][, usable_columns]
    mat_cpmlength[[table_name]] <- mat_cpmlength[[table_name]][, usable_columns]

    if (isTRUE(read_table)) {
      mat_cpm[[table_name]] <- try(edgeR::cpm(y=as.matrix(mat_cpm[[table_name]]),
                                              lib.sizes=reads_per_sample), silent=TRUE)
      mat_over_counts[[table_name]] <- mat_over_counts[[table_name]] / reads_per_sample
    } else {
      mat_cpm[[table_name]] <- try(edgeR::cpm(y=as.matrix(mat_cpm[[table_name]]),
                                              lib.sizes=indexes_per_sample), silent=TRUE)
      mat_over_counts[[table_name]] <- mat_over_counts[[table_name]] / indexes_per_sample
    }
    mat_over_countslength[[table_name]] <- mat_over_counts[[table_name]] / sequence_length
    if (class(mat_cpm[[table_name]]) == "try-error") {
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
