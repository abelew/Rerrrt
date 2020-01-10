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
                            prune_n=TRUE, verbose=TRUE) {
  meta <- hpgltools::extract_metadata(metadata=sample_sheet)
  samples <- list()
  for (s in 1:nrow(meta)) {
    if (verbose) {
      message("Starting sample: ", s, ".")
    }
    identical <- meta[s, ident_column]
    changed <- meta[s, mut_column]
    samples[[s]] <- quantify_parsed(changed=changed, identical=identical, min_reads=min_reads,
                                    min_indexes=min_indexes, min_sequencer=min_sequencer,
                                    min_position=min_position, max_position=max_position,
                                    max_mutations_per_read=max_mutations_per_read,
                                    prune_n=prune_n, verbose=verbose)
  }
  names(samples) <- meta[["sampleid"]]

  tables <- c("miss_reads_by_position", "miss_indexes_by_position", "miss_sequencer_by_position",
              "miss_reads_by_string", "miss_indexes_by_string", "miss_sequencer_by_string",
              "miss_reads_by_ref_nt", "miss_indexes_by_ref_nt", "miss_sequencer_by_ref_nt",
              "miss_reads_by_hit_nt", "miss_indexes_by_hit_nt", "miss_sequencer_by_hit_nt",
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

  normalization_data <- normalize_matrices(matrices, reads_per_sample, indexes_per_sample)
  normalized <- normalization_data[["normalized"]]
  normalized_by_counts <- normalization_data[["normalized_by_counts"]]
  matrices_by_counts <- normalization_data[["matrices_by_counts"]]

  retlist <- list(
    "samples" = samples,
    "reads_per_sample" = reads_per_sample,
    "indexes_per_sample" = indexes_per_sample,
    "matrices" = matrices,
    "matrices_by_counts" = matrices_by_counts,
    "normalized" = normalized,
    "normalized_by_counts" = normalized_by_counts)
  return(retlist)
}

#' Combine the columns describing a mutation into a single column and categorize.
#'
#' Given a table with columns including the position, mutation type, reference
#' nt at that position, and product nt at that position; create a single column
#' from it, standardize it, and make some categories describing each mutation.
#' These columns currently include:  'mt': X_Y telling that this is a mutation
#' from X to Y, 'transition_transversion': this is a transition or transversion
#' mismatch, 'strong_weak': this mismatch went from strong->weak, weak->strong,
#' weak->weak, or strong->strong.  undef is used in the case of indels.
#'
#' @param t Table from readr.
#' @return The same table with some new columns describing the mutations therein.
expand_mutation_string <- function(t) {
  t[["string"]] <- paste0(t[["position"]], "_",
                          t[["type"]], "_",
                          t[["reference"]], "_",
                          t[["hit"]])
  t[["string"]] <- gsub(x=t[["string"]],
                        pattern="^([[:digit:]])_",
                        replacement="00\\1_")
  t[["string"]] <- gsub(x=t[["string"]],
                        pattern="^([[:digit:]][[:digit:]])_",
                        replacement="0\\1_")
  t[["string"]] <- as.factor(t[["string"]])
  t[["mt"]] <- paste0(t[["reference"]], "_", t[["hit"]])
  t[["transition_transversion"]] <- "undef"
  transition_idx <- t[["mt"]] == "A_G" | t[["mt"]] == "G_A" |
    t[["mt"]] == "T_C" | t[["mt"]] == "C_T"
  t[transition_idx, "transition_transversion"] <- "transition"
  transversion_idx <- t[["mt"]] == "A_T" | t[["mt"]] == "A_C" |
    t[["mt"]] == "T_A" | t[["mt"]] == "C_A" |
    t[["mt"]] == "G_C" | t[["mt"]] == "G_T" |
    t[["mt"]] == "C_G" | t[["mt"]] == "T_G"
  t[transversion_idx, "transition_transversion"] <- "transversion"
  t[["transition_transversion"]] <- as.factor(t[["transition_transversion"]])
  t[["strong_weak"]] <- "undef"
  weak_strong_idx <- t[["mt"]] == "A_G" | t[["mt"]] == "T_G" |
    t[["mt"]] == "A_C" | t[["mt"]] == "T_C"
  t[weak_strong_idx, "strong_weak"] <- "weak_strong"
  strong_weak_idx <- t[["mt"]] == "G_A" | t[["mt"]] == "G_T" |
    t[["mt"]] == "C_A" | t[["mt"]] == "C_T"
  t[strong_weak_idx, "strong_weak"] <- "strong_weak"
  weak_weak_idx <- t[["mt"]] == "A_T" | t[["mt"]] == "T_A"
  t[weak_weak_idx, "strong_weak"] <- "weak_weak"
  strong_strong_idx <- t[["mt"]] == "G_C" | t[["mt"]] == "C_G"
  t[strong_strong_idx, "strong_weak"] <- "strong_strong"
  t[["strong_weak"]] <- as.factor(t[["strong_weak"]])
  return(t)
}

#' Quantify the tables of reads identical/different to/from the template sequence.
#'
#' The heavy lifting of locating reads which are/not identical to the template was
#' performed by errrt.pm and returns a series of compressed tables of identical read
#' ids and non-identical insertions, deletions, and mismatches.  This function is
#' intended to read those two files, gather the resulting data, and make
#' some sense of it.
#'
#' @param changed File containing the set of changed reads/indexes, by default
#'   named 'step4.txt.xz' by errrt.pm.
#' @param identical File containing the set of identical reads/indexes, by
#'   default named 'step2_identical_reads.txt.xz' by errrt.pm.
#' @param min_reads Minimum number of reads for each index required to include
#'   each index in the final result.
#' @param min_indexes Minimum number of indexes required to include a given
#'   mutation in the final result.
#' @param min_sequencer Minimum number of total reads to consider an index as
#'   associated with a sequencer-based error.
#' @param prune_n Remove mutations of a base to 'N'?
#' @param min_position Filter mutations before this position?
#' @param max_position Filter mutations after this position?
#' @param verbose Print information describing what is happening while this runs.
#' @return List containing a bunch of summary information about the data.
#' @export
quantify_parsed <- function(changed=NULL, identical=NULL, min_reads=NULL,
                            min_indexes=NULL, min_sequencer=10, prune_n=TRUE,
                            min_position=24, max_position=176,
                            max_mutations_per_read=NULL, verbose=TRUE) {
  if (verbose) {
    message("Reading the file containing mutations: ", changed)
  }
  chng <- readr::read_tsv(changed,
                          ##col_names=c("readid", "index", "direciton"
                          ##            "position", "type", "ref", "hit"))
                          col_types=c("ncfnccc"), col_names=TRUE)
  if (verbose) {
    message("Reading the file containing the identical reads: ", identical)
  }
  ident <- readr::read_tsv(identical,
                           col_types=c("ncc"),
                           col_names=c("readid", "direction", "index"))
  chng <- data.table::as.data.table(chng)
  ident <- data.table::as.data.table(ident)
  start_rows <- nrow(chng)
  end_rows <- start_rows

  mutations_by_read <- chng %>%
    group_by(readid) %>%
    summarise("read_mutants" = n())
  mutations_by_read[["readid"]] <- as.character(mutations_by_read[["readid"]])

  if (is.numeric(min_position)) {
    chng <- filter_min_position(chng, min_position=min_position, verbose=verbose)
  }
  if (is.numeric(max_position)) {
    chng <- filter_max_position(chng, max_position=max_position, verbose=verbose)
  }
  if (isTRUE(prune_n)) {
    chng <- filter_ns(chng, verbose=verbose)
  }
  if (is.numeric(max_mutations_per_read)) {
    chng <- filter_max_mutations(chng, mutations_by_read,
                                 max_mutations_per_read=max_mutations_per_read,
                                 verbose=verbose)
  }

  post_rows <- nrow(chng)
  sum_removed <- start_rows - post_rows
  pct_removed <- 1 - (post_rows / start_rows)
  sum_pct <- scales::percent(x=pct_removed, accuracy=0.01)
  if (verbose) {
    message("Mutation data: all filters removed ",
            sum_removed, " reads, or ", sum_pct, ".")
  }

  ## Get some information about the numbers of indexes in the data.
  if (verbose) {
    message("All data: gathering information about the indexes observed, this is slow.")
  }
  pruned <- prune_indexes(chng, ident, min_reads=min_reads, verbose=verbose)

  ## If it is not 'all', then we want to subset chng and ident accordingly.
  reads_per_index_summary <- pruned[["index_summary"]]
  wanted_indexes <- pruned[["kept_indexes"]]
  if (wanted_indexes[1] != "all") {
    chng_ident <- filter_pruned_indexes(chng, ident, wanted_indexes, min_reads=min_reads,
                                        verbose=verbose)
  }
  chng <- chng_ident[["chng"]]
  ident <- chng_ident[["ident"]]

  classifications <- classify_sequences(chng, ident, reads_per_index_summary,
                                        min_sequencer=min_sequencer, verbose=verbose)
  all_changed <- classifications[["all_changed"]]
  strict_identical <- classifications[["strict_identical"]]
  strict_rt <- classifications[["rt_changed"]]
  strict_sequencer <- classifications[["strict_sequencer"]]
  all_identical <- classifications[["all_identical"]]

  mutation_directions <- count_mutation_direction(strict_identical, strict_rt,
                                                  strict_sequencer, verbose=verbose)
  mutation_counts <- count_mutation_types(strict_rt, strict_sequencer,
                                          min_indexes=min_indexes, verbose=verbose)

  retlist <- list(
    ## Big tables
    "mutations_by_read" = mutations_by_read,
    "changed_table" = all_changed,
    "ident_table" = all_identical,
    "strict_rt" = strict_rt,
    "strict_identical" = strict_identical,
    "strict_sequencer" = strict_sequencer,
    ## Index information.
    "all_index_count" = pruned[["unfiltered_index_num"]],
    "all_index_table" = pruned[["index_summary"]],
    "filtered_index_count" = pruned[["filtered_index_num"]],
    "filtered_indexes" = pruned[["kept_indexes"]],
    ## Mutations by sequencer direction.
    "identical_reads_by_direction" = mutation_directions[["identical"]],
    "changed_reads_by_direction" = mutation_directions[["changed"]],
    "sequencer_mutations_by_direction" = mutation_directions[["sequencer"]],
    ## In case we want to normalize by the nt content of the template:
    "template_content" = mutation_counts[["nt_in_template"]],
    ## Counting various interesting categories in the data.
    ##   String identity
    "miss_reads_by_string" = mutation_counts[["miss_reads_by_string"]],
    "miss_indexes_by_string" = mutation_counts[["miss_indexes_by_string"]],
    "miss_sequencer_by_string" = mutation_counts[["miss_sequencer_by_string"]],
    ##   Position in the template.
    "miss_reads_by_position" = mutation_counts[["miss_reads_by_position"]],
    "miss_indexes_by_position" = mutation_counts[["miss_indexes_by_position"]],
    "miss_sequencer_by_position" = mutation_counts[["miss_sequencer_by_position"]],
    ##   Nucleotide of the reference.
    "miss_reads_by_ref_nt" = mutation_counts[["miss_reads_by_refnt"]],
    "miss_indexes_by_ref_nt" = mutation_counts[["miss_indexes_by_refnt"]],
    "miss_sequencer_by_ref_nt" = mutation_counts[["miss_sequencer_by_refnt"]],
    ##   Nucleotide of the hit.
    "miss_reads_by_hit_nt" = mutation_counts[["miss_reads_by_hitnt"]],
    "miss_indexes_by_hit_nt" = mutation_counts[["miss_indexes_by_hitnt"]],
    "miss_sequencer_by_hit_nt" = mutation_counts[["miss_sequencer_by_hitnt"]],
    ##   Overall type of mutation.
    "miss_reads_by_type" = mutation_counts[["miss_reads_by_type"]],
    "miss_indexes_by_type" = mutation_counts[["miss_indexes_by_type"]],
    "miss_sequencer_by_type" = mutation_counts[["miss_sequencer_by_type"]],
    ##   Transitions/transversions.
    "miss_reads_by_trans" = mutation_counts[["miss_reads_by_trans"]],
    "miss_indexes_by_trans" = mutation_counts[["miss_indexes_by_trans"]],
    "miss_sequencer_by_trans" = mutation_counts[["miss_sequencer_by_trans"]],
    ##   Strong->Weak, Weak->Weak, etc.
    "miss_reads_by_strength" = mutation_counts[["miss_reads_by_strength"]],
    "miss_indexes_by_strength" = mutation_counts[["miss_indexes_by_strength"]],
    "miss_sequencer_by_strength" = mutation_counts[["miss_sequencer_by_strength"]],
    ##   Insertions by position.
    "insert_reads_by_position" = mutation_counts[["insert_reads_by_position"]],
    "insert_indexes_by_position" = mutation_counts[["insert_indexes_by_position"]],
    "insert_sequencer_by_position" = mutation_counts[["insert_sequencer_by_position"]],
    ##   Insertions by nucleotide.
    "insert_reads_by_nt" = mutation_counts[["insert_reads_by_nt"]],
    "insert_indexes_by_nt" = mutation_counts[["insert_indexes_by_nt"]],
    "insert_sequencer_by_nt" = mutation_counts[["insert_sequencer_by_nt"]],
    ##   Deletions by position.
    "delete_reads_by_position" = mutation_counts[["delete_reads_by_position"]],
    "delete_indexes_by_position" = mutation_counts[["delete_indexes_by_position"]],
    "delete_sequencer_by_position" = mutation_counts[["delete_sequencer_by_position"]],
    ##   Deletions by nucleotide.
    "delete_reads_by_nt" = mutation_counts[["delete_reads_by_nt"]],
    "delete_indexes_by_nt" = mutation_counts[["delete_indexes_by_nt"]],
    "delete_sequencer_by_nt" = mutation_counts[["delete_sequencer_by_nt"]])
  return(retlist)
}

matrices_from_tables <- function(tables, samples, verbose=FALSE) {
  matrices <- list()
  indexes_per_sample <- c()
  reads_per_sample <- c()
  for (t in tables) {
    inc <- 0
    col_names <- c()
    if (verbose) {
      message("Making a matrix of ", t, ".")
    }
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

normalize_matrices <- function(matrices, reads_per_sample, indexes_per_sample) {
  ## I still need to figure out how I want to normalize these numbers.
  normalized <- matrices
  normalized_by_counts <- matrices
  matrices_by_counts <- matrices
  for (t in 1:length(names(normalized))) {
    table_name <- names(normalized)[t]
    if (nrow(matrices[[table_name]]) == 0) {
      message("Skipping table: ", table_name)
      next
    }
    read_table <- stringr::str_detect(string=table_name, pattern="reads")

    ## For a few tables, the control sample is all 0s, so let us take that into account.
    usable_columns <- colSums(normalized[[table_name]]) > 0
    normalized[[table_name]] <- normalized[[table_name]][, usable_columns]
    normalized_by_counts[[table_name]] <- normalized_by_counts[[table_name]][, usable_columns]

    if (isTRUE(read_table)) {
      normalized[[table_name]] <- try(edgeR::cpm(y=as.matrix(normalized[[table_name]]),
                                                 lib.sizes=reads_per_sample), silent=TRUE)
      matrices_by_counts[[table_name]] <- matrices_by_counts[[table_name]] / reads_per_sample
    } else {
      normalized[[table_name]] <- try(edgeR::cpm(y=as.matrix(normalized[[table_name]]),
                                                 lib.sizes=indexes_per_sample), silent=TRUE)
      matrices_by_counts[[table_name]] <- matrices_by_counts[[table_name]] / indexes_per_sample
    }
    if (class(normalized[[table_name]]) == "try-error") {
      message("Normalization failed for table: ", table_name, ", leaving it alone.")
      normalized[[table_name]] <- matrices[[table_name]]
    } else {
      if (isTRUE(read_table)) {
        normalized_by_counts[[table_name]] <- normalized[[table_name]] /
          reads_per_sample[usable_columns]
      } else {
        normalized_by_counts[[table_name]] <- normalized[[table_name]] /
          indexes_per_sample[usable_columns]
      }
    }
  }

  retlist <- list(
    "normalized" = normalized,
    "normalized_by_counts" = normalized_by_counts,
    "matrices_by_counts" = matrices_by_counts)
  return(retlist)
}
