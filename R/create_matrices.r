#' Given a samples sheet with some metadata, create a big pile of matrices
#' describing the data.
#'
#' @param sample_sheet xlsx/csv/whatever file of metadata.
#' @param ident_column Column containing the files of reads identical to the
#'  template.
#' @param mut_column Column containing the files of reads not identical to the
#'  template.
#' @param min_reads Filter for the minimum number of reads / tag.
#' @param min_tags Filter for the minimum numer of tags / mutation.
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
#' @param position_offset Change the position numbers to match previous work.
#' @param savefile Save the timeconsuming portion to this file.
#' @param overwrite If the savefile exists, overwrite it?
#' @return List with lots of matrices of the resulting data.
#' @export
create_matrices <- function(sample_sheet="sample_sheets/all_samples.xlsx",
                            ident_column="identtable",
                            mut_column="mutationtable",
                            preprocess_column="preprocessingtags",
                            min_reads=NULL, min_tags=NULL,
                            min_sequencer=10, min_position=NULL,
                            max_position=NULL, max_mutations_per_read=NULL,
                            prune_n=TRUE, excel=NULL, verbose=TRUE,
                            position_offset=-16,
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
  filtered_matrix <- data.frame()
  reads_remaining <- data.frame()
  tags_remaining <- data.frame()
  sample_names <- as.character(meta[["sampleid"]])

  preprocess_tags <- list()
  preprocess_mtrx <- data.frame()
  message("Counting reads/tag.")
  preprocessed <- read_preprocessing(meta, preprocess_column)
  preprocessed_hist_df <- preprocessed[["tag_counts"]]
  preprocessed_hist <- tag_histogram(sample_sheet = sample_sheet, min = 1, max = 10)
  preprocessed_hist <- preprocessed_hist[["histogram"]]

  ## 1. Mutant tags before filtering.
  pre_chng_tag_density_dt <- data.table()
  ## 2. Mutant tags after filtering.
  post_chng_tag_density_dt <- data.table()
  ## 3. Table identical tags before filtering.
  pre_ident_tag_density_dt <- data.table()
  ## 4. identical tags after filtering.
  post_ident_tag_density_dt <- data.table()
  ## 5. Table of both, before filtering.
  pre_tag_density_dt <- data.table()
  ## 6. Table of both, after filtering.
  post_tag_density_dt <- data.table()

  for (s in 1:nrow(meta)) {
    identical <- meta[s, ident_column]
    changed <- meta[s, mut_column]
    preprocess <- meta[s, preprocess_column]
    sname <- as.character(meta[s, "sampleid"])
    if (verbose) {
      message("Starting sample: ", sname, ".")
    }

    samples[[s]] <- quantify_parsed(changed=changed, identical=identical, min_reads=min_reads,
                                    min_tags=min_tags, min_sequencer=min_sequencer,
                                    min_position=min_position, max_position=max_position,
                                    max_mutations_per_read=max_mutations_per_read,
                                    prune_n=prune_n, verbose=verbose,
                                    position_offset=position_offset,
                                    pre_tags=preprocessed[["all_tags"]][[s]])
    if (s == 1) {
      filtered_matrix <- as.data.frame(samples[[s]][["reads_filtered"]])
      colnames(filtered_matrix) <- sname
      reads_remaining <- as.data.frame(samples[[s]][["reads_remaining"]])
      colnames(filtered_matrix) <- sname
      tags_remaining <- as.data.frame(samples[[s]][["tags_remaining"]])
      colnames(tags_remaining) <- sname
    } else {
      filtered_matrix <- cbind(filtered_matrix, samples[[s]][["reads_filtered"]])
      reads_remaining <- cbind(reads_remaining, samples[[s]][["reads_remaining"]])
      tags_remaining <- cbind(tags_remaining, samples[[s]][["tags_remaining"]])
    }

    ## Get tag density before filters, mutant tags
    ## suck it r cmd check
    tag <- NULL
    tag_chng_density <- data.table(
        "tag" = samples[[s]][["starting_chng_tags"]],
        "sample" = sname) %>%
      group_by(tag, sample) %>%
      summarise(tag_count = n(), .groups = "keep")

    ## Get tag density after filters, mutant tags
    changed_dt <- data.table(
        "tag" = samples[[s]][["changed_table"]][["tag"]],
        "sample" = sname) %>%
      group_by(tag, sample) %>%
      summarise(tag_count = n(), .groups = "drop")

    tag_chng_density[["tag"]] <- NULL
    ## Table 1 above, e.g. a table with rows:reads/tag, columns:samples, cells:counts.
    pre_chng_tag_density_dt <- make_tag_table(pre_chng_tag_density_dt,
                                              tag_chng_density, sname)

    changed_dt[["tag"]] <- NULL
    ## Table 2 above.
    post_chng_tag_density_dt <- make_tag_table(post_chng_tag_density_dt,
                                               changed_dt, sname)

    ## Get tag density before filters, indentical tags
    tag_ident_density <- data.table(
        "tag" = samples[[s]][["starting_ident_tags"]],
        "sample" = sname) %>%
      group_by(tag, sample) %>%
      summarise(tag_count = n(), .groups = "drop")
    tag_ident_density[["tag"]] <- NULL
    ## Table 3 above.
    pre_ident_tag_density_dt <- make_tag_table(pre_ident_tag_density_dt,
                                               tag_ident_density, sname)

    ## Get tag density after filter, identical tags
    ident_tags <- data.table(
        "tag" = samples[[s]][["ident_table"]][["tag"]],
        "sample" = sname) %>%
      group_by(tag, sample) %>%
      summarise(tag_count = n(), .groups = "drop")
    ident_tags[["tag"]] <- NULL
    ## Table 4 above.
    post_ident_tag_density_dt <- make_tag_table(post_ident_tag_density_dt,
                                                ident_tags, sname)

    ## Table 5 above, the sum of pre_ident and pre_chng density tables.
    pre_tag_density_dt <- add_tag_densities(pre_tag_density_dt,
                                            pre_ident_tag_density_dt,
                                            pre_chng_tag_density_dt, sname)

    ## Table 6 above, sum of post densities
    post_tag_density_dt <- add_tag_densities(post_tag_density_dt,
                                             post_ident_tag_density_dt,
                                             post_chng_tag_density_dt, sname)
  } ## End iterating over every sample.

  names(samples) <- sample_names
  colnames(filtered_matrix) <- sample_names
  colnames(reads_remaining) <- sample_names
  colnames(tags_remaining) <- sample_names

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
  used_rows <- rowSums(tags_remaining) > 0
  tags_remaining <- tags_remaining[used_rows, ]

  total_ident_reads <- reads_remaining["ident_tag_pruned", ]
  total_mut_reads <- reads_remaining["mut_tag_pruned", ]
  total_reads <- total_ident_reads + total_mut_reads
  reads_remaining <- rbind(reads_remaining, total_reads)

  tables <- c("miss_reads_by_position", "miss_tags_by_position", "miss_sequencer_by_position",
              "miss_reads_by_string", "miss_tags_by_string", "miss_sequencer_by_string",
              "miss_reads_by_refnt", "miss_tags_by_refnt", "miss_sequencer_by_refnt",
              "miss_reads_by_hitnt", "miss_tags_by_hitnt", "miss_sequencer_by_hitnt",
              "miss_reads_by_type", "miss_tags_by_type", "miss_sequencer_by_type",
              "miss_reads_by_trans", "miss_tags_by_trans", "miss_sequencer_by_trans",
              "miss_reads_by_strength", "miss_tags_by_strength", "miss_sequencer_by_strength",
              "insert_reads_by_position", "insert_tags_by_position", "insert_sequencer_by_position",
              "insert_reads_by_nt", "insert_tags_by_nt", "insert_sequencer_by_nt",
              "delete_reads_by_position", "delete_tags_by_position", "delete_sequencer_by_position",
              "delete_reads_by_nt", "delete_tags_by_nt", "delete_sequencer_by_nt")

  ## Now build up the matrices of data where each matrix should have rownames which make
  ## sense for its name and the column names should be by sample.  Then the cells should
  ## be filled in with the data for each sample.
  pre_normalization_data <- matrices_from_tables(tables, samples, verbose=verbose)
  matrices <- pre_normalization_data[["matrices"]]
  tags_per_sample <- pre_normalization_data[["tags_per_sample"]]
  reads_per_sample <- pre_normalization_data[["reads_per_sample"]]

  ## There are a few tables which get out of order when merged.
  ## Let us fix that now.  They are the _by_string tables
  potentially_misordered <- c("miss_reads_by_string", "miss_tags_by_string",
                              "miss_sequencer_by_string")
  for (tab in potentially_misordered) {
    fixme <- matrices[[tab]]
    fixme[["names"]] <- as.character(fixme[["names"]])
    order_idx <- order(fixme[["names"]])
    fixme <- fixme[order_idx, ]
    matrices[[tab]] <- fixme
  }

  ## Perform some normalizations of the data.
  normalized <- normalize_matrices(matrices, reads_per_sample, tags_per_sample)

  retlist <- list(
      "metadata" = meta,
      "samples" = samples,
      "filtered" = filtered_matrix,
      "reads_remaining" = reads_remaining,
      "tags_remaining" = tags_remaining,
      "reads_per_sample" = reads_per_sample,
      "tags_per_sample" = tags_per_sample,
      "matrices" = matrices,
      "matrices_cpm" = normalized[["mat_cpm"]],
      "matrices_cpmlength" = normalized[["mat_cpmlength"]],
      "matrices_counts" = normalized[["mat_over_counts"]],
      "matrices_countslength" = normalized[["mat_over_countslength"]],
      "pre_chng_tag_density_dt" = pre_chng_tag_density_dt,
      "post_chng_tag_density_dt" = post_chng_tag_density_dt,
      "pre_ident_tag_density_dt" = pre_ident_tag_density_dt,
      "post_ident_tag_density_dt" = post_ident_tag_density_dt,
      "pre_tag_density_dt" = pre_tag_density_dt,
      "post_tag_density_dt" = post_tag_density_dt,
      "preprocessed_hist" = preprocessed_hist
  )
  retlist[["plots"]] <- barplot_matrices(retlist)
  if (!is.null(excel)) {
    retlist[["excel"]] <- write_matrices(retlist, excel, verbose=verbose)
  }
  return(retlist)
}
