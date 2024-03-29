#' Quantify the tables of reads identical/different to/from the template sequence.
#'
#' The heavy lifting of locating reads which are/not identical to the template was
#' performed by errrt.pm and returns a series of compressed tables of identical read
#' ids and non-identical insertions, deletions, and mismatches.  This function is
#' intended to read those two files, gather the resulting data, and make
#' some sense of it.
#'
#' @param changed File containing the set of changed reads/tags, by default
#'  named 'step4.txt.xz' by errrt.pm.
#' @param identical File containing the set of identical reads/tags, by
#'  default named 'step2_identical_reads.txt.xz' by errrt.pm.
#' @param min_reads Minimum number of reads for each tag required to include
#'  each tag in the final result.
#' @param min_tags Minimum number of tags required to include a given
#'  mutation in the final result.
#' @param min_sequencer Minimum number of total reads to consider an tag as
#'  associated with a sequencer-based error.
#' @param prune_n Remove mutations of a base to 'N'?
#' @param min_position Filter mutations before this position?
#' @param max_position Filter mutations after this position?
#' @param max_mutations_per_read If a read has more than this number of mutations, drop it.
#' @param reencode Rewrite the tags on a different scale to save memory (6 bit instead of 2)
#' @param pre_tags Dataframe of the number of reads for each tag during preprocessing.
#' @param verbose Print information describing what is happening while this runs.
#' @return List containing a bunch of summary information about the data.
#' @export
quantify_parsed <- function(changed=NULL, identical=NULL, min_reads=NULL,
                            min_tags=NULL, min_sequencer=10, prune_n=TRUE,
                            min_position=24, max_position=176,
                            max_mutations_per_read=10, reencode=FALSE,
                            position_offset=-16,
                            verbose=TRUE, pre_tags=NULL) {
  if (verbose) {
    message("  Reading the file containing mutations: ", changed)
  }
  chng <- readr::read_tsv(changed,
                          ##col_names=c("readid", "tag", "direciton"
                          ##            "position", "type", "ref", "hit"))
                          col_types=c("ncfnccc"), col_names=TRUE)
  ## IMPORTANT: Deletions set the hit column to NA and lead to problems later down.
  na_idx <- is.na(chng[["hit"]])
  chng[na_idx, "hit"] <- ""
  colnames(chng) <- c("readid", "tag", "direction", "position",
                      "type", "reference", "hit")
  chng[["position"]] <- chng[["position"]] + position_offset
  if (verbose) {
    message("  Reading the file containing the identical reads: ", identical)
  }
  ident <- readr::read_tsv(identical,
                           col_types=c("ncc"), col_names=TRUE)
  ## Check for tables which were written before I added column names.
  if (colnames(ident)[1] == "1") {
    ident <- rbind(colnames(ident), ident)
  }
  ## I have been switching all incidences of 'index' to 'tag', but
  ## this one is in the preprocessing tool and so will take a little
  ## more effort.
  colnames(ident) <- c("read_num", "direction", "tag")

  ## The reads in chng and ident should sum to the reads in pre_tags
  ## The problem with doing this in the following fashion lies in the following:
  ## Lots of reads have >1 mutation, as a result table(chng_table) will over-represent the
  ## tags for those reads.  For ident_table, this is fine.
  ##ident_table <- as.data.frame(table(ident[["tag"]]))
  ##colnames(ident_table) <- c("tag", "ident_num")
  ##chng_table <- as.data.frame(table(chng[["tag"]]))
  ##colnames(chng_table) <- c("tag", "chng_num")
  ##query_table <- merge(data.table::as.data.table(ident_table),
  ##                     data.table::as.data.table(chng_table), by = "tag", all = TRUE)
  ##query_table <- merge(data.table::as.data.table(query_table),
  ##                     data.table::as.data.table(pre_tags), by = "tag", all = TRUE)
  ##na_idx <- is.na(query_table)
  ##query_table[na_idx] <- 0

  ## Ditto here vis a vis changing the column name.
  ident_table <- as.data.frame(table(ident[["tag"]]))
  colnames(ident_table) <- c("tag", "ident_num")
  ident_table <- data.table::as.data.table(ident_table)

  ## Instead of using table() as before, I will count by readid.
  chng_table_order <- order(chng[["readid"]])
  chng_table <- chng[chng_table_order, ]
  chng_unique_reads_idx <- !duplicated(chng_table[["readid"]])
  chng_tag_table <- chng_table[chng_unique_reads_idx, ]
  chng_tags_order <- order(chng_tag_table[["tag"]])
  chng_tags <- chng_tag_table[["tag"]][chng_tags_order]
  chng_tag_df <- as.data.frame(table(chng_tag_table[["tag"]]))
  colnames(chng_tag_df) <- c("tag", "chng_num")
  chng_table <- data.table::as.data.table(chng_tag_df)

  query_table <- merge(ident_table, chng_table, by = "tag", all = TRUE)
  na_idx <- is.na(query_table)
  query_table[na_idx] <- 0
  query_table[["num"]] <- query_table[["chng_num"]] + query_table[["ident_num"]]
  query_table <- query_table[, c("tag", "num")]

  ## This was an effort to lower the memory usage.
  ## It did work, but takes significant time.
  ## Changes in how I count tags saved much more memory and so I think this is not needed.
  if (isTRUE(reencode)) {
    message("Re-encoding tags.")
    tag <- NULL
    chng <- data.table::as.data.table(chng)
    chng[, "tag" := dna_to_6bit(tag), by="tag"]
    ident <- data.table::as.data.table(ident)
    ident[, "tag" := dna_to_6bit(tag), by="tag"]
  }
  start_rows <- nrow(chng)
  end_rows <- start_rows

  starting_chng_tags <- chng[["tag"]]
  starting_ident_tags <- ident[["tag"]]

  ## go away r cmd check
  readid <- NULL
  ## Note that this is called 'mutations_by_read' as opposed to 'tags_by_read' because there may
  ## be multiple mutations in a single read.
  mutations_by_read <- chng %>%
    group_by(readid) %>%
    summarise("read_mutants" = n())
  mutations_by_read[["readid"]] <- as.character(mutations_by_read[["readid"]])
  if (verbose) {
    message("Here is a table of the number of mutants observed per read.")
    print(table(mutations_by_read[["read_mutants"]]))
  }

  filtered <- c(0, 0, 0, 0, 0, 0)
  names(filtered) <- c("min_position", "max_position", "prune_n",
                       "max_mutations_per_read", "mut_tags", "ident_tags")
  reads_remaining <- c(nrow(ident), nrow(chng), 0, 0, 0, 0, 0, 0)
  names(reads_remaining) <- c("identicals", "mutants", "min_position",
                              "max_position", "prune_n", "max_mutations_per_read",
                              "ident_tag_pruned", "mut_tag_pruned")

  if (verbose) {
    message("  Counting tags before filtering.")
  }
  chng_tags <- sort(unique(chng[["tag"]]))
  ident_tags <- sort(unique(ident[["tag"]]))
  all_tags <- sort(unique(c(chng_tags, ident_tags)))
  start_tags <- length(all_tags)
  tags_remaining <- c(start_tags, length(chng_tags), length(ident_tags),
                         0, 0)
  names(tags_remaining) <- c("all_start", "all_mutants", "all_identicals",
                                "post_read_filters", "post_pruning")

  if (is.numeric(min_position)) {
    filt <- filter_min_position(chng, min_position=min_position, verbose=verbose)
    filtered["min_position"] <- filt[["removed"]]
    reads_remaining["min_position"] <- filt[["remaining"]]
    chng <- filt[["chng"]]
  }
  if (is.numeric(max_position)) {
    filt <- filter_max_position(chng, max_position=max_position, verbose=verbose)
    filtered["max_position"] <- filt[["removed"]]
    reads_remaining["max_position"] <- filt[["remaining"]]
    chng <- filt[["chng"]]
  }
  if (isTRUE(prune_n)) {
    filt <- filter_ns(chng, verbose=verbose)
    filtered["prune_n"] <- filt[["removed"]]
    reads_remaining["prune_n"] <- filt[["remaining"]]
    chng <- filt[["chng"]]
  }

  min_proportion <- FALSE
  if (is.numeric(min_proportion)) {
    filt <- filter_proportion(chng, ident, min_proportion=min_proportion, verbose=verbose)
    filtered["proportion_mutant"] <- filt[["removed"]]
    reads_remaining["proportion_mutant"] <- filt[["remaining"]]
    chng <- filt[["chng"]]
  }

  if (is.numeric(max_mutations_per_read)) {
    filt <- filter_max_mutations(chng, mutations_by_read,
                                 max_mutations_per_read=max_mutations_per_read,
                                 verbose=verbose)
    filtered["max_mutations_per_read"] <- filt[["removed"]]
    reads_remaining["max_mutations_per_read"] <- filt[["remaining"]]
    chng <- filt[["chng"]]
  }

  post_rows <- nrow(chng)
  sum_removed <- start_rows - post_rows
  pct_removed <- 1 - (post_rows / start_rows)
  sum_pct <- scales::percent(x=pct_removed, accuracy=0.01)
  if (verbose) {
    message("  Mutation data: all filters removed ",
            sum_removed, " reads, or ", sum_pct, ".")
  }

  ## Get some information about the numbers of tags in the data.
  pruned <- prune_tags(chng, ident, min_reads=min_reads, verbose=verbose)
  tags_remaining["post_read_filters"] <- pruned[["unfiltered_tag_num"]]
  tags_remaining["post_pruning"] <- length(pruned[["kept_tags"]])
  ## If it is not 'all', then we want to subset chng and ident accordingly.
  reads_per_tag_summary <- pruned[["tag_summary"]]
  wanted_tags <- pruned[["kept_tags"]]
  if (wanted_tags[1] != "all") {
    chng_ident <- filter_pruned_tags(chng, ident, wanted_tags, min_reads=min_reads,
                                        verbose=verbose)
  }
  chng <- chng_ident[["chng"]]
  ident <- chng_ident[["ident"]]
  filtered["mut_tags"] <- chng_ident[["chng_removed"]]
  filtered["ident_tags"] <- chng_ident[["ident_removed"]]
  reads_remaining["ident_tag_pruned"] <- chng_ident[["ident_remaining"]]
  reads_remaining["mut_tag_pruned"] <- chng_ident[["chng_remaining"]]

  classifications <- classify_sequences(chng, ident, reads_per_tag_summary,
                                        min_sequencer=min_sequencer, verbose=verbose)
  all_changed <- classifications[["all_changed"]]
  strict_identical <- classifications[["strict_identical"]]
  strict_rt <- classifications[["rt_changed"]]
  strict_sequencer <- classifications[["strict_sequencer"]]
  all_identical <- classifications[["all_identical"]]

  mutation_directions <- count_mutation_direction(strict_identical, strict_rt,
                                                  strict_sequencer, verbose=verbose)
  mutation_counts <- count_mutation_types(strict_rt, strict_sequencer,
                                          min_tags=min_tags, verbose=verbose)

  retlist <- list(
    ## Big tables
    "starting_chng_tags" = starting_chng_tags,
    "starting_ident_tags" = starting_ident_tags,
    "reads_filtered" = filtered,
    "reads_remaining" = reads_remaining,
    "tags_remaining" = tags_remaining,
    "mutations_by_read" = mutations_by_read,
    "changed_table" = all_changed,
    "ident_table" = all_identical,
    "strict_rt" = strict_rt,
    "strict_identical" = strict_identical,
    "strict_sequencer" = strict_sequencer,
    ## Tag information.
    "all_tag_count" = pruned[["unfiltered_tag_num"]],
    "all_tag_table" = pruned[["tag_summary"]],
    "filtered_tag_count" = pruned[["filtered_tag_num"]],
    "filtered_tags" = pruned[["kept_tags"]],
    ## Mutations by sequencer direction.
    "identical_reads_by_direction" = mutation_directions[["identical"]],
    "changed_reads_by_direction" = mutation_directions[["changed"]],
    "sequencer_mutations_by_direction" = mutation_directions[["sequencer"]],
    ## In case we want to normalize by the nt content of the template:
    "template_content" = mutation_counts[["nt_in_template"]],
    ## Counting various interesting categories in the data.
    ##   String identity
    "miss_reads_by_string" = mutation_counts[["miss_reads_by_string"]],
    "miss_tags_by_string" = mutation_counts[["miss_tags_by_string"]],
    "miss_sequencer_by_string" = mutation_counts[["miss_sequencer_by_string"]],
    ##   Position in the template.
    "miss_reads_by_position" = mutation_counts[["miss_reads_by_position"]],
    "miss_tags_by_position" = mutation_counts[["miss_tags_by_position"]],
    "miss_sequencer_by_position" = mutation_counts[["miss_sequencer_by_position"]],
    ##   Nucleotide of the reference.
    "miss_reads_by_refnt" = mutation_counts[["miss_reads_by_refnt"]],
    "miss_tags_by_refnt" = mutation_counts[["miss_tags_by_refnt"]],
    "miss_sequencer_by_refnt" = mutation_counts[["miss_sequencer_by_refnt"]],
    ##   Nucleotide of the hit.
    "miss_reads_by_hitnt" = mutation_counts[["miss_reads_by_hitnt"]],
    "miss_tags_by_hitnt" = mutation_counts[["miss_tags_by_hitnt"]],
    "miss_sequencer_by_hitnt" = mutation_counts[["miss_sequencer_by_hitnt"]],
    ##   Overall type of mutation.
    "miss_reads_by_type" = mutation_counts[["miss_reads_by_type"]],
    "miss_tags_by_type" = mutation_counts[["miss_tags_by_type"]],
    "miss_sequencer_by_type" = mutation_counts[["miss_sequencer_by_type"]],
    ##   Transitions/transversions.
    "miss_reads_by_trans" = mutation_counts[["miss_reads_by_trans"]],
    "miss_tags_by_trans" = mutation_counts[["miss_tags_by_trans"]],
    "miss_sequencer_by_trans" = mutation_counts[["miss_sequencer_by_trans"]],
    ##   Strong->Weak, Weak->Weak, etc.
    "miss_reads_by_strength" = mutation_counts[["miss_reads_by_strength"]],
    "miss_tags_by_strength" = mutation_counts[["miss_tags_by_strength"]],
    "miss_sequencer_by_strength" = mutation_counts[["miss_sequencer_by_strength"]],
    ##   Insertions by position.
    "insert_reads_by_position" = mutation_counts[["insert_reads_by_position"]],
    "insert_tags_by_position" = mutation_counts[["insert_tags_by_position"]],
    "insert_sequencer_by_position" = mutation_counts[["insert_sequencer_by_position"]],
    ##   Insertions by nucleotide.
    "insert_reads_by_nt" = mutation_counts[["insert_reads_by_nt"]],
    "insert_tags_by_nt" = mutation_counts[["insert_tags_by_nt"]],
    "insert_sequencer_by_nt" = mutation_counts[["insert_sequencer_by_nt"]],
    ##   Deletions by position.
    "delete_reads_by_position" = mutation_counts[["delete_reads_by_position"]],
    "delete_tags_by_position" = mutation_counts[["delete_tags_by_position"]],
    "delete_sequencer_by_position" = mutation_counts[["delete_sequencer_by_position"]],
    ##   Deletions by nucleotide.
    "delete_reads_by_nt" = mutation_counts[["delete_reads_by_nt"]],
    "delete_tags_by_nt" = mutation_counts[["delete_tags_by_nt"]],
    "delete_sequencer_by_nt" = mutation_counts[["delete_sequencer_by_nt"]])
  return(retlist)
}
