#' Given a set of sequences containing mutations and identical reads along with
#' a summary of how many reads/tag are identical/mutant, classify them as
#' identical tags, RT tags, or sequencer tags.
#'
#' The reads_per_tag_summary provides a table of the number of mutant and
#' identical reads associated with every tag in the data.  Identical reads are
#' (obviously) ones with 0 mutations.  Sequencer tag/reads are when there is 1
#' mutant out of n reads for an tag.  RT tag/reads are when all reads for an
#' tag are mutant.
#'
#' @param chng Set of reads with mutants.
#' @param ident Set of reads which are identical.
#' @param reads_per_tag_summary Table containing how many reads/tag are
#'  identical/mutant.
#' @param min_sequencer Define how stringent one must be to deem a given
#'  read/tag group as a sequencer-based error.  Too small and gets too many
#'  false positives, too large and one gets insufficient data to play with.
#' @param verbose Print some information while running.
#' @return Groups of reads which are RT, identical, and sequencer; along with
#'   summary information.
#' @export
classify_sequences <- function(chng, ident, reads_per_tag_summary,
                               min_sequencer = 10, verbose = verbose) {
  ## We can merge the information from prune_tags() to the changed and
  ## identical data in order to have the read summaries along with the mutation
  ## information.  The reason this is important is that we may want eventually
  ## to keep tags for which 1 read out of n is wild-type (where n is pretty
  ## big).  But for the moment we are just looking for tags with 0 identical
  ## reads associated with them.
  if (verbose) {
    message("Gathering identical, mutant, and sequencer reads/tags.")
  }
  chng_with_summary <- merge(chng, reads_per_tag_summary,
                             by = "tag", all.x = TRUE)
  ident_with_summary <- merge(ident, reads_per_tag_summary,
                              by = "tag", all.x = TRUE)
  ## This extracts the subset of tags with only mutant reads.
  only_mutant_idx <- chng_with_summary[["ident_reads"]] == 0
  rt_mutants <- chng_with_summary[only_mutant_idx, ]
  ## This does the opposite, getting the tags with only identical reads.
  ## I am not certain that this is something in which we are actually interested.
  only_ident_idx <- ident_with_summary[["mut_reads"]] == 0
  only_ident <- ident_with_summary[only_ident_idx, ]
  ## I think the way to calculate the most stringent set of sequencer-based
  ## errors is to count up the set where chng_with_summary mutant reads is 1 and
  ## the total reads is >= n.
  strict_sequencer_idx <- chng_with_summary[["mut_reads"]] == 1 &
    chng_with_summary[["all_reads"]] >= min_sequencer
  strict_sequencer <- chng_with_summary[strict_sequencer_idx, ]

  ## I was thinking to potentially quantify the reads which are neither
  ## deemed RT nor sequencer.
  ## uncharacterized_idx <- !only_mutant_idx
  ## uncharacterized_mutants <- chng_with_summary[uncharacterized_idx, ]
  ## uncharacterized_idx <- chng_with_summary[["mut_reads"]] != 1
  ## uncharacterized_mutants <- chng_with_summary[uncharacterized_idx, ]

  all_changed <- expand_mutation_string(chng_with_summary)
  strict_sequencer <- expand_mutation_string(strict_sequencer)
  rt_mutants <- expand_mutation_string(rt_mutants)

  if (verbose) {
    message("Before classification, there are ", nrow(ident_with_summary), " identical reads.")
    message("Before classification, there are ", nrow(all_changed), " reads with mutations.")
    message("After classification, there are ",
            nrow(only_ident), " reads/tags which are only identical.")
    message("After classification, there are ",
            nrow(strict_sequencer), " reads/tags which are strictly sequencer.")
    message("After classification, there are ",
            nrow(rt_mutants), " reads/tags which are consistently repeated.")
  }

  retlist <- list(
    "all_identical" = ident_with_summary,
    "all_changed" = all_changed,
    "rt_changed" = rt_mutants,
    "strict_sequencer" = strict_sequencer,
    "strict_identical" = only_ident)
  return(retlist)
}

#' Count up how many reads are on the forward and reverse strands.
#'
#' @param identical Set of reads identical to the template.
#' @param changed Reads deemed to be from the RT.
#' @param sequencer Reads deemed to be from the sequencer.
#' @param verbose Print some information while running.
#' @return List of how many forward and reverse reads are in the data.
#' @export
count_mutation_direction <- function(identical, changed, sequencer, verbose = FALSE) {
  ## Suck it, R CMD CHECK
  direction <- ident_reads <- mut_reads <- NULL
  ident_by_direction <- identical %>%
    group_by(direction) %>%
    summarise("dir_reads" = sum(ident_reads))
  changed_by_direction <- changed[, c("readid", "direction", "mut_reads")] %>%
    distinct() %>%
    group_by(direction) %>%
    summarise("dir_reads" = sum(mut_reads))
  sequencer_by_direction <- sequencer[, c("readid", "direction", "mut_reads")] %>%
    distinct() %>%
    group_by(direction) %>%
    summarise("dir_reads" = sum(mut_reads))

  forward_reads <- ident_by_direction[1, "dir_reads"] +
    changed_by_direction[1, "dir_reads"] +
    sequencer_by_direction[1, "dir_reads"]
  reverse_reads <- ident_by_direction[2, "dir_reads"] +
    changed_by_direction[2, "dir_reads"] +
    sequencer_by_direction[2, "dir_reads"]
  if (verbose) {
    message("Counted by direction: ", forward_reads,
            " forward reads and ", reverse_reads, " reverse_reads.")
  }
  retlist <- list(
    "forward" = forward_reads,
    "reverse" = reverse_reads,
    "identical" = ident_by_direction,
    "changed" = changed_by_direction,
    "sequencer" = sequencer_by_direction)
  return(retlist)
}
