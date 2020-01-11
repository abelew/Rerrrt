#' Given a set of sequences containing mutations and identical reads along with
#' a summary of how many reads/index are identical/mutant, classify them as
#' identical indexes, RT indexes, or sequencer indexes.
#'
#' The reads_per_index_summary provides a table of the number of mutant and
#' identical reads associated with every index in the data.  Identical reads are
#' (obviously) ones with 0 mutations.  Sequencer index/reads are when there is 1
#' mutant out of n reads for an index.  RT index/reads are when all reads for an
#' index are mutant.
#'
#' @param chng Set of reads with mutants.
#' @param ident Set of reads which are identical.
#' @param reads_per_index_summary Table containing how many reads/index are
#'   identical/mutant.
#' @param min_sequencer Define how stringent one must be to deem a given
#'   read/index group as a sequencer-based error.  Too small and gets too many
#'   false positives, too large and one gets insufficient data to play with.
#' @param verbose Print some information while running.
#' @return Groups of reads which are RT, identical, and sequencer; along with
#'   summary information.
#' @export
classify_sequences <- function(chng, ident, reads_per_index_summary,
                               min_sequencer=10, verbose=verbose) {
  ## We can merge the information from prune_indexes() to the changed and
  ## identical data in order to have the read summaries along with the mutation
  ## information.  The reason this is important is that we may want eventually
  ## to keep indexes for which 1 read out of n is wild-type (where n is pretty
  ## big).  But for the moment we are just looking for indexes with 0 identical
  ## reads associated with them.
  if (verbose) {
    message("Gathering identical, mutant, and sequencer reads/indexes.")
  }
  chng_with_summary <- merge(chng, reads_per_index_summary,
                             by="index", all.x=TRUE)
  ident_with_summary <- merge(ident, reads_per_index_summary,
                              by="index", all.x=TRUE)
  ## This extracts the subset of indexes with only mutant reads.
  only_mutant_idx <- chng_with_summary[["ident_reads"]] == 0
  rt_mutants <- chng_with_summary[only_mutant_idx, ]
  ## This does the opposite, getting the indexes with only identical reads.
  ## I am not certain that this is something in which we are actually interested.
  only_ident_idx <- ident_with_summary[["mut_reads"]] == 0
  only_ident <- ident_with_summary[only_ident_idx, ]
  ## I think the way to calculate the most stringent set of sequencer-based
  ## errors is to count up the set where chng_with_summary mutant reads is 1 and
  ## the total reads is >= n.
  strict_sequencer_idx <- chng_with_summary[["mut_reads"]] == 1 &
    chng_with_summary[["all_reads"]] >= min_sequencer
  strict_sequencer <- chng_with_summary[strict_sequencer_idx, ]

  all_changed <- expand_mutation_string(chng_with_summary)
  strict_sequencer <- expand_mutation_string(strict_sequencer)
  rt_mutants <- expand_mutation_string(rt_mutants)

  if (verbose) {
    message("Before classification, there are ", nrow(ident_with_summary), " identical reads.")
    message("Before classification, there are ", nrow(all_changed), " reads with mutations.")
    message("After classification, there are ",
            nrow(only_ident), " reads/indexes which are only identical.")
    message("After classification, there are ",
            nrow(strict_sequencer), " reads/indexes which are strictly sequencer.")
    message("After classification, there are ",
            nrow(rt_mutants), " reads/indexes which are deemed from reverse transcriptase.")
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
count_mutation_direction <- function(identical, changed, sequencer, verbose=FALSE) {
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

#' Count up all sorts of fun mutation types in the data!
#'
#' This is my favorite function in this!  It counts up all sorts of fun mutation
#' types in the data and gives back tables of the results.
#'
#' @param strict_rt Indexes from RT errors.
#' @param strict_equencer Indexes from sequencer errors.
#' @param min_indexes The minimum number of indexes required to deem a given
#'   error as 'real'.
#' @param verbose Print some information about the run.
#' @return Tables of the various mutation types.
#' @export
count_mutation_types <- function(strict_rt, strict_sequencer,
                                 min_indexes=3, verbose=FALSE) {

  ## Pull out the number of indexes for each position, type, reference, hit.
  ## We can use this to subset only those mutations with >= x indexes.
  ##  Counting by string.
  num_indexes_by_string <- strict_rt %>%
    group_by(string, position, type, reference, hit, mt,
             transition_transversion, strong_weak) %>%
    summarise("index_count" = n())
  num_indexes_by_sequencer <- strict_sequencer %>%
    group_by(string, position, type, reference, hit, mt,
             transition_transversion, strong_weak) %>%
    summarise("index_count" = n())
  if (is.numeric(min_indexes)) {
    if (verbose) {
      message("Subsetting based on mutations with at least ", min_indexes, " indexes.")
    }
    min_idx <- num_indexes_by_string[["index_count"]] >= min_indexes
    num_indexes_by_string <- num_indexes_by_string[min_idx, ]
    min_idx <- num_indexes_by_sequencer[["index_count"]] >= min_indexes
    num_indexes_by_sequencer <- num_indexes_by_sequencer[min_idx, ]
  }

  ## Get information by position in the template
  ##  Counting by reference position.
  miss_reads_by_position <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(position, reference) %>%
    summarise("position_count" = sum(all_reads))
  miss_indexes_by_position <- num_indexes_by_string %>%
    filter(type == "mis") %>%
    group_by(position, reference) %>%
    summarise("position_count" = sum(index_count))
  miss_sequencer_by_position <- num_indexes_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(position, reference) %>%
    summarise("position_count" = sum(index_count))

  ## Given the reads by position, we should be able to recapitulate the sequence
  ## of the original template.
  nt_in_template <- miss_reads_by_position %>%
    group_by(reference) %>%
    summarise("ref_nt" = n())

  ## Count mutations by identity.
  miss_reads_by_string <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(string) %>%
    summarise("string_count" = sum(all_reads))
  miss_indexes_by_string <- num_indexes_by_string %>%
    filter(type == "mis") %>%
    group_by(string) %>%
    summarise("string_count" = sum(index_count))
  miss_sequencer_by_string <- num_indexes_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(string) %>%
    summarise("string_count" = sum(index_count))

  ## Count mutations by template nucleotide.
  miss_reads_by_refnt <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(reference) %>%
    summarise("nt_count" = sum(all_reads))
  miss_indexes_by_refnt <- num_indexes_by_string %>%
    filter(type == "mis") %>%
    group_by(reference) %>%
    summarise("nt_count" = sum(index_count))
  miss_sequencer_by_refnt <- num_indexes_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(reference) %>%
    summarise("nt_count" = sum(index_count))

  ## Count mutations by hit nucleotide.
  miss_reads_by_hitnt <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(hit) %>%
    summarise("nt_count" = sum(all_reads))
  miss_indexes_by_hitnt <- num_indexes_by_string %>%
    filter(type == "mis") %>%
    group_by(hit) %>%
    summarise("nt_count" = sum(index_count))
  miss_sequencer_by_hitnt <- num_indexes_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(hit) %>%
    summarise("nt_count" = sum(index_count))

  ## Count mutations by overall type A->G, G->A, etc.
  miss_reads_by_type <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(mt) %>%
    summarise("mt_count" = sum(all_reads))
  miss_indexes_by_type <- num_indexes_by_string %>%
    filter(type == "mis") %>%
    group_by(mt) %>%
    summarise("mt_count" = sum(index_count))
  miss_sequencer_by_type <- num_indexes_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(mt) %>%
    summarise("mt_count" = sum(index_count))

  ## Count transitions/transversions
  miss_reads_by_trans <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(transition_transversion) %>%
    summarise("trans_count" = sum(all_reads))
  miss_indexes_by_trans <- num_indexes_by_string %>%
    filter(type == "mis") %>%
    group_by(transition_transversion) %>%
    summarise("trans_count" = sum(index_count))
  miss_sequencer_by_trans <- num_indexes_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(transition_transversion) %>%
    summarise("trans_count" = sum(index_count))

  ## Strong/weak
  miss_reads_by_strength <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(strong_weak) %>%
    summarise("strength_count" = sum(all_reads))
  miss_indexes_by_strength <- num_indexes_by_string %>%
    filter(type == "mis") %>%
    group_by(strong_weak) %>%
    summarise("strength_count" = sum(index_count))
  miss_sequencer_by_strength <- num_indexes_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(strong_weak) %>%
    summarise("strength_count" = sum(index_count))

  ## Insertions by position
  insert_reads_by_position <- strict_rt %>%
    filter(type == "ins") %>%
    group_by(position) %>%
    summarise("position_count" = sum(all_reads))
  insert_indexes_by_position <- num_indexes_by_string %>%
    filter(type == "ins") %>%
    group_by(position) %>%
    summarise("position_count" = sum(index_count))
  insert_sequencer_by_position <- num_indexes_by_sequencer %>%
    filter(type == "ins") %>%
    group_by(position) %>%
    summarise("position_count" = sum(index_count))

  ## Inserts by reference nucleotide
  insert_reads_by_nt <- strict_rt %>%
    filter(type == "ins") %>%
    group_by(hit) %>%
    summarise("insert_count" = sum(all_reads))
  insert_indexes_by_nt <- num_indexes_by_string %>%
    filter(type == "ins") %>%
    group_by(hit) %>%
    summarise("insert_count" = sum(index_count))
  insert_sequencer_by_nt <- num_indexes_by_sequencer %>%
    filter(type == "ins") %>%
    group_by(hit) %>%
    summarise("insert_count" = sum(index_count))

  ## Deletions by position.
  delete_reads_by_position <- strict_rt %>%
    filter(type == "del") %>%
    group_by(position) %>%
    summarise("del_count" = sum(all_reads))
  delete_indexes_by_position <- num_indexes_by_string %>%
    filter(type == "del") %>%
    group_by(position) %>%
    summarise("position_count" = sum(index_count))
  delete_sequencer_by_position <- num_indexes_by_sequencer %>%
    filter(type == "del") %>%
    group_by(position) %>%
    summarise("position_count" = sum(index_count))

  ## Deletions by nucleotide.
  delete_reads_by_nt <- strict_rt %>%
    filter(type == "del") %>%
    group_by(reference) %>%
    summarise("del_count" = sum(all_reads))
  delete_indexes_by_nt <- num_indexes_by_string %>%
    filter(type == "del") %>%
    group_by(reference) %>%
    summarise("del_count" = sum(index_count))
  delete_sequencer_by_nt <- num_indexes_by_sequencer %>%
    filter(type == "del") %>%
    group_by(reference) %>%
    summarise("del_count" = sum(index_count))

  message("Classified mutation strings according to various queries.")
  retlist <- list(
    ## In case we want to normalize by the nt content of the template:
    "template_content" = nt_in_template,
    ## Counting various interesting categories in the data.
    ##   String identity
    "miss_reads_by_string" = miss_reads_by_string,
    "miss_indexes_by_string" = miss_indexes_by_string,
    "miss_sequencer_by_string" = miss_sequencer_by_string,
    ##   Position in the template.
    "miss_reads_by_position" = miss_reads_by_position,
    "miss_indexes_by_position" = miss_indexes_by_position,
    "miss_sequencer_by_position" = miss_sequencer_by_position,
    ##   Nucleotide of the reference.
    "miss_reads_by_refnt" = miss_reads_by_refnt,
    "miss_indexes_by_refnt" = miss_indexes_by_refnt,
    "miss_sequencer_by_refnt" = miss_sequencer_by_refnt,
    ##   Nucleotide of the hit.
    "miss_reads_by_hitnt" = miss_reads_by_hitnt,
    "miss_indexes_by_hitnt" = miss_indexes_by_hitnt,
    "miss_sequencer_by_hitnt" = miss_sequencer_by_hitnt,
    ##   Overall type of mutation.
    "miss_reads_by_type" = miss_reads_by_type,
    "miss_indexes_by_type" = miss_indexes_by_type,
    "miss_sequencer_by_type" = miss_sequencer_by_type,
    ##   Transitions/transversions.
    "miss_reads_by_trans" = miss_reads_by_trans,
    "miss_indexes_by_trans" = miss_indexes_by_trans,
    "miss_sequencer_by_trans" = miss_sequencer_by_trans,
    ##   Strong->Weak, Weak->Weak, etc.
    "miss_reads_by_strength" = miss_reads_by_strength,
    "miss_indexes_by_strength" = miss_indexes_by_strength,
    "miss_sequencer_by_strength" = miss_sequencer_by_strength,
    ##   Insertions by position.
    "insert_reads_by_position" = insert_reads_by_position,
    "insert_indexes_by_position" = insert_indexes_by_position,
    "insert_sequencer_by_position" = insert_sequencer_by_position,
    ##   Insertions by nucleotide.
    "insert_reads_by_nt" = insert_reads_by_nt,
    "insert_indexes_by_nt" = insert_indexes_by_nt,
    "insert_sequencer_by_nt" = insert_sequencer_by_nt,
    ##   Deletions by position.
    "delete_reads_by_position" = delete_reads_by_position,
    "delete_indexes_by_position" = delete_indexes_by_position,
    "delete_sequencer_by_position" = delete_sequencer_by_position,
    ##   Deletions by nucleotide.
    "delete_reads_by_nt" = delete_reads_by_nt,
    "delete_indexes_by_nt" = delete_indexes_by_nt,
    "delete_sequencer_by_nt" = delete_sequencer_by_nt)
  return(retlist)
}
