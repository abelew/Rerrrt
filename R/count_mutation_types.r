#' Count up all sorts of fun mutation types in the data!
#'
#' This is my favorite function in this!  It counts up all sorts of fun mutation
#' types in the data and gives back tables of the results.
#'
#' @param strict_rt Tags from RT errors.
#' @param strict_sequencer Tags from sequencer errors.
#' @param min_tags The minimum number of tags required to deem a given
#'  error as 'real'.
#' @param verbose Print some information about the run.
#' @return Tables of the various mutation types.
#' @export
count_mutation_types <- function(strict_rt, strict_sequencer,
                                 min_tags = 3, verbose = FALSE) {

  ## Pull out the number of tags for each position, type, reference, hit.
  ## We can use this to subset only those mutations with >= x tags.
  ##  Counting by string.

  ## Suck it r cmd check
  string <- position <- type <- reference <- hit <- mt <- NULL
  transition_transversion <- strong_weak <- NULL
  ## previous_nt <- next_nt <- previous_strength <- next_strength <- NULL

  num_tags_by_string <- strict_rt %>%
    group_by(string, position, type, reference, hit, mt,
             transition_transversion, strong_weak) %>%
    summarise("tag_count" = n())
  num_tags_by_sequencer <- strict_sequencer %>%
    group_by(string, position, type, reference, hit, mt,
             transition_transversion, strong_weak) %>%
    summarise("tag_count" = n())
  if (is.numeric(min_tags)) {
    if (verbose) {
      message("Subsetting based on mutations with at least ", min_tags, " tags.")
    }
    min_idx <- num_tags_by_string[["tag_count"]] >= min_tags
    num_tags_by_string <- num_tags_by_string[min_idx, ]
    min_idx <- num_tags_by_sequencer[["tag_count"]] >= min_tags
    num_tags_by_sequencer <- num_tags_by_sequencer[min_idx, ]
  }

  ## Get information by position in the template
  ##  Counting by reference position.
  all_reads <- tag_count <- NULL
  miss_reads_by_position <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(position, reference) %>%
    summarise("position_count" = sum(all_reads))
  miss_tags_by_position <- num_tags_by_string %>%
    filter(type == "mis") %>%
    group_by(position, reference) %>%
    summarise("position_count" = sum(tag_count))
  miss_sequencer_by_position <- num_tags_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(position, reference) %>%
    summarise("position_count" = sum(tag_count))

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
  miss_tags_by_string <- num_tags_by_string %>%
    filter(type == "mis") %>%
    group_by(string) %>%
    summarise("string_count" = sum(tag_count))
  miss_sequencer_by_string <- num_tags_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(string) %>%
    summarise("string_count" = sum(tag_count))

  ## Count mutations by template nucleotide.
  miss_reads_by_refnt <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(reference) %>%
    summarise("nt_count" = sum(all_reads))
  miss_tags_by_refnt <- num_tags_by_string %>%
    filter(type == "mis") %>%
    group_by(reference) %>%
    summarise("nt_count" = sum(tag_count))
  miss_sequencer_by_refnt <- num_tags_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(reference) %>%
    summarise("nt_count" = sum(tag_count))

  ## Count mutations by hit nucleotide.
  miss_reads_by_hitnt <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(hit) %>%
    summarise("nt_count" = sum(all_reads))
  miss_tags_by_hitnt <- num_tags_by_string %>%
    filter(type == "mis") %>%
    group_by(hit) %>%
    summarise("nt_count" = sum(tag_count))
  miss_sequencer_by_hitnt <- num_tags_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(hit) %>%
    summarise("nt_count" = sum(tag_count))

  ## Count mutations by overall type A->G, G->A, etc.
  miss_reads_by_type <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(mt) %>%
    summarise("mt_count" = sum(all_reads))
  miss_tags_by_type <- num_tags_by_string %>%
    filter(type == "mis") %>%
    group_by(mt) %>%
    summarise("mt_count" = sum(tag_count))
  miss_sequencer_by_type <- num_tags_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(mt) %>%
    summarise("mt_count" = sum(tag_count))

  ## Count transitions/transversions
  miss_reads_by_trans <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(transition_transversion) %>%
    summarise("trans_count" = sum(all_reads))
  miss_tags_by_trans <- num_tags_by_string %>%
    filter(type == "mis") %>%
    group_by(transition_transversion) %>%
    summarise("trans_count" = sum(tag_count))
  miss_sequencer_by_trans <- num_tags_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(transition_transversion) %>%
    summarise("trans_count" = sum(tag_count))

  ## Strong/weak
  miss_reads_by_strength <- strict_rt %>%
    filter(type == "mis") %>%
    group_by(strong_weak) %>%
    summarise("strength_count" = sum(all_reads))
  miss_tags_by_strength <- num_tags_by_string %>%
    filter(type == "mis") %>%
    group_by(strong_weak) %>%
    summarise("strength_count" = sum(tag_count))
  miss_sequencer_by_strength <- num_tags_by_sequencer %>%
    filter(type == "mis") %>%
    group_by(strong_weak) %>%
    summarise("strength_count" = sum(tag_count))

  ## Insertions by position
  insert_reads_by_position <- strict_rt %>%
    filter(type == "ins") %>%
    group_by(position) %>%
    summarise("position_count" = sum(all_reads))
  insert_tags_by_position <- num_tags_by_string %>%
    filter(type == "ins") %>%
    group_by(position) %>%
    summarise("position_count" = sum(tag_count))
  insert_sequencer_by_position <- num_tags_by_sequencer %>%
    filter(type == "ins") %>%
    group_by(position) %>%
    summarise("position_count" = sum(tag_count))

  ## Inserts by reference nucleotide
  insert_reads_by_nt <- strict_rt %>%
    filter(type == "ins") %>%
    group_by(hit) %>%
    summarise("insert_count" = sum(all_reads))
  insert_tags_by_nt <- num_tags_by_string %>%
    filter(type == "ins") %>%
    group_by(hit) %>%
    summarise("insert_count" = sum(tag_count))
  insert_sequencer_by_nt <- num_tags_by_sequencer %>%
    filter(type == "ins") %>%
    group_by(hit) %>%
    summarise("insert_count" = sum(tag_count))

  ## Deletions by position.
  delete_reads_by_position <- strict_rt %>%
    filter(type == "del") %>%
    group_by(position) %>%
    summarise("del_count" = sum(all_reads))
  delete_tags_by_position <- num_tags_by_string %>%
    filter(type == "del") %>%
    group_by(position) %>%
    summarise("position_count" = sum(tag_count))
  delete_sequencer_by_position <- num_tags_by_sequencer %>%
    filter(type == "del") %>%
    group_by(position) %>%
    summarise("position_count" = sum(tag_count))

  ## Deletions by nucleotide.
  delete_reads_by_nt <- strict_rt %>%
    filter(type == "del") %>%
    group_by(reference) %>%
    summarise("del_count" = sum(all_reads))
  delete_tags_by_nt <- num_tags_by_string %>%
    filter(type == "del") %>%
    group_by(reference) %>%
    summarise("del_count" = sum(tag_count))
  delete_sequencer_by_nt <- num_tags_by_sequencer %>%
    filter(type == "del") %>%
    group_by(reference) %>%
    summarise("del_count" = sum(tag_count))

  ## Previous nucleotide
  ## rt_reads_by_previous_nt <- strict_rt %>%
  ##   group_by(previous_nt) %>%
  ##   summarise("previous_count" = sum(all_reads))
  ## rt_tags_by_previous_nt <- num_tags_by_string %>%
  ##   group_by(previous_nt) %>%
  ##   summarise("previous_count" = sum(tag_count))
  ## sequencer_by_previous_nt <- num_tags_by_sequencer %>%
  ##   group_by(previous_nt) %>%
  ##   summarise("previous_count" = sum(tag_count))
  ## Next nucleotide
  ## rt_reads_by_next_nt <- strict_rt %>%
  ##   group_by(next_nt) %>%
  ##   summarise("next_count" = sum(all_reads))
  ## rt_tags_by_next_nt <- num_tags_by_string %>%
  ##  group_by(next_nt) %>%
  ##   summarise("next_count" = sum(tag_count))
  ## sequencer_by_next_nt <- num_tags_by_sequencer %>%
  ##   group_by(next_nt) %>%
  ##   summarise("next_count" = sum(tag_count))
  ## Previous strength
  ## rt_reads_by_prevstr <- strict_rt %>%
  ##   group_by(previous_strength) %>%
  ##   summarise("previous_strength" = sum(all_reads))
  ## rt_tags_by_prevstr <- num_tags_by_string %>%
  ##   group_by(previous_strength) %>%
  ##   summarise("previous_strength" = sum(tag_count))
  ## sequencer_by_prevstr <- num_tags_by_sequencer %>%
  ##   group_by(previous_strength) %>%
  ##   summarise("previous_strength" = sum(tag_count))
  ## Next strength
  ## rt_reads_by_nextstr <- strict_rt %>%
  ##   group_by(next_strength) %>%
  ##   summarise("next_strength" = sum(all_reads))
  ## rt_tags_by_nextstr <- num_tags_by_string %>%
  ##   group_by(next_strength) %>%
  ##   summarise("next_strength" = sum(tag_count))
  ## sequencer_by_nextstr <- num_tags_by_sequencer %>%
  ##   group_by(next_strength) %>%
  ##   summarise("next_strength" = sum(tag_count))

  message("Classified mutation strings according to various queries.")
  retlist <- list(
    ## In case we want to normalize by the nt content of the template:
    "template_content" = nt_in_template,
    ## Counting various interesting categories in the data.
    ##   String identity
    "miss_reads_by_string" = miss_reads_by_string,
    "miss_tags_by_string" = miss_tags_by_string,
    "miss_sequencer_by_string" = miss_sequencer_by_string,
    ##   Position in the template.
    "miss_reads_by_position" = miss_reads_by_position,
    "miss_tags_by_position" = miss_tags_by_position,
    "miss_sequencer_by_position" = miss_sequencer_by_position,
    ##   Nucleotide of the reference.
    "miss_reads_by_refnt" = miss_reads_by_refnt,
    "miss_tags_by_refnt" = miss_tags_by_refnt,
    "miss_sequencer_by_refnt" = miss_sequencer_by_refnt,
    ##   Nucleotide of the hit.
    "miss_reads_by_hitnt" = miss_reads_by_hitnt,
    "miss_tags_by_hitnt" = miss_tags_by_hitnt,
    "miss_sequencer_by_hitnt" = miss_sequencer_by_hitnt,
    ##   Overall type of mutation.
    "miss_reads_by_type" = miss_reads_by_type,
    "miss_tags_by_type" = miss_tags_by_type,
    "miss_sequencer_by_type" = miss_sequencer_by_type,
    ##   Transitions/transversions.
    "miss_reads_by_trans" = miss_reads_by_trans,
    "miss_tags_by_trans" = miss_tags_by_trans,
    "miss_sequencer_by_trans" = miss_sequencer_by_trans,
    ##   Strong->Weak, Weak->Weak, etc.
    "miss_reads_by_strength" = miss_reads_by_strength,
    "miss_tags_by_strength" = miss_tags_by_strength,
    "miss_sequencer_by_strength" = miss_sequencer_by_strength,
    ##   Insertions by position.
    "insert_reads_by_position" = insert_reads_by_position,
    "insert_tags_by_position" = insert_tags_by_position,
    "insert_sequencer_by_position" = insert_sequencer_by_position,
    ##   Insertions by nucleotide.
    "insert_reads_by_nt" = insert_reads_by_nt,
    "insert_tags_by_nt" = insert_tags_by_nt,
    "insert_sequencer_by_nt" = insert_sequencer_by_nt,
    ##   Deletions by position.
    "delete_reads_by_position" = delete_reads_by_position,
    "delete_tags_by_position" = delete_tags_by_position,
    "delete_sequencer_by_position" = delete_sequencer_by_position,
    ##   Deletions by nucleotide.
    "delete_reads_by_nt" = delete_reads_by_nt,
    "delete_tags_by_nt" = delete_tags_by_nt,
    "delete_sequencer_by_nt" = delete_sequencer_by_nt)
    ##   Counting by previous/next nucleotides
    ## "rt_reads_by_previous_nt" = rt_reads_by_previous_nt,
    ## "rt_tags_by_previous_nt" = rt_tags_by_previous_nt,
    ## "sequencer_by_previous_nt" = sequencer_by_previous_nt,
    ## "rt_reads_by_next_nt" = rt_reads_by_next_nt,
    ## "rt_tags_by_next_nt" = rt_tags_by_next_nt,
    ## "sequencer_by_next_nt" = sequencer_by_next_nt,
    ## "rt_reads_by_prevstr" = rt_reads_by_prevstr,
    ## "rt_tags_by_prevstr" = rt_tags_by_prevstr,
    ## "sequencer_by_prevstr" = sequencer_by_prevstr,
    ## "rt_reads_by_nextstr" = rt_reads_by_nextstr,
    ## "rt_tags_by_nextstr" = rt_tags_by_nextstr,
    ## "sequencer_by_nextstr" = sequencer_by_nextstr
    ## )

  return(retlist)
}
