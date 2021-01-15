#' Calculate a lower limit error rate for the sequencer.
#'
#' I say 'lower limit' because I assume my set of sequencer-based errors did not
#' detect the full set of error actually from the sequencer.  I chose to take
#' only the set with >= 10 reads / index for which only n (1) read(s) is a
#' mutant.  My hope is that this will avoid false positives, but it will also
#' limit the perceived sequencer error rate.
#'
#' @param data_summary Result from create_matrices().
#' @return Estimate of errors/nucleotide deemed to originate from the sequencer.
#' @export
sequencer_error <- function(data_summary) {
  sequencer_errors <- c()
  sequencer_idents <- c()
  idents <- c()
  only_muts <- c()
  ## Get the number of nucleotides analyzed / read.
  template_length <- nrow(data_summary[["samples"]][[1]][["miss_sequencer_by_position"]])

  for (s in 1:length(data_summary[["samples"]])) {
    sname <- names(data_summary[["samples"]])[s]
    ## Each row of the strict_sequencer table represents 1 error / index for n reads.
    sequencer_errors <- c(nrow(data_summary[["samples"]][[sname]][["strict_sequencer"]]), sequencer_errors)
    ## Therefore for each of those errors, the ident_reads column tells us how many nucleotide in the identical reads for those indexes were read.
    sequencer_idents <- c(sum(data_summary[["samples"]][[sname]][["strict_sequencer"]][["ident_reads"]]), sequencer_idents)
    ## We also have the population of identical reads.
    idents <- c(sum(data_summary[["samples"]][[sname]][["ident_table"]][["all_reads"]]), idents)
    ## Reads from indexes with only mutations
    only_muts <- c(sum(data_summary[["samples"]][[sname]][["only_mutants"]][["all_reads"]]), only_muts)
  }

  ## All nucleotides read in this set:
  all_nucleotides <- (sequencer_errors + sequencer_idents + idents + only_muts) * template_length
  error_nucleotides <- sequencer_errors + only_muts

  sequencer_rate <- error_nucleotides / all_nucleotides
  return(sequencer_rate)
}
