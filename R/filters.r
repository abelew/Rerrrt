filter_min_position <- function(chng, min_position=24, verbose=FALSE) {
  pre_rows <- nrow(chng)
  post_rows <- pre_rows
  if (verbose) {
    message("Mutation data: removing any differences before position: ",
            min_position, ".")
    message("Mutation data: before pruning, there are: ", pre_rows, " reads.")
  }
  min_idx <- chng[["position"]] >= min_position
  chng <- chng[min_idx, ]
  post_rows <- nrow(chng)
  pct <- 1 - (post_rows / pre_rows)
  delta <- pre_rows - post_rows
  pct_diff <- scales::percent(x=pct, accuracy=0.01)
  if (verbose) {
    message("Mutation data: after min-position pruning, there are: ",
            post_rows, " reads: ", delta, " lost or ", pct_diff, ".")
  }
  return(chng)
}

filter_max_position <- function(chng, max_position=176, verbose=FALSE) {
  pre_rows <- nrow(chng)
  post_rows <- pre_rows
  if (verbose) {
    message("Mutation data: removing any differences after position: ",
            max_position, ".")
    message("Mutation data: before pruning, there are: ", pre_rows, " reads.")
  }
  min_idx <- chng[["position"]] <= max_position
  chng <- chng[min_idx, ]
  post_rows <- nrow(chng)
  pct <- 1 - (post_rows / pre_rows)
  delta <- pre_rows - post_rows
  pct_diff <- scales::percent(x=pct, accuracy=0.01)
  if (verbose) {
    message("Mutation data: after max-position pruning, there are: ",
            post_rows, " reads: ", delta, " lost or ", pct_diff, ".")
  }
  return(chng)
}

filter_ns <- function(chng, verbose=FALSE) {
  pre_rows <- nrow(chng)
  if (verbose) {
    message("Mutation data: removing any reads with 'N' as the hit.")
    n_idx <- chng[["hit"]] != "N"
    chng <- chng[n_idx, ]
    post_rows <- nrow(chng)
    pct <- 1 - (post_rows / pre_rows)
    delta <- pre_rows - post_rows
    pct_diff <- scales::percent(x=pct, accuracy=0.01)
    message("Mutation data: after N pruning, there are: ",
            post_rows, " reads: ", delta, " lost or ", pct_diff, ".")
  }
  return(chng)
}

filter_max_mutations <- function(chng, mutation_df,
                                 max_mutations_per_read=10, verbose=FALSE) {
  pre_rows <- nrow(chng)
  message("Mutation data: removing reads with greater than ",
          max_mutations_per_read, " mutations.")
  excluded_readid_idx <- mutation_df[["read_mutants"]] > max_mutations_per_read
  excluded_readids <- mutation_df[excluded_readid_idx, "readid"][["readid"]]
  excluded_idx <- chng[["readid"]] %in% excluded_readids
  chng <- chng[!excluded_idx, ]
  post_rows <- nrow(chng)
  pct <- 1 - (post_rows / pre_rows)
  delta <- pre_rows - post_rows
  pct_diff <- scales::percent(x=pct, accuracy=0.01)
  message("Mutation data: after max_mutation pruning, there are: ",
          post_rows, " reads: ", delta, " lost or ", pct_diff, ".")
  return(chng)
}

filter_pruned_indexes <- function(chng, ident, wanted_indexes, min_reads=3, verbose=FALSE) {
  pre_rows <- nrow(chng)
  pre_ident <- nrow(ident)
  if (verbose) {
    message("All data: removing indexes with fewer than ", min_reads, " reads/index.")
    message("All data: before reads/index pruning, there are: ", pre_rows, " changed reads.")
    message("All data: before reads/index pruning, there are: ", pre_ident, " identical reads.")
  }
  kept_idx <- chng[["index"]] %in% wanted_indexes
  chng <- chng[kept_idx, ]
  kept_idx <- ident[["index"]] %in% wanted_indexes
  ident <- ident[kept_idx, ]
  post_rows <- nrow(chng)
  post_ident <- nrow(ident)
  pct_rows <- scales::percent(x=post_rows / pre_rows, accuracy=0.01)
  pct_ident <- scales::percent(x=post_ident / pre_ident, accuracy=0.01)
  if (verbose) {
    message("All data: after index pruning, there are: ",
            post_rows, " changed reads: ", pct_rows, ".")
    message("All data: after index pruning, there are: ",
            post_ident, " identical reads: ", pct_ident, ".")
  }
  retlist <- list(
    "chng" = chng,
    "ident" = ident)
  return(retlist)
}

#' Summarize the mutant/identical data with respect to the number of
#' reads/index.
#'
#' This should provide a table of how many reads/index are (a)identical, (b)contain
#' mutations, and the sum of (a + b).  If a minimum number of reads is requested
#' (e.g. min_reads is a number), return the list of indexes which have at least
#' that many reads.  This set of indexes may be used in other contexts to limit
#' the data.
#'
#' @param chng The result of read_tsv() on a file containing a mutation table.
#' @param ident The result of read_tsv() on a file containing identical reads.
#' @param min_reads Minimum read / index filter.
#' @param verbose Print information about this while it runs?
#' @return List with the summary of the numbers of reads observed and the
#'   indexes kept.  The list of indexes kept may just be 'all'.
prune_indexes <- function(chng, ident, min_reads=3, verbose=TRUE) {
  ## Get a matrix of how many reads/index are identical to the template.
  sum_ident <- ident %>%
    group_by(index) %>%
    summarise("ident_reads" = n())
  ## Get a matrix of how many reads/index exist for each mutation type.
  sum_mut <- chng[, c("readid", "index")] %>%
    distinct() %>%
    group_by(index) %>%
    summarise("mut_reads" = n())
  ## Convert them to data tables, merge them, and set any NAs to 0.
  sum_all <- merge(sum_ident, sum_mut, by="index", all=TRUE)
  na_idx <- is.na(sum_all)
  sum_all[na_idx] <- 0
  ## Sum up the reads.
  sum_all[["all_reads"]] <- sum_all[["ident_reads"]] + sum_all[["mut_reads"]]
  pre_rows <- nrow(sum_all)
  ## Extract indexes which have >= the minimum number of reads desired.
  kept_indexes <- "all"
  if (is.numeric(min_reads)) {
    if (verbose) {
      message("Before reads/index pruning, there are: ", pre_rows, " indexes in all the data.")
    }
    kept_idx <- sum_all[["all_reads"]] >= min_reads
    kept_indexes <- sum_all[["index"]][kept_idx]
    filtered_index_num <- length(kept_indexes)
    pct <- 1 - (filtered_index_num / pre_rows)
    pct_diff <- scales::percent(x=pct, accuracy=0.01)
    delta <- pre_rows - filtered_index_num
    if (verbose) {
      message("After reads/index pruning, there are: ",
              filtered_index_num, " indexes: ", delta, " lost or ", pct_diff, ".")
    }
  }
  ## Return the table and indexes.
  retlist <- list(
    "unfiltered_index_num" = nrow(sum_all),
    "filtered_index_num" = filtered_index_num,
    "index_summary" = sum_all,
    "kept_indexes" = kept_indexes)
  return(retlist)
}
