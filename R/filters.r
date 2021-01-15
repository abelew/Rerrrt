#' Filter away reads before position x.
#'
#' The error rate in the data is crazy-high in the immediate vicinity of the
#' RT/PCR primers.  Thus we exclude reads with mutations before this position.
#'
#' @param chng Data set of reads which are not identical to the template.
#' @param min_position Position before which to filter.
#' @param verbose Print some information while running?
#' @return The remaining reads and a little summary information.
#' @export
filter_min_position <- function(chng, min_position=24, verbose=FALSE) {
  pre_rows <- nrow(chng)
  post_rows <- pre_rows
  if (verbose) {
    message("   Mutation data: removing any differences before position: ",
            min_position, ".")
    message("   Mutation data: before pruning, there are: ", pre_rows, " reads.")
  }
  min_idx <- chng[["position"]] >= min_position
  chng <- chng[min_idx, ]
  post_rows <- nrow(chng)
  pct <- 1 - (post_rows / pre_rows)
  delta <- pre_rows - post_rows
  pct_diff <- scales::percent(x=pct, accuracy=0.01)
  if (verbose) {
    message("   Mutation data: after min-position pruning, there are: ",
            post_rows, " reads: ", delta, " lost or ", pct_diff, ".")
  }
  retlist <- list(
    "start" = pre_rows,
    "removed" = delta,
    "remaining" = post_rows,
    "chng" = chng)
  return(retlist)
}

#' Filter away reads after position x.
#'
#' The error rate in the data is crazy-high in the immediate vicinity of the
#' RT/PCR primers.  Thus we exclude reads with mutations before this position.
#'
#' @param chng Data set of reads which are not identical to the template.
#' @param max_position Position after which to filter.
#' @param verbose Print some information while running?
#' @return The remaining reads and a little summary information.
#' @export
filter_max_position <- function(chng, max_position=176, verbose=FALSE) {
  pre_rows <- nrow(chng)
  post_rows <- pre_rows
  if (verbose) {
    message("   Mutation data: removing any differences after position: ",
            max_position, ".")
    message("   Mutation data: before pruning, there are: ", pre_rows, " reads.")
  }
  min_idx <- chng[["position"]] <= max_position
  chng <- chng[min_idx, ]
  post_rows <- nrow(chng)
  pct <- 1 - (post_rows / pre_rows)
  delta <- pre_rows - post_rows
  pct_diff <- scales::percent(x=pct, accuracy=0.01)
  if (verbose) {
    message("   Mutation data: after max-position pruning, there are: ",
            post_rows, " reads: ", delta, " lost or ", pct_diff, ".")
  }
  retlist <- list(
    "start" = pre_rows,
    "removed" = delta,
    "remaining" = post_rows,
    "chng" = chng)
  return(retlist)
}

#' Filter away reads which contain Ns.
#'
#' It is unlikely that we will make much sense out of reads with N.  So toss
#' em!
#'
#' @param chng Data set of reads which are not identical to the template.
#' @param verbose Print some information while running?
#' @return The remaining reads and a little summary information.
#' @export
filter_ns <- function(chng, verbose=FALSE) {
  pre_rows <- nrow(chng)
  if (verbose) {
    message("   Mutation data: removing any reads with 'N' as the hit.")
    n_idx <- chng[["hit"]] != "N"
    chng <- chng[n_idx, ]
    post_rows <- nrow(chng)
    pct <- 1 - (post_rows / pre_rows)
    delta <- pre_rows - post_rows
    pct_diff <- scales::percent(x=pct, accuracy=0.01)
    message("   Mutation data: after N pruning, there are: ",
            post_rows, " reads: ", delta, " lost or ", pct_diff, ".")
  }
  retlist <- list(
    "start" = pre_rows,
    "removed" = delta,
    "remaining" = post_rows,
    "chng" = chng)
  return(retlist)
}

#' Filter away reads which have > n mutations in them.
#'
#' A subset of the reads have a weirdly large number of differences to the
#' template.  These add confusing noise to the data, however removing them may
#' take away meaningful noise from the data.  So this function may prove to be
#' awesome or terrible, depending on one's perspective.
#'
#' @param chng Data set of reads which are not identical to the template.
#' @param mutation_df Data frame of readids and their number of mutations.
#' @param max_mutations_per_read What it says on the tin.
#' @param verbose Print some information while running?
#' @return The remaining reads and a little summary information.
#' @export
filter_max_mutations <- function(chng, mutation_df,
                                 max_mutations_per_read=10, verbose=FALSE) {
  pre_rows <- nrow(chng)
  message("   Mutation data: removing reads with greater than ",
          max_mutations_per_read, " mutations.")
  excluded_readid_idx <- mutation_df[["read_mutants"]] > max_mutations_per_read
  excluded_readids <- mutation_df[excluded_readid_idx, "readid"][["readid"]]
  excluded_idx <- chng[["readid"]] %in% excluded_readids
  chng <- chng[!excluded_idx, ]
  post_rows <- nrow(chng)
  pct <- 1 - (post_rows / pre_rows)
  delta <- pre_rows - post_rows
  pct_diff <- scales::percent(x=pct, accuracy=0.01)
  message("   Mutation data: after max_mutation pruning, there are: ",
          post_rows, " reads: ", delta, " lost or ", pct_diff, ".")
  retlist <- list(
    "start" = pre_rows,
    "removed" = delta,
    "remaining" = post_rows,
    "chng" = chng)
  return(retlist)
}

#' Filter out thing which are of an insufficient proportion.
#'
#' This is not currently used, but was something I wanted to consider.
#'
#' @param chng Table of changed reads.
#' @param ident Table of identical reads.
#' @param min_proportion Minimum proportion of changed/identical reads/index.
#' @param verbose Print while running.
filter_proportion <- function(chng, ident, min_proportion=0.5, verbose=FALSE) {
  pre_rows <- nrow(chng)
  pre_ident <- nrow(ident)
  if (verbose) {
    message("All data: removing indexes with less than ", min_proportion, " mutant/identical reads/index.")
    message("All data: before proportion filtering, there are: ", pre_rows, " changed reads.")
    message("All data: before proportion filtering, there are: ", pre_ident, " identical reads.")
  }

  chng_table <- table(chng[["index"]])
  gt_zero <- ident[["index"]] %in% names(chng_table)
  ident_subset <- ident[gt_zero, "index"]
  ident_table <- table(ident_subset)
  chng_table_dt <- data.table::as.data.table(chng_table)
  colnames(chng_table_dt) <- c("index", "chng")
  ident_table_dt <- data.table::as.data.table(ident_table)
  colnames(ident_table_dt) <- c("index", "ident")
  both_dt <- merge(chng_table_dt, ident_table_dt, by="index")
  both_dt[["prop"]] <- both_dt[["chng"]] / (both_dt[["chng"]] + both_dt[["ident"]])
  kept_idx <- both_dt[["prop"]] >= min_proportion
  wanted_indexes_dt <- both_dt[kept_idx, ]
  kept_idx <- chng[["index"]] %in% wanted_indexes_dt[["index"]]
  chng <- chng[kept_idx, ]

  post_rows <- nrow(chng)
  post_ident <- nrow(ident)
  pct_rows <- scales::percent(x=post_rows / pre_rows, accuracy=0.01)
  pct_ident <- scales::percent(x=post_ident / pre_ident, accuracy=0.01)
  if (verbose) {
    message("All data: after proportion filtering, there are: ",
            post_rows, " changed reads: ", pct_rows, ".")
    message("All data: after proportion filtering, there are: ",
            post_ident, " identical reads: ", pct_ident, ".")
  }
  retlist <- list(
    "chng_start" = pre_rows,
    "chng_removed" = pre_rows - post_rows,
    "chng_remaining" = post_rows,
    "chng" = chng,
    "ident_start" = pre_ident,
    "ident_removed" = pre_ident - post_ident,
    "ident_remaining" = post_ident,
    "ident" = ident)
}

#' Filter away reads which are associated with indexes that have been deemed
#' insufficient.  These 'bad' indexes were defined by prune_indexes().
#'
#' At least in theory we should have at least x reads per index before that
#' index is deemed interesting.  prune_indexes() creates a list of indexes which
#' have more than this x reads per index and provides it to this.  This then
#' removes the rest.  Caveat: this filter works on both the mutated and
#' identical data.
#'
#' @param chng Data set of reads which are not identical to the template.
#' @param ident Data set of reads which are identical to the template.
#' @param wanted_indexes Character vector of indexes which have sufficient
#'  evidence.
#' @param min_reads Used when verbose to remind us of the number of
#'  reads/index.
#' @param verbose Print some information while running?
#' @return The remaining reads and a little summary information.
#' @export
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
    "chng_start" = pre_rows,
    "chng_removed" = pre_rows - post_rows,
    "chng_remaining" = post_rows,
    "chng" = chng,
    "ident_start" = pre_ident,
    "ident_removed" = pre_ident - post_ident,
    "ident_remaining" = post_ident,
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
#'  indexes kept.  The list of indexes kept may just be 'all'.
#' @export
prune_indexes <- function(chng, ident, min_reads=3, verbose=TRUE) {
  ## Get a matrix of how many reads/index are identical to the template.
  if (verbose) {
    message("Gathering information about the number of reads per index.")
  }

  ## Suck it r cmd check
  index <- NULL
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
