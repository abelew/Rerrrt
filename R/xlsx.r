#' Make a pretty xlsx file of the various summary statistics and data from an RT
#' error rate experiment.
#'
#' This takes the results of create_matrices() and hopefully prints some pretty
#' data about them.
#'
#' @param retlst Result from create_matrices()
#' @param excel Excel filename.
#' @param min_reads Written to a summary table.
#' @param min_indexes Written to a summary table.
#' @param min_sequencer Written to a summary table.
#' @param min_position Written to a summary table.
#' @param max_position Written to a summary table.
#' @param max_mutations_per_read Written to a summary table.
#' @param prune_n Written to a summary table.
#' @param verbose Print!
#' @export
write_matrices <- function(retlst, excel, min_reads=NULL, min_indexes=NULL,
                           min_sequencer=10, min_position=NULL,
                           max_position=NULL, max_mutations_per_read=NULL,
                           prune_n=TRUE, verbose=FALSE) {
  excel_basename <- gsub(pattern="\\.xlsx", replacement="", x=excel)
  excel_dir <- dirname(as.character(excel))
  if (!file.exists(excel_dir)) {
    dir.create(excel_dir, recursive=TRUE)
  }
  if (file.exists(excel) & verbose) {
    message("Deleting the file ", excel, " before writing the tables.")
    file.remove(excel)
  }
  wb <- openxlsx::createWorkbook(creator="Rerrrt")

  if (verbose) {
    message("  Writing a legend.")
  }
  legend_written <- write_legend_sheet(wb=wb, sheet="legend", min_reads=min_reads,
                                       min_indexes=min_indexes, min_sequencer=min_sequencer,
                                       min_position=min_position, max_position=max_position,
                                       max_mutations_per_read=max_mutations_per_read,
                                       prune_n=prune_n)
  sheet <- "summary"
  start_row <- 1
  start_col <- 1
  ## Write the metadata
  meta_written <- hpgltools::write_xlsx(wb=wb, data=retlst[["metadata"]], sheet=sheet,
                                        title="Sample Metadata.",
                                        start_row=start_row, start_col=start_col)
  plot_col <- meta_written[["end_col"]] + 6

  ## Write the set of remaining reads
  start_row <- meta_written[["end_row"]] + 2
  reads_written <- hpgltools::write_xlsx(wb=wb, data=retlst[["reads_remaining"]], sheet=sheet,
                                         title="Reads remaining after each filter.",
                                         start_row=start_row, start_col=start_col)

  ## Write the set of pruned reads
  start_row <- reads_written[["end_row"]] + 2
  filtered_written <- hpgltools::write_xlsx(wb=wb, data=retlst[["filtered"]], sheet=sheet,
                                            title="Reads removed after each filter.",
                                            start_row=start_row, start_col=start_col)

  ## Write the set of remaining indexes
  start_row <- filtered_written[["end_row"]] + 2
  indexes_written <- hpgltools::write_xlsx(wb=wb, data=retlst[["indexes_remaining"]], sheet=sheet,
                                           title="Indexes remaining after each filter.",
                                           start_row=start_row, start_col=start_col)

  ## Write out the total reads and indexes on a new set of rows.
  start_row <- indexes_written[["end_row"]] + 2
  reads_per_sample <- t(as.data.frame(retlst[["reads_per_sample"]]))
  rownames(reads_per_sample) <- "reads"
  reads_written <- hpgltools::write_xlsx(wb=wb, data=reads_per_sample, sheet=sheet,
                                         title="Total individual mutants and identical reads.",
                                         start_row=start_row, start_col=start_col)

  start_row <- reads_written[["end_row"]] + 2
  indexes <- t(as.data.frame(retlst[["indexes_per_sample"]]))
  rownames(indexes) <- "indexes"
  indexes_written <- hpgltools::write_xlsx(wb=wb, data=indexes, sheet=sheet,
                                           title="Final indexes per Sample.",
                                           start_row=start_row, start_col=start_col)

  ## Plot index distributions
  plot_row <- 1
  title <- "Index density for mutant reads before filtering."
  message("Plotting ", title)
  a_plot <- plot_index_density(data=retlst[["pre_chng_index_density_dt"]],
                               metadata=retlst[["metadata"]], title=title)
  plotted <- hpgltools::xlsx_plot_png(a_plot=a_plot, wb=wb, plotname="pre_chng_density",
                                      sheet=sheet, start_row=plot_row, start_col=plot_col)
  title <- "Index density for identical reads before filtering."
  message("Plotting ", title)
  a_plot <- plot_index_density(data=retlst[["pre_ident_index_density_dt"]],
                               metadata=retlst[["metadata"]], title=title)
  plotted <- hpgltools::xlsx_plot_png(a_plot=a_plot, wb=wb, plotname="pre_ident_density",
                                      sheet=sheet, start_row=plot_row, start_col=plot_col + 10)
  title <- "Index density for all reads before filtering."
  message("Plotting ", title)
  a_plot <- plot_index_density(data=retlst[["pre_index_density_dt"]],
                               metadata=retlst[["metadata"]], title=title)
  plotted <- hpgltools::xlsx_plot_png(a_plot=a_plot, wb=wb, plotname="pre_all_index_density",
                                      sheet=sheet, start_row=plot_row, start_col=plot_col + 20)
  plot_row <- plot_row + 31
  title <- "Index density for mutant reads after filtering."
  message("Plotting ", title)
  dt <- retlst[["post_chng_index_density_dt"]]
  a_plot <- plot_index_density(data=retlst[["post_chng_index_density_dt"]],
                               metadata=retlst[["metadata"]], title=title)
  plotted <- hpgltools::xlsx_plot_png(a_plot=a_plot, wb=wb, plotname="post_chng_density",
                                      sheet=sheet, start_row=plot_row, start_col=plot_col)
  title <- "Index density for identical reads after filtering."
  message("Plotting ", title)
  a_plot <- plot_index_density(data=retlst[["post_ident_index_density_dt"]],
                               metadata=retlst[["metadata"]], title=title)
  plotted <- hpgltools::xlsx_plot_png(a_plot=a_plot, wb=wb, plotname="post_ident_density",
                                      sheet=sheet, start_row=plot_row, start_col=plot_col + 10)
  title <- "Index density for all reads after filtering."
  message("Plotting ", title)
  a_plot <- plot_index_density(data=retlst[["post_index_density_dt"]],
                               metadata=retlst[["metadata"]], title=title)
  plotted <- hpgltools::xlsx_plot_png(a_plot=a_plot, wb=wb, plotname="post_all_index_density",
                                      sheet=sheet, start_row=plot_row, start_col=plot_col + 20)

  ## Print out raw matrices
  table_names <- list(
      "reads" = c("miss_reads_by_refnt", "miss_reads_by_hitnt", "miss_reads_by_type",
                  "miss_reads_by_trans", "miss_reads_by_strength", "miss_reads_by_position",
                  "insert_reads_by_nt", "delete_reads_by_nt", "miss_reads_by_position",
                  "insert_reads_by_position", "delete_reads_by_position",
                  "miss_reads_by_string"),
      "rt" = c("miss_indexes_by_refnt", "miss_indexes_by_hitnt", "miss_indexes_by_type",
               "miss_indexes_by_trans", "miss_indexes_by_strength", "miss_indexes_by_position",
               "insert_indexes_by_nt", "delete_indexes_by_nt", "miss_indexes_by_position",
               "insert_indexes_by_position", "delete_indexes_by_position",
               "miss_indexes_by_string"),
      "sequencer" = c("miss_sequencer_by_refnt", "miss_sequencer_by_hitnt",
                      "miss_sequencer_by_type", "miss_sequencer_by_trans",
                      "miss_sequencer_by_strength", "miss_sequencer_by_position",
                      "insert_sequencer_by_nt", "delete_sequencer_by_nt",
                      "miss_sequencer_by_position", "insert_sequencer_by_position",
                      "delete_sequencer_by_position", "miss_sequencer_by_string"))

  if (verbose) {
    message("  Writing raw data.")
  }
  raw_column <- write_matrix_sheet(retlst, wb, table_names,
                                   sheet="raw", table_type="matrices")
  if (verbose) {
    message("  Writing cpm data.")
  }
  cpm_column <-  write_matrix_sheet(retlst, wb, table_names,
                                    sheet="cpm", table_type="matrices_cpm")
  if (verbose) {
    message("  Writing data normalized by reads/indexes.")
  }
  counts_column <- write_matrix_sheet(retlst, wb, table_names,
                                      sheet="counts", table_type="matrices_counts")
  if (verbose) {
    message("  Writing data normalized by reads/indexes and length.")
  }
  countslength_column <- write_matrix_sheet(retlst, wb, table_names,
                                            sheet="countslength",
                                            table_type="matrices_countslength")

  if (verbose) {
    message("  Writing data normalized by cpm(reads/indexes) and length.")
  }
  countslength_column <- write_matrix_sheet(retlst, wb, table_names,
                                            sheet="cpmlength",
                                            table_type="matrices_cpmlength")

  finished <- openxlsx::saveWorkbook(wb=wb, file=excel, overwrite=TRUE)
  return(finished)
}

#' Write out a worksheet containing one group of error rate matrices
#'
#' This function is responsible for writing out the actual data.  Most of the
#' logic in it is in an attempt to properly space the various tables and plots
#' apart.
#'
#' @param retlst Set of results from create_matrices()
#' @param wb Workbook to write.
#' @param table_names The set of tables to write.
#' @param sheet Name of the sheet to write.
#' @param table_type I have a few different tables for each metric, raw, cpm,
#'  cpmlength, etc.  Choose one.
write_matrix_sheet <- function(retlst, wb, table_names, sheet="raw", table_type="matrices") {
  start_row <- 1
  read_col <- 1
  padding <- 3
  end_col <- 0
  first_plot_row <- 1
  second_plot_row <- 1
  first_plot_col <- ((nrow(retlst[["metadata"]]) + 3) * 3) + 2
  second_plot_col <- 1
  plot_width_per_inch <- 10 / 6
  plot_height <- 33
  read_tables <- table_names[["reads"]]
  rt_tables <- table_names[["rt"]]
  sequencer_tables <- table_names[["sequencer"]]
  for (t in 1:length(read_tables)) {
    read_name <- read_tables[t]
    read_table <- retlst[[table_type]][[read_name]]
    if (is.null(read_table)) {
      read_table <- data.frame()
    }
    rt_name <- rt_tables[t]
    rt_table <- retlst[[table_type]][[rt_name]]
    if (is.null(rt_table)) {
      rt_table <- data.frame()
    }
    sequencer_name <- sequencer_tables[t]
    sequencer_table <- retlst[[table_type]][[sequencer_name]]
    if (is.null(sequencer_table)) {
      sequencer_table <- data.frame()
    }
    read_rows <- nrow(read_table)
    rt_rows <- nrow(rt_table)
    sequencer_rows <- nrow(sequencer_table)
    max_rows <- max(read_rows, rt_rows, sequencer_rows)
    rt_col <- read_col + ncol(read_table) + padding
    sequencer_col <- rt_col + ncol(rt_table) + padding
    end_col <- sequencer_col + ncol(sequencer_table) + padding
    print_reads <- TRUE
    usable <- c("matrix", "data.frame", "data.table")
    if (sum(usable %in% class(read_table)) == 0) {
      print_reads <- FALSE
    } else if (nrow(read_table) == 0) {
      print_reads <- FALSE
    }

    if (isTRUE(print_reads)) {
      read_written <- hpgltools::write_xlsx(
                                   wb=wb, sheet=sheet,
                                   start_row=start_row, start_col=read_col,
                                   data=read_table, title=read_name)
    }

    print_reads <- TRUE
    if (sum(usable %in% class(read_table)) == 0) {
      print_reads <- FALSE
    } else if (nrow(read_table) == 0) {
      print_reads <- FALSE
    }
    if (isTRUE(print_reads)) {
      rt_written <- hpgltools::write_xlsx(
                                 wb=wb, sheet=sheet,
                                 start_row=start_row, start_col=rt_col,
                                 data=rt_table, title=rt_name)
      ##message("Adding plot ", rt_name)
      a_plot <- retlst[["plots"]][[table_type]][[rt_name]]
      x_length <- length(levels(as.factor(a_plot[["data"]][["category"]])))
      padding <- 3
      plot_width <- padding + (x_length / 4)
      rt_plot <- hpgltools::xlsx_plot_png(a_plot=a_plot, wb=wb, sheet=sheet,
                                          plotname=paste0(table_type, "_", rt_name),
                                          width=plot_width,
                                          start_col=first_plot_col,
                                          start_row=first_plot_row)
      plot_columns <- plot_width_per_inch * plot_width
      second_plot_col <- first_plot_col + plot_columns
      first_plot_row <- first_plot_row + plot_height
    }

    print_reads <- TRUE
    if (sum(usable %in% class(sequencer_table)) == 0) {
      print_reads <- FALSE
    } else if (nrow(sequencer_table) == 0) {
      print_reads <- FALSE
    }
    if (isTRUE(print_reads)) {
      seq_written <- hpgltools::write_xlsx(
                                    wb=wb, sheet=sheet,
                                    start_row=start_row,
                                    start_col=sequencer_col,
                                    data=sequencer_table,
                                    title=sequencer_name)
      ##message("Adding plot ", sequencer_name)
      a_plot <- retlst[["plots"]][[table_type]][[sequencer_name]]
      x_length <- length(levels(as.factor(a_plot[["data"]][["category"]])))
      first_plot_col <- end_col + 2
      padding <- 3
      plot_width <- padding + (x_length / 4)
      rt_plot <- hpgltools::xlsx_plot_png(a_plot=a_plot, wb=wb, sheet=sheet,
                                          plotname=paste0(table_type, "_", sequencer_name),
                                          width=plot_width,
                                          start_col=second_plot_col,
                                          start_row=second_plot_row)
      second_plot_row <- second_plot_row + plot_height
    }
    if (max_rows > 0) {
      start_row <- start_row + max(read_rows, rt_rows, sequencer_rows) + padding
    }

  } ## End of the for loop
  return(end_col)
}

#' Write an initial sheet describing what to expect in the rest of the workbook.
#'
#' @param wb Workbook to write
#' @param sheet Worksheet to write
#' @param min_reads Passed down from create_matrices, printed in a table.
#' @param min_indexes Passed down from create_matrices, printed in a table.
#' @param min_sequencer Passed down from create_matrices, printed in a table.
#' @param min_position Passed down from create_matrices, printed in a table.
#' @param max_position Passed down from create_matrices, printed in a table.
#' @param max_mutations_per_read Passed down from create_matrices, printed in a table.
#' @param prune_n Passed down from create_matrices, printed in a table.
write_legend_sheet <- function(wb, sheet="legend", min_reads=NULL, min_indexes=NULL,
                               min_sequencer=10, min_position=NULL, max_position=NULL,
                               max_mutations_per_read=NULL, prune_n=NULL) {
  legend <- data.frame(rbind(
    c("Sheet Name", "Table", "Purpose"),
    c("Summary", "Sample Metadata", "Describe the samples used in the analysis."),
    c("Summary", "Reads remaining at each step", "How many reads remain after each filter."),
    c("Summary", "Reads filtered at each step", "How many reads were removed after each filter."),
    c("Summary", "Indexes remaining at each step", "How many indexes remain after each filter."),
    c("Summary", "Reads per sample", "Raw reads counted in each sample."),
    c("Summary", "Indexes per sample", "Final index count after all filters."),
    c("Raw", "miss_reads_by_refnt", "Mismatch reads observed, counted by template nucleotide."),
    c("Raw", "miss_indexes_by_refnt", "Mismatch indexes observed, counted by template nucleotide."),
    c("Raw", "miss_sequencer_by_refnt", "Mismatches indexes from the sequencer, counted by template nt."),
    c("Raw", "miss_reads_by_hitnt", "Mismatch reads observed, counted by the new nucleotide."),
    c("Raw", "miss_indexes_by_hitnt", "Mismatch indexes observed, counted by the new nucleotide."),
    c("Raw", "miss_sequencer_by_hitnt", "Mismatches indexes from the sequencer, counted by new nt."),
    c("Raw", "miss_reads_by_type", "Mismatch reads observed, counted by x->y."),
    c("Raw", "miss_indexes_by_hitnt", "Mismatch indexes observed, counted by x->y."),
    c("Raw", "miss_sequencer_by_hitnt", "Mismatches indexes from the sequencer, counted by x->y."),
    c("Raw", "miss_reads_by_type", "Mismatch reads observed, counted by x->y."),
    c("Raw", "miss_indexes_by_hitnt", "Mismatch indexes observed, counted by x->y."),
    c("Raw", "miss_sequencer_by_hitnt", "Mismatches indexes from the sequencer, counted by x->y."),
    c("Raw", "miss_reads_by_trans", "Mismatch reads observed, counted by transition/transversion."),
    c("Raw", "miss_indexes_by_trans", "Mismatch indexes observed, counted by transition/transversion."),
    c("Raw", "miss_sequencer_by_hitnt", "Mismatches indexes from the sequencer, counted by transition/transversion."),
    c("Raw", "miss_reads_by_strength", "Mismatch reads observed, counted by nucleotide H2 bond strength."),
    c("Raw", "miss_indexes_by_trans", "Mismatch indexes observed, counted by strength."),
    c("Raw", "miss_sequencer_by_hitnt", "Mismatches indexes from the sequencer, counted by strength."),
    c("Raw", "miss_reads_by_position", "Mismatch reads observed, counted by template position."),
    c("Raw", "miss_indexes_by_position", "Mismatch indexes observed, counted by template position."),
    c("Raw", "miss_sequencer_by_position", "Mismatches indexes from the sequencer, counted position."),
    c("Raw", "insert_reads_by_nt", "Insertion reads observed, counted by inserted nucleotide."),
    c("Raw", "insert_indexes_by_nt", "Insertion indexes observed, counted by inserted nucleotide."),
    c("Raw", "insert_sequencer_by_nt", "Insertion indexes from the sequencer, counted by insertion."),
    c("Raw", "delete_reads_by_nt", "Deletion reads observed, counted by deleted nucleotide."),
    c("Raw", "delete_indexes_by_nt", "Deletion indexes observed, counted by deleted nucleotide."),
    c("Raw", "delete_sequencer_by_nt", "Deletion indexes from the sequencer, counted by nt."),
    c("Raw", "insert_reads_by_position", "Insertion reads observed, counted by inserted position."),
    c("Raw", "insert_indexes_by_position", "Insertion indexes observed, counted by inserted position."),
    c("Raw", "insert_sequencer_by_positino", "Insertion indexes from the sequencer, counted by position."),
    c("Raw", "delete_reads_by_position", "Deletion reads observed, counted by position."),
    c("Raw", "delete_indexes_by_position", "Deletion indexes observed, counted by position."),
    c("Raw", "delete_sequencer_by_positino", "Deletion indexes from the sequencer, counted by position."),
    c("Raw", "miss_reads_by_string", "Reads counted by all categories of mismatch in one string."),
    c("Raw", "miss_indexes_by_string", "Mismatch indexes observed, counted by string."),
    c("Raw", "delete_sequencer_by_positino", "Mismatch indexes from the sequencer, counted by string."),
    c("cpm", "Same tables as 'raw'.", "The tables from 'raw' rewritten as counts per million."),
    c("counts", "Same tables as 'raw'.", "The tables from 'raw' divided by the number of reads or indexes in each sample."),
    c("countslength", "Same tables as 'raw'.", "The tables from 'counts' divided by the number of nucleotides queried.")),
    stringsAsFactors=FALSE)
  xlsx_result <- hpgltools::write_xlsx(
                              wb, data=legend, sheet=sheet,
                              title="Summary of the sheets and tables used in this workbook.")

  params <- data.frame(
    rbind(
      c("Parameter", "Purpose", "Setting"),
      c("min_reads", "Minimum number of reads for an index to be deemed 'real'", min_reads),
      c("min_indexes", "Minimum number of index groups for a mutation to be deemed 'real'", min_indexes),
      c("min_sequencer", "For a group's mutation to be deemed from the sequencer, it must be 1 read out of min_sequencer for an index group", min_sequencer),
      c("min_position", "Filter out errors observed before this position", min_position),
      c("max_position", "Filter out errors observed after this position", max_position),
      c("max_mutations_per_read", "Filter out reads with more than this number of mutations", max_mutations_per_read),
      c("prune_n", "Filter out mutations from Ns in the read", prune_n)))
  second_xlsx_result <- hpgltools::write_xlsx(
                                     wb, data=params, sheet=sheet, start_row=2, start_col=7,
                                     title="Parameters used to generate this workbook.")
  return(xlsx_result)
}
