#' Make a pretty xlsx file of the various summary statistics and data from an RT
#' error rate experiment.
#'
#' This takes the results of create_matrices() and hopefully prints some pretty
#' data about them.
#'
#' @param retlist Result from create_matrices()
#' @param excel Excel filename.
#' @export
write_matrices <- function(retlist, excel) {
  excel_basename <- gsub(pattern="\\.xlsx", replacement="", x=excel)
  excel_dir <- dirname(as.character(excel))
  if (!file.exists(excel_dir)) {
    dir.create(excel_dir, recursive=TRUE)
  }
  if (file.exists(excel)) {
    message("Deleting the file ", excel, " before writing the tables.")
    file.remove(excel)
  }
  wb <- openxlsx::createWorkbook(creator="Rerrrt")

  sheet <- "summary"
  start_row <- 1
  start_col <- 1
  ## Write the metadata
  meta_written <- hpgltools::write_xlsx(wb=wb, data=retlist[["metadata"]], sheet=sheet,
                                        title="Sample Metadata.",
                                        start_row=start_row, start_col=start_col)
  plot_column <- meta_written[["end_col"]] + 6

  ## Write the set of remaining reads
  start_row <- meta_written[["end_row"]] + 2
  reads_written <- hpgltools::write_xlsx(wb=wb, data=retlist[["reads_remaining"]], sheet=sheet,
                                         title="Reads remaining at each step.",
                                         start_row=start_row, start_col=start_col)

  ## Write the set of pruned reads
  start_row <- reads_written[["end_row"]] + 2
  filtered_written <- hpgltools::write_xlsx(wb=wb, data=retlist[["filtered"]], sheet=sheet,
                                            title="Reads Filtered at Each Step.",
                                            start_row=start_row, start_col=start_col)

  ## Write the set of remaining indexes
  start_row <- filtered_written[["end_row"]] + 2
  indexes_written <- hpgltools::write_xlsx(wb=wb, data=retlist[["indexes_remaining"]], sheet=sheet,
                                           title="Reads remaining at each step.",
                                           start_row=start_row, start_col=start_col)

  ## Write out the total reads and indexes on a new set of rows.
  start_row <- indexes_written[["end_row"]] + 2
  reads_per_sample <- t(as.data.frame(retlist[["reads_per_sample"]]))
  rownames(reads_per_sample) <- "reads"
  reads_written <- hpgltools::write_xlsx(wb=wb, data=reads_per_sample, sheet=sheet,
                                         title="Reads per Sample.",
                                         start_row=start_row, start_col=start_col)

  start_row <- reads_written[["end_row"]] + 2
  indexes <- t(as.data.frame(retlist[["indexes_per_sample"]]))
  rownames(indexes) <- "indexes"
  indexes_written <- hpgltools::write_xlsx(wb=wb, data=indexes, sheet=sheet,
                                           title="Indexes per Sample.",
                                           start_row=start_row, start_col=start_col)

  ## Plot index distributions
  plot_row <- 1
  plotted <- hpgltools::xlsx_plot_png(a_plot=retlist[["pre_index_density_plot"]],
                                      sheet=sheet, start_row=plot_row, start_col=plot_col)
  plot_row <- plot_row + 20
  plotted <- hpgltools::xlsx_plot_png(a_plot=retlist[["post_index_density_plot"]],
                                      sheet=sheet, start_row=plot_row, start_col=plot_col)

  ## Print out raw matrices
  read_tables <- c("miss_reads_by_refnt", "miss_reads_by_hitnt", "miss_reads_by_type",
                   "miss_reads_by_trans", "miss_reads_by_strength", "miss_reads_by_position",
                   "insert_reads_by_nt", "delete_reads_by_nt", "miss_reads_by_position",
                   "insert_reads_by_position", "delete_reads_by_position",
                   "miss_reads_by_string")
  rt_tables <- c("miss_indexes_by_refnt", "miss_indexes_by_hitnt", "miss_indexes_by_type",
                 "miss_indexes_by_trans", "miss_indexes_by_strength", "miss_indexes_by_position",
                 "insert_indexes_by_nt", "delete_indexes_by_nt", "miss_indexes_by_position",
                 "insert_indexes_by_position", "delete_indexes_by_position",
                 "miss_indexes_by_string")
  sequencer_tables <- c("miss_sequencer_by_refnt", "miss_sequencer_by_hitnt", "miss_sequencer_by_type",
                        "miss_sequencer_by_trans", "miss_sequencer_by_strength", "miss_sequencer_by_position",
                        "insert_sequencer_by_nt", "delete_sequencer_by_nt", "miss_sequencer_by_position",
                        "insert_sequencer_by_position", "delete_sequencer_by_position",
                        "miss_sequencer_by_string")

  raw_column <- write_matrix_sheet(retlist, wb, "raw", "matrices",
                                   read_tables, rt_tables, sequencer_tables)
  cpm_column <-  write_matrix_sheet(retlist, wb, "cpm", "matrices_cpm",
                                    read_tables, rt_tables, sequencer_tables)
  cpmlength_column <- write_matrix_sheet(retlist, wb, "cpmlength", "matrices_cpmlength",
                                         read_tables, rt_tables, sequencer_tables)
  counts_column <- write_matrix_sheet(retlist, wb, "counts", "matrices_counts",
                                      read_tables, rt_tables, sequencer_tables)
  countslength_column <- write_matrix_sheet(retlist, wb, "countslength", "matrices_countslength",
                                            read_tables, rt_tables, sequencer_tables)

  finished <- openxlsx::saveWorkbook(wb=wb, file=excel, overwrite=TRUE)
  return(finished)
}


write_matrix_sheet <- function(retlist, wb, sheet, type, read_tables,
                               rt_tables, sequencer_tables) {
  start_row <- 1
  read_col <- 1
  padding <- 3
  end_col <- 0
  for (t in 1:length(read_tables)) {
    read_name <- read_tables[t]
    read_table <- retlist[[type]][[read_name]]
    rt_name <- rt_tables[t]
    rt_table <- retlist[[type]][[rt_name]]
    sequencer_name <- sequencer_tables[t]
    sequencer_table <- retlist[[type]][[sequencer_name]]
    read_rows <- nrow(read_table)
    rt_rows <- nrow(rt_table)
    sequencer_rows <- nrow(sequencer_table)
    max_rows <- max(read_rows, rt_rows, sequencer_rows)
    rt_col <- read_col + ncol(read_table) + padding
    sequencer_col <- rt_col + ncol(rt_table) + padding
    end_col <- sequencer_col + ncol(sequencer_table) + padding
    if (nrow(read_table) > 0) {
      read_written <- hpgltools::write_xlsx(
                                   wb=wb, sheet=sheet,
                                   start_row=start_row, start_col=read_col,
                                   data=read_table,
                                   title=read_name)
    }
    if (nrow(rt_table) > 0) {
      rt_written <- hpgltools::write_xlsx(
                                 wb=wb, sheet=sheet,
                                 start_row=start_row, start_col=rt_col,
                                 data=rt_table,
                                 title=rt_name)
    }
    if (nrow(sequencer_table) > 0) {
      seq_written <- hpgltools::write_xlsx(
                                  wb=wb, sheet=sheet,
                                  start_row=start_row, start_col=sequencer_col,
                                  data=sequencer_table,
                                  title=sequencer_name)
    }
    if (max_rows > 0) {
      start_row <- start_row + max(read_rows, rt_rows, sequencer_rows) + padding
    }
  }
  return(end_col)
}
