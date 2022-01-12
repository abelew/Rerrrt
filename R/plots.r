#' Make a bar plot for every type of data returned by create_matrices.
#'
#' The function create_matrices() creates two lists containing the various
#' categorizations of the data.  This takes that information and makes a bar
#' plot for every element in those lists.
#'
#' @param retlist Result from create_matrices() with all the fun stuff to plot.
#' @param bar_title Title for the plot.
barplot_matrices <- function(retlist, bar_title = "") {
  meta <- retlist[["metadata"]]
  matrices <- retlist[["matrices"]]
  matrices_cpm <- retlist[["matrices_cpm"]]
  matrices_cpmlength <- retlist[["matrices_cpmlength"]]
  matrices_counts <- retlist[["matrices_counts"]]
  matrices_countslength <- retlist[["matrices_countslength"]]
  matrix_plots <- list()
  cpm_plots <- list()
  cpmlength_plots <- list()
  counts_plots <- list()
  countslength_plots <- list()
  for (t in 1:length(names(matrices))) {
    tname <- names(matrices)[t]
    matrix_plots[[tname]] <- categorical_barplot(matrices[[tname]], meta, bar_title = tname)
    cpm_plots[[tname]] <- categorical_barplot(matrices_cpm[[tname]], meta, bar_title = tname)
    cpmlength_plots[[tname]] <- categorical_barplot(matrices_cpmlength[[tname]], meta, bar_title = tname)
    counts_plots[[tname]] <- categorical_barplot(matrices_counts[[tname]], meta, bar_title = tname)
    countslength_plots[[tname]] <- categorical_barplot(matrices_countslength[[tname]],
                                                       meta, bar_title = tname)
  }
  retlist <- list(
    "matrices" = matrix_plots,
    "matrices_cpm" = cpm_plots,
    "matrices_cpmlength" = cpmlength_plots,
    "matrices_counts" = counts_plots,
    "matrices_countslength" = countslength_plots)
  return(retlist)
}

categorical_barplot <- function(mtrx, meta, bar_title = "") {
  if (! "matrix" %in% class(mtrx) & ! "data.frame" %in% class(mtrx)) {
    return(NULL)
  }
  if (length(mtrx) == 0) {
    return(NULL)
  } else if (is.null(mtrx)) {
    return(NULL)
  } else if (nrow(mtrx) < 4) {
    return(NULL)
  }
  mtrx_df <- as.data.frame(mtrx)
  mtrx_df[["names"]] <- as.factor(mtrx_df[["names"]])
  matrix_melted <- reshape2::melt(data = mtrx_df, value.name = "norm", id.vars = "names")
  colnames(matrix_melted) <- c("category","sampleid", "number")
  matrix_melted <- merge(meta, matrix_melted, by = "sampleid")
  max_chars <- max(nchar(as.character(matrix_melted[["category"]])))
  angle <- 0
  if (max_chars > 3) {
    angle <- 90
  }

  plt <- ggplot(data = matrix_melted,
                mapping = ggplot2::aes_string(x = "category", y = "number", fill = "sampletype")) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::ggtitle(bar_title) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = angle, hjust = 1))
  return(plt)
}

#' Make a density plot of how many tags are associated with each sample.
#'
#' @param data Table of tags/sample.
#' @param metadata Experiment metadata, used to acquire colors.
#' @param title Plot title.
#' @return Density plot!
#' @export
plot_tag_density <- function(data, metadata, title = "") {
  rnames <- as.numeric(data[["num"]])
  dt <- as.matrix(data[, -1])
  na_idx <- is.na(dt)
  dt[na_idx] <- 0
  rownames(dt) <- rnames
  reordered <- order(rnames, decreasing = FALSE)
  dt <- dt[reordered, ]
  dt <- log2(dt + 1)
  dt <- make_frequency_matrix(dt)
  long <- reshape2::melt(dt, id.vars = "num", value.name = "num")
  colnames(long) <- c("tag_count", "sampleid", "frequency")
  long[["tag_count"]] <- log2(long[["tag_count"]])

  ## Merge the metadata and sample data.
  long <- merge(metadata, long, by = "sampleid")
  long[["violinwidth"]] <- long[["frequency"]]
  plt <- ggplot2::ggplot(data = long,
                         mapping = aes_string(x = "sampleid", y = "tag_count",
                                            fill = "sampletype")) +
    ggplot2::ggtitle(label = title) +
    ggplot2::ylab(label = "log2(Reads per tag)") +
    suppressWarnings(ggplot2::geom_violin(stat = "identity",
                                          mapping = aes_string(violinwidth = "frequency"))) +  ## Something is screwy with this.
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  return(plt)
}

tag_histogram <- function(sample_sheet = NULL, melted = NULL, min = 1, max = 10) {
  hist <- NULL
  melted <- NULL
  if (is.null(sample_sheet) & is.null(melted)) {
    stop("This requires either a sample sheet from which to acquire data or a melted data frame.")
  } else if (is.null(sample_sheet)) {
    ## We already have the melted data, so plot it.
    hist <- plot_tag_histogram(melted, min = min, max = max)
  } else {
    melted <- gather_reads_per_tag(sample_sheet)
    hist <- plot_tag_histogram(melted, min = min, max = max)
  }
  retlist <- list(
      "histogram" = hist,
      "melted" = melted)
  return(retlist)
}

gather_reads_per_tag <- function(sample_sheet) {
  meta <- extract_metadata(sample_sheet)
  hist_files <- paste0("preprocessing/", rownames(meta), "/idx_count.txt.xz")
  hist_df <- data.frame()
  for (f in 1:length(hist_files)) {
    file <- hist_files[f]
    name <- rownames(meta)[f]
    message("Starting ", name, ".")
    datum <- read.table(file)
    hist_data <- table(datum[["V2"]])
    start <- as.data.frame(hist_data)
    rownames(start) <- start[["Var1"]]
    colnames(start) <- c("num", name)
    if (f == 1) {
      hist_df <- start
    } else {
      hist_df <- merge(hist_df, start, by = "num", all = TRUE)
    }
  }
  na_idx <- is.na(hist_df)
  hist_df[na_idx] <- 0

  hist_melted <- reshape2::melt(hist_df)
  return(hist_melted)
}

plot_tag_histogram <- function(melted, min = 1, max = 10) {
  keep_idx <- as.numeric(melted[["num"]]) <= max
  melted <- melted[keep_idx, ]
  keep_idx <- as.numeric(melted[["num"]]) >= min
  melted <- melted[keep_idx, ]
  reordered_idx <- order(as.numeric(melted[["num"]]))
  melted <- melted[reordered_idx, ]
  new_levels <- levels(as.factor(melted[["num"]]))
  levels_idx <- order(as.numeric(new_levels))
  new_levels <- new_levels[levels_idx]
  melted[["num"]] <- factor(melted[["num"]], levels = new_levels)
  hist_plot <- ggplot(melted, aes_string(x = "num", y = "value",
                                         color = "variable", fill = "variable")) +
    suppressWarnings(geom_col(stat = "identity", position = "dodge")) +
    ggtitle("Number of tags with x reads/tag from 1-10")
  return(hist_plot)
}

