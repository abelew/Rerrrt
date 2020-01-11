#' Make a bar plot for every type of data returned by create_matrices.
#'
#' The function create_matrices() creates two lists containing the various
#' categorizations of the data.  This takes that information and makes a bar
#' plot for every element in those lists.
#'
#' @param summary Result from create_matrices()
#' @return List full of bar plots!
#' @export
barplot_matrices <- function(summary) {
  matrices <- summary[["matrices"]]
  matrices_cpm <- summary[["matrices_cpm"]]
  matrices_cpmlength <- summary[["matrices_cpmlength"]]
  matrices_counts <- summary[["matrices_counts"]]
  matrices_countslength <- summary[["matrices_countslength"]]
  matrix_plots <- list()
  cpm_plots <- list()
  cpmlength_plots <- list()
  counts_plots <- list()
  countslength_plots <- list()
  for (t in 1:length(names(matrices))) {
       tname <- names(matrices)[t]
       matrix <- as.matrix(matrices[[tname]])
       matrix_melted <- reshape2::melt(data=matrix, value.name="norm")
       colnames(matrix_melted) <- c("position","sample", "number")
       matrix_plots[[tname]] <- ggplot(
         data=matrix_melted,
         mapping=ggplot2::aes_string(x="position", y="number", fill="sample")) +
         ggplot2::geom_col(position="dodge") +
         ggplot2::theme_bw()

       mcpm <- as.matrix(matrices_cpm[[tname]])
       mcpm_melted <- reshape2::melt(data=mcpm, value.name="norm")
       colnames(mcpm_melted) <- c("position","sample", "number")
       cpm_plots[[tname]] <- ggplot(
         data=mcpm_melted,
         mapping=ggplot2::aes_string(x="position", y="number", fill="sample")) +
         ggplot2::geom_col(position="dodge") +
         ggplot2::theme_bw()

       mcpml <- as.matrix(matrices_cpmlength[[tname]])
       mcpml_melted <- reshape2::melt(data=mcpml, value.name="norm")
       colnames(mcpml_melted) <- c("position","sample", "number")
       cpmlength_plots[[tname]] <- ggplot(
         data=mcpml_melted,
         mapping=ggplot2::aes_string(x="position", y="number", fill="sample")) +
         ggplot2::geom_col(position="dodge") +
         ggplot2::theme_bw()

       mcount <- as.matrix(matrices_counts[[tname]])
       mcount_melted <- reshape2::melt(data=mcount, value.name="norm")
       colnames(mcount_melted) <- c("position","sample", "number")
       counts_plots[[tname]] <- ggplot(
         data=mcount_melted,
         mapping=ggplot2::aes_string(x="position", y="number", fill="sample")) +
         ggplot2::geom_col(position="dodge") +
         ggplot2::theme_bw()

       mcountl <- as.matrix(matrices_countslength[[tname]])
       mcountl_melted <- reshape2::melt(data=mcountl, value.name="norm")
       colnames(mcountl_melted) <- c("position","sample", "number")
       countslength_plots[[tname]] <- ggplot(
         data=mcountl_melted,
         mapping=ggplot2::aes_string(x="position", y="number", fill="sample")) +
         ggplot2::geom_col(position="dodge") +
         ggplot2::theme_bw()
  }
  retlist <- list(
    "matrices" = matrix_plots,
    "cpm" = cpm_plots,
    "cpmlength" = cpmlength_plots,
    "counts" = counts_plots,
    "countslength" = countslength_plots)
  return(retlist)
}

#' Make a density plot of how many indexes are associated with each sample.
#'
#' @param df Table of indexes/sample.
#' @param max Maximum number of indexes to plot.
#' @return Density plot!
#' @export
plot_index_density <- function(df, max=20) {
  df[["sample"]] <- as.factor(df[["sample"]])
  df_count <- df %>%
    group_by(sample, index) %>%
    summarise(index_count=n())
  plt <- ggplot2::ggplot(data=df_count,
                         mapping=ggplot2::aes_string(x="index_count",
                                                     alpha=0.4,
                                                     fill="sample")) +
    ggplot2::geom_density(adjust=4) +
    ggplot2::xlim(c(1, max))
  return(plt)
}
