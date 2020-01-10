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
  normalized <- summary[["normalized"]]
  mc <- summary[["matrices_by_counts"]]
  nc <- summary[["normalized_by_counts"]]
  matrix_plots <- list()
  normal_plots <- list()
  matrix_count_plots <- list()
  normal_count_plots <- list()
  for (t in 1:length(names(matrices))) {
       tname <- names(matrices)[t]
       matrix <- as.matrix(matrices[[tname]])
       normal <- as.matrix(normalized[[tname]])
       m <- as.matrix(mc[[tname]])
       n <- as.matrix(nc[[tname]])
       matrix_melted <- reshape2::melt(data=matrix, value.name="norm")
       colnames(matrix_melted) <- c("position","sample", "number")
       normal_melted <- reshape2::melt(data=normal, value.name="norm")
       colnames(normal_melted) <- c("position","sample", "norm")
       mm <- reshape2::melt(data=m, value.name="norm")
       colnames(mm) <- c("position","sample", "norm")
       nm <- reshape2::melt(data=n, value.name="norm")
       colnames(nm) <- c("position","sample", "norm")

       matrix_plots[[tname]] <- ggplot2::ggplot(data=matrix_melted,
                                                mapping=ggplot2::aes_string(x="position", y="number", fill="sample")) +
         ggplot2::geom_col(position="dodge") +
         ggplot2::theme_bw()
       normal_plots[[tname]] <- ggplot2::ggplot(data=normal_melted,
                                                mapping=ggplot2::aes_string(x="position", y="norm", fill="sample")) +
         ggplot2::geom_col(position="dodge") +
         ggplot2::theme_bw()
       matrix_count_plots[[tname]] <- ggplot2::ggplot(data=mm,
                                                      mapping=ggplot2::aes_string(x="position", y="norm", fill="sample")) +
         ggplot2::geom_col(position="dodge") +
         ggplot2::theme_bw()
       normal_count_plots[[tname]] <- ggplot2::ggplot(data=nm,
                                                mapping=ggplot2::aes_string(x="position", y="norm", fill="sample")) +
         ggplot2::geom_col(position="dodge") +
         ggplot2::theme_bw()
  }
  retlist <- list(
    "matrices" = matrix_plots,
    "normal" = normal_plots,
    "matrix_count" = matrix_count_plots,
    "normal_count" = normal_count_plots
  )
  return(retlist)
}
