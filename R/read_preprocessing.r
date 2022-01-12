read_preprocessing <- function(meta, preprocess_column) {
  input_files <- meta[[preprocess_column]]
  hist_df <- data.frame()
  tag_lst <- list()
  if (is.null(input_files)) {
    return(NULL)
  }
  for (f in 1:length(input_files)) {
    file <- input_files[f]
    name <- rownames(meta)[f]
    datum <- read.table(file)
    colnames(datum) <- c("tag", "num")
    tag_lst[[name]] <- datum
    hist_data <- table(datum[["num"]])
    start <- as.data.frame(hist_data)
    rownames(start) <- start[["Var1"]]
    colnames(start) <- c("num", name)
    if (f == 1) {
      hist_df <- start
    } else {
      hist_df <- merge(hist_df, start, by = "num", all = TRUE)
    }
  } ## End of the for loop
  na_idx <- is.na(hist_df)
  hist_df[na_idx] <- 0
  retlist <- list(
      "all_tags" = tag_lst,
      "tag_counts" = hist_df)
  return(retlist)
}
