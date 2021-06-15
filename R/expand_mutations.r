#' Combine the columns describing a mutation into a single column and categorize.
#'
#' Given a table with columns including the position, mutation type, reference
#' nt at that position, and product nt at that position; create a single column
#' from it, standardize it, and make some categories describing each mutation.
#' These columns currently include:  'mt': X_Y telling that this is a mutation
#' from X to Y, 'transition_transversion': this is a transition or transversion
#' mismatch, 'strong_weak': this mismatch went from strong->weak, weak->strong,
#' weak->weak, or strong->strong.  undef is used in the case of indels.
#'
#' @param t Table from readr.
#' @return The same table with some new columns describing the mutations therein.
expand_mutation_string <- function(t) {
  t[["string"]] <- paste0(t[["position"]], "_",
                          t[["type"]], "_",
                          t[["reference"]], "_",
                          t[["hit"]])
  t[["string"]] <- gsub(x=t[["string"]],
                        pattern="^([[:digit:]])_",
                        replacement="00\\1_")
  t[["string"]] <- gsub(x=t[["string"]],
                        pattern="^([[:digit:]][[:digit:]])_",
                        replacement="0\\1_")
  t[["string"]] <- as.factor(t[["string"]])
  t[["mt"]] <- paste0(t[["reference"]], "_", t[["hit"]])
  t[["transition_transversion"]] <- "undef"
  transition_idx <- t[["mt"]] == "A_G" | t[["mt"]] == "G_A" |
    t[["mt"]] == "T_C" | t[["mt"]] == "C_T"
  t[transition_idx, "transition_transversion"] <- "transition"
  transversion_idx <- t[["mt"]] == "A_T" | t[["mt"]] == "A_C" |
    t[["mt"]] == "T_A" | t[["mt"]] == "C_A" |
    t[["mt"]] == "G_C" | t[["mt"]] == "G_T" |
    t[["mt"]] == "C_G" | t[["mt"]] == "T_G"
  t[transversion_idx, "transition_transversion"] <- "transversion"
  t[["transition_transversion"]] <- as.factor(t[["transition_transversion"]])
  t[["strong_weak"]] <- "undef"
  weak_strong_idx <- t[["mt"]] == "A_G" | t[["mt"]] == "T_G" |
    t[["mt"]] == "A_C" | t[["mt"]] == "T_C"
  t[weak_strong_idx, "strong_weak"] <- "weak_strong"
  strong_weak_idx <- t[["mt"]] == "G_A" | t[["mt"]] == "G_T" |
    t[["mt"]] == "C_A" | t[["mt"]] == "C_T"
  t[strong_weak_idx, "strong_weak"] <- "strong_weak"
  weak_weak_idx <- t[["mt"]] == "A_T" | t[["mt"]] == "T_A"
  t[weak_weak_idx, "strong_weak"] <- "weak_weak"
  strong_strong_idx <- t[["mt"]] == "G_C" | t[["mt"]] == "C_G"
  t[strong_strong_idx, "strong_weak"] <- "strong_strong"
  t[["strong_weak"]] <- as.factor(t[["strong_weak"]])
  ## t <- get_previous_nts(t)
  ## t[["previous_nt"]] <- as.factor(t[["previous_nt"]])
  ## t[["previous_strength"]] <- as.factor(t[["previous_strength"]])
  ## t <- get_next_nts(t)
  ## t[["next_nt"]] <- as.factor(t[["next_nt"]])
  ## t[["next_strength"]] <- as.factor(t[["next_strength"]])
  return(t)
}

make_template_df <- function(t) {
  subset <- unique(t[, c("position", "reference")])
  ## FIXME Why are there NAs in this?
  keep_idx <- is.na(subset[["reference"]])
  subset <- subset[!keep_idx, ]
  subset_idx <- order(subset[["position"]], decreasing=FALSE)
  subset <- subset[subset_idx, ]
  subset[["position"]] <- as.numeric(subset[["position"]])
  rownames(subset) <- subset[["position"]]
  start <- min(subset[["position"]])
  end <- max(subset[["position"]])

  ## I don't think the following for loop is necessary, but I'm playing defensively.
  for (pos in start:end) {
    check <- is.na(subset[as.character(pos), "position"])
    if (isTRUE(check)) {
      subset[pos, "position"] <- position
      subset[pos, "reference"] <- "N"
    }
  }
  missing_idx <- subset[["reference"]] == "N"
  if (sum(missing_idx) > 0) {
    message("There are ", sum(missing_idx), " missing positions in this data.")
  }
  retlist <- list(
    "start" = start,
    "end" = end,
    "template_df" = subset)
  return(retlist)
}

get_previous_nts <- function(t) {
  t[["previous_nt"]] <- "N"
  t[["previous_strength"]] <- "template"
  template <- make_template_df(t)
  start <- template[["start"]]
  end <- template[["end"]]
  subset <- template[["template_df"]]
  for (pos in start:end) {
    previous_idx <- t[["position"]] == pos
    ## Add logic to handle the first nucleotide
    if (pos == start) {
      t[previous_idx, "previous_nt"] <- "fiveprime"
      t[previous_idx, "previous_strength"] <- "template"
    } else {
      template_idx <- as.character(pos - 1)
      previous_nt <- subset[template_idx, "reference"]
      previous_strength <- "weak"
      if (previous_nt == "G" | previous_nt == "C") {
        previous_strength <- "strong"
      }
      t[previous_idx, "previous_nt"] <- previous_nt
      t[previous_idx, "previous_strength"] <- previous_strength
    }
  }
  return(t)
}

get_next_nts <- function(t) {
  t[["next_nt"]] <- "N"
  t[["next_strength"]] <- "template"
  template <- make_template_df(t)
  start <- template[["start"]]
  end <- template[["end"]]
  subset <- template[["template_df"]]
  for (pos in start:end) {
    next_idx <- t[["position"]] == pos
    ## Add logic to handle the first nucleotide
    if (pos == end) {
      t[next_idx, "next_nt"] <- "threeprime"
      t[next_idx, "next_strength"] <- "template"
    } else {
      template_idx <- as.character(pos + 1)
      next_nt <- subset[template_idx, "reference"]
      next_strength <- "weak"
      if (next_nt == "G" | next_nt == "C") {
        next_strength <- "strong"
      }
      t[next_idx, "next_nt"] <- next_nt
      t[next_idx, "next_strength"] <- next_strength
    }
  }
  return(t)
}
