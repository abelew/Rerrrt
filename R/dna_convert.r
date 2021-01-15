dna_to_bit <- function(start) {
  lookup_table <- list(
      "A" = "00", "T" = "01", "G" = "10", "C" = "11")
}

dna_to_2bit <- function(start) {
  ## If there is an odd number of characters, prefix it with A.
  evenp <- nchar(start) %% 2
  ## I am not sure why, but I am getting:
  ## "The condition has length > 1 and only the first element will be used."
  ## I do not see where this will get > 1 element.
  if (evenp[1] == 1) {
    start <- paste0("A", start)
  }
  split <- strsplit(x=start, split="")[[1]]
  pairs <- paste0(split[c(TRUE, FALSE)], split[c(FALSE, TRUE)])
  lookup_table <- list(
      "AA" = "0", "AT" = "1", "AG" = "2", "AC" = "3",
      "TA" = "4", "TT" = "5", "TG" = "6", "TC" = "7",
      "GA" = "8", "GT" = "9", "GG" = "a", "GC" = "b",
      "CA" = "c", "CT" = "d", "CG" = "e", "CC" = "f")
  pieces <- ""
  for (c in 1:length(pairs)) {
    pieces <- paste0(pieces, lookup_table[[pairs[c]]])
  }
  return(pieces)
}

dna_to_6bit <- function(start) {
  ## If there is an odd number of characters, prefix it with A.
  char_mod <- nchar(start) %% 3
  if (char_mod[1] == 2) {
    start <- paste0("A", start)
  } else  if (char_mod[1] == 1) {
    start <- paste0("AA", start)
  }

  split <- strsplit(x=start, split="")[[1]]
  triplets <- c()
  while (length(split > 0)) {
    triplets <- c(triplets, paste0(split[1], split[2], split[3]))
    split <- split[-1]
    split <- split[-1]
    split <- split[-1]
  }
  lookup_table <- list(
      "AAA" = "0", "AAT" = "1", "AAG" = "2", "AAC" = "3",
      "ATA" = "4", "ATT" = "5", "ATG" = "6", "ATC" = "7",
      "AGA" = "8", "AGT" = "9", "AGG" = ":", "AGC" = ";",
      "ACA" = "<", "ACT" = "=", "ACG" = ">", "ACC" = "?",

      "TAA" = "@", "TAT" = "A", "TAG" = "B", "TAC" = "C",
      "TTA" = "D", "TTT" = "E", "TTG" = "F", "TTC" = "G",
      "TGA" = "H", "TGT" = "I", "TGG" = "J", "TGC" = "K",
      "TCA" = "L", "TCT" = "M", "TCG" = "N", "TCC" = "O",

      "GAA" = "P", "GAT" = "Q", "GAG" = "R", "GAC" = "S",
      "GTA" = "T", "GTT" = "U", "GTG" = "V", "GTC" = "W",
      "GGA" = "X", "GGT" = "Y", "GGG" = "Z", "GGC" = "[",
      "GCA" = "\\", "GCT" = "]", "GCG" = "^", "GCC" = "_",

      "CAA" = "`", "CAT" = "a", "CAG" = "b", "CAC" = "c",
      "CTA" = "d", "CTT" = "e", "CTG" = "f", "CTC" = "g",
      "CGA" = "h", "CGT" = "i", "CGG" = "j", "CGC" = "k",
      "CCA" = "l", "CCT" = "m", "CCG" = "n", "CCC" = "o")

  pieces <- ""
  for (c in 1:length(triplets)) {
    pieces <- paste0(pieces, lookup_table[[triplets[c]]])
  }
  return(pieces)
}
