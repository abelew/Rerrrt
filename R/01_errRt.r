## The following was taken from ggplot2's ggplot2.r
## I presume it is a blanket importer cue for roxygen2 to add
## import statements to the NAMESPACE file so that when ggplot2 is used
## it will ensure that these libraries are available.
## I checked the roxygen documentation and it appears that
## imports are saved as the exclusive set, as a result repeating these
## at each function declaration serves to make explicit what each function
## requires while not (I think) adding excessive cruft to the NAMESPACE

#' Pipe operator
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' errRt: Some functions for analyzing sequencing data intended to quantify RT error rates.
#'
#' This package (along with its perl counterpart) attempts to make it easier to
#' analyze high-throughput experiments which quantify error rates by using a
#' sample metadata worksheet to automate the filtering, classification, and
#' counting of errors in the sequencing data of a known template.
#'
#' @docType package
#' @name errRt
#' @import methods
#' @import data.table
#' @import glue
#' @import ggplot2
#' @import dplyr
#' @import utils
NULL
