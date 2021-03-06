% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify.r
\name{classify_sequences}
\alias{classify_sequences}
\title{Given a set of sequences containing mutations and identical reads along with
a summary of how many reads/index are identical/mutant, classify them as
identical indexes, RT indexes, or sequencer indexes.}
\usage{
classify_sequences(
  chng,
  ident,
  reads_per_index_summary,
  min_sequencer = 10,
  verbose = verbose
)
}
\arguments{
\item{chng}{Set of reads with mutants.}

\item{ident}{Set of reads which are identical.}

\item{reads_per_index_summary}{Table containing how many reads/index are
identical/mutant.}

\item{min_sequencer}{Define how stringent one must be to deem a given
read/index group as a sequencer-based error.  Too small and gets too many
false positives, too large and one gets insufficient data to play with.}

\item{verbose}{Print some information while running.}
}
\value{
Groups of reads which are RT, identical, and sequencer; along with
  summary information.
}
\description{
The reads_per_index_summary provides a table of the number of mutant and
identical reads associated with every index in the data.  Identical reads are
(obviously) ones with 0 mutations.  Sequencer index/reads are when there is 1
mutant out of n reads for an index.  RT index/reads are when all reads for an
index are mutant.
}
