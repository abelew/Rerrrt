Rerrrt
------

# Overview

This provides a series of functions in R intended to make it easier to
categorize and visualize the data generated as per the paper:
"Accurate fidelity analysis of the reverse transcriptase by a modified
next-generation sequencing."

It therefore provides a series of filters to remove likely spurious
variants.  Once those are removed, it provides metrics to classify
variants which originated from the experiment as opposed to those
which originated from the sequencer.  Given these variants, it
classifies them according to mutation type (strong to weak, puring to
pyrimidine, etc).  Finally, it contains a series of functions to plot
the results and print them to an xlsx file.

A complete run of these tasks is provided as a vignette, accessible
via:

```{r vignette}
library(Rerrrt)
browseVignettes("Rerrrt")
```

For the interested, the full catalog of raw data is available here:

https://rawerrors.umiacs.io/error_rate_raw_and_processed_data.tar

and should be available via SRA shortly.
