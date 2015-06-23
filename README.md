CNsolidate
==========

CNV detection using 12 independent change point detection algorithms and an expert voting system

All execution methods and parameter defintions are controlled by a single configuration object

library(CNsolidate)

# CONFIGURATION

Get a config object with defualt parameters

set = settings()

# EXECUTION

Run all desired normalisation steps:

`for (x in 1:length(set$norm)) {
      if (set$norm[x] == 1) {
        result <- try( do.call( names(set$norm[x]), list(set) ) )
      }
}`


Run all desired CNV detection algorithms:

for (x in 1:length(set$algorithms)) {
      if (set$algorithms[x] == 1) {
        result <- try( do.call( names(set$algorithms[x]), list(set) ) )
      }
}

Run all desired combination, algorithm weighting, post processing and annotation steps:

for (x in 1:length(set$combine)) {
      if (set$algorithms[x] == 1) {
        result <- try( do.call( names(set$combine[x]), list(set) ) )
      }
}


