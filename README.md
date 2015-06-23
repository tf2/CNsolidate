CNsolidate
==========

CNV detection using 12 independent change point detection algorithms and an expert voting system

library(CNsolidate)
set = settings()

# EXECUTION

Run all desired normalisation steps:
for (x in 1:length(set$norm)) {
      if (set$algorithms[x] == 1) {
        result <- try( do.call( names(set$algorithms[x]), list(set) ) )
      }
}


Run all desired CNV detection algorithms:
for (x in 1:length(set$algorithms)) {
      if (set$algorithms[x] == 1) {
        result <- try( do.call( names(set$algorithms[x]), list(set) ) )
      }
}

Run all desired combination, algorithm weighting, post processing and annotation steps:
for (x in 1:length(set$algorithms)) {
      if (set$algorithms[x] == 1) {
        result <- try( do.call( names(set$algorithms[x]), list(set) ) )
      }
}


