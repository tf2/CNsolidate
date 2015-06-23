CNsolidate
==========

CNV detection using 12 independent change point detection algorithms and an expert voting system

See example at the bottom for complete walkthrough using synthetic data

# INSTALL
Download zip file or pull master branch
```bash
R CMD build CNsolidate-master
R CMD install CNsolidate_1.2.tar.gz
```

Load library with R in the usual way:

```R
library(CNsolidate)
```

# CONFIGURATION
All execution methods and parameter defintions are controlled by a single configuration object

Get a config object with defualt parameters
```R
set = settings()
```
# INPUTS, OUTPUTS and FORMATS

```R
set$files$file = <YOUR_INPUT_FILE>
set$files$odir = <YOUR_OUTPUT_DIRECTORY>
set$files$format = "fe"
```

Note: format "fe" is for the Agilent feature extraction file format.
Simple "bed" format also supported see walkthrough

# EXECUTION

Run all desired normalisation steps:
```R
for (x in 1:length(set$norm)) {
      if (set$norm[x] == 1) {
        result <- try( do.call( names(set$norm[x]), list(set) ) )
      }
}
```
Run all desired CNV detection algorithms:
```R
for (x in 1:length(set$algorithms)) {
      if (set$algorithms[x] == 1) {
        result <- try( do.call( names(set$algorithms[x]), list(set) ) )
      }
}
```
Run all desired combination, algorithm weighting, post processing and annotation steps:
```R
for (x in 1:length(set$combine)) {
      if (set$algorithms[x] == 1) {
        result <- try( do.call( names(set$combine[x]), list(set) ) )
      }
}
```
# OUTPUT
The final CNV calls are contained in a file with "_FinalReport.txt" extension which will be created in the ouput directory.

# Example
Make temp directory, cd and start R:
```bash
mytmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'` && cd $mytmpdir
R
```
Run a test is R using synthetic data example
```R
# load library
library(CNsolidate)

# generate some test data
outputdir = getwd()
inputname = paste(outputdir, "/testData.txt", sep="")
testchr = 1
testlength = 20000
testmean = 0.001
testsd = 0.2
testproportion = 0.2
testData = syn.genome(data.frame(testchr,testlength,testmean,testsd,testproportion))

# write test data to file note: last column is a weight used for each probe (set all to 1 for this example)
write.table(data.frame(testData$data,1), file=inputname, sep="\t", row.names=F, col.names=F, quote=F)

# set default parameters
set = settings()
set$files$file = inputname
set$files$odir = outputdir
set$files$format = "bed"

# run CNV callers
for (x in 1:length(set$algorithms)) {
      if (set$algorithms[x] == 1) {
        result <- try( do.call( names(set$algorithms[x]), list(set) ) )
      }
}

# combine CNV callers
for (x in 1:length(set$combine)) {
      if (set$combine[x] == 1) {
        result <- try( do.call( names(set$combine[x]), list(set) ) )
      }
}

# read final CNV call list
finalreport = paste(substr(set$files$file, 1, nchar(set$files$file)-4), "_FinalReport.txt", sep="")
cnv_calls = read.table(finalreport)
colnames(cnv_calls) = c("chr", "start", "stop", "meanl2r", "probes", "start_index", "stop_index", "algorithms", "wscore", "adj_wscore", "p_value")

# make a crude plot
plot(testData$data[,4], pch=20, ylim=c(-3,3))
abline(v=cnv_calls$start_index, col="green")
abline(v=cnv_calls$stop_index, col="red")
```



