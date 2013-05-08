
.First.lib <- function(lib,pkg) {
   library.dynam("cnvSource",pkg,lib)
   cat("roots 0.1-1 loaded\n")
}

