
############ First.lib ###############

.onLoad <- function(lib, pkg){
   library.dynam("CNsolidate", pkg, lib)
}

.onAttach <- function (lib, pkg){
  data(genomicInfo)  ##RPR automatic load of the genome info...  
}

.onUnload <- function(libpath)
    library.dynam.unload("CNsolidate", libpath)


############ End of .First.lib ###############


