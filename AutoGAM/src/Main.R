#!/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yangkaixin1/software/miniconda/bin/Rscript
if(T){ # prepare environment
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(mgcv)
  library(purrr)
  library(caret)
  library(MuMIn)
  library(stringi)
  library(parallel)
  library(DataCombine)
  library(lubridate)
  library(tidyverse)
  P <<- rprojroot::find_rstudio_root_file
  R.utils::sourceDirectory(P('src/R'))
}

if(T){ # trim data and fitting GAM
  library(getopt)
  spec = matrix(c(
    'var', 'v', 2, 'character',
    'posneg', 'p', 2, 'character',
    'requantification', 'r', 2, 'character'
  ), byrow=TRUE, ncol=4)
  opt = getopt(spec)
  
  if(is.null(opt$var)){
    inputTrimer(raw_sample_path=posneg, raw_virus_path=requantification)
    getGAMinput()
  } else {
    fitGAM(var)
  }
}

