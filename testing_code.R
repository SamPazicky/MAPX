# library(data.table)
# library(tidyverse)
library(roxygen2)
roxygenise()
# rm(list=ls())
library(devtools)

usethis::use_vignette("my-vignette")

install(build_vignettes = TRUE)
devtools::build()
# for installation: install.packages(path_to_file, repos = NULL, type="source")