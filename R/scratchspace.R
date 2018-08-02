# library(googlesheets)
# beetles <- googlesheets::gs_read(googlesheets::gs_title("Coleoptera_fams_subfams_01Aug2018_NPL"))
# beetles$ScholarCount <- 0
# for(i in sequence(nrow(beetles))) {
#   if(!is.na(beetles$Family[i])) {
#     try(beetles$ScholarCount[i] <- get_counts_from_scholar(beetles$Family[i]))
#     print(paste(beetles$Family[i], beetles$ScholarCount[i]))
#     Sys.sleep(10)
#     write.csv(beetles, file="~/Downloads/Beetles.csv")
#   }
# }

source("packages.R")
source("functions.R")
get_figures_from_plosone(dir="plosone")
