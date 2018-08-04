source("R/packages.R")  # Load all the packages you need.
source("R/functions.R") # Load all the functions into your environment.
source("R/plan.R")      # Build your workflow plan data frame.

pkgconfig::set_config("drake::strings_in_dots" = "literals")
knitr::opts_knit$set(root.dir = "docs")
knitr::opts_knit$set(base.dir = "docs")

#make(scholarplan)
#make(plan)
make(familyrun)
