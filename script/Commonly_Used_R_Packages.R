# Define a function to check and install missing packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

# Define the list of required R packages
packages = c(
  "ggplot2", "vegan", "reshape2", "xlsx", "readxl", 
  "dplyr", "stringr", "xlsx", "multcompView", "tidyr", 
  "ggpubr"
)