# Use the renv package to make the environment reproducible
source("renv/activate.R")
# Configure renv to use a different library path for each branch in the git repository
# branch <- system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
# Sys.setenv(RENV_PATHS_LIBRARY = file.path("renv/library/branches", branch))

# Increase the maximum printed elements to a number > RNA + SystemControl
options(max.print = 1500)

## DECLARE GLOBAL VARIABLES
# These variables are declared here so it is possible to read, and therefore share, them between different R scripts or Quarto files
gl_obj_dir <- "Analysis/Objects"
gl_img_dir <- "Analysis/Images"
gl_img_ext <- c(".png")
