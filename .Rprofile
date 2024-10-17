source("renv/activate.R")
# Increase the maximum printed elements to a number > RNA + SystemControl
options(max.print = 1500)
# Configure renv to use a different library path for each branch in the git repository
branch <- system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
Sys.setenv(RENV_PATHS_LIBRARY = file.path("renv/library/branches", branch))
