## REFERENCE PROFILES
# Get cell reference profile data from NanoString
# Use this reference profile as it is the only one available from CosMx data, originally from:
# https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Cell-Profiles/main/Human/IO/IO.profiles.csv
# I have added a name to the first column
io_profiles <- read_csv(
  file = here("Analysis", "metadata", "NanoString.CosMx.Human.IO.profiles.csv"),
  col_types = cols(
    `Gene` = col_character(),
    `B cell` = col_double(),
    `Dendritic cell` = col_double(),
    Endothelial = col_double(),
    Fibroblast = col_double(),
    Macrophage = col_double(),
    `Mast cell` = col_double(),
    Monocyte = col_double(),
    Neutrophil = col_double(),
    `NK cell` = col_double(),
    Plasma = col_double(),
    Plasmablast = col_double(),
    `Plasmacytoid dendritic cell` = col_double(),
    `T cell CD4` = col_double(),
    `T cell CD8` = col_double(),
    `T cell regulatory` = col_double()
  )
)
row_name = io_profiles$Gene
io_profiles %<>% dplyr::select(-Gene) %>% as.matrix
rownames(io_profiles) = row_name
