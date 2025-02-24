## REFERENCE PROFILES
# Get cell reference profile data from NanoString
# Use these reference profiles as they have been created from a CosMx dataset
# I have added the name "Gene" to the first column of these reference profiles

# CURRENT DATASET ---------------------------------------------------------
# # Get the list of genes from the current dataset
# seurat_obj <- breast_cancer_patients_filt_cells
# # Extract the count assay data
# counts.mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
# # Extract the list of genes
# current_dataset_genes <- rownames(counts.mat)
# # Length of the list of genes
# current_dataset_genes_length <- length(breast_cancer_patients_genes)

# IO PROFILE ------------------------------------------------------
# https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Cell-Profiles/main/Human/IO/IO.profiles.csv
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
# Extract the list of genes
io_genes_list = io_profiles$Gene
# Length of the list of genes
io_genes_list_length = length(io_genes_list)
# Convert the profiles to a matrix with named rows
io_profiles_matrix <- io_profiles %>% dplyr::select(-Gene) %>% as.matrix
rownames(io_profiles_matrix) = io_genes_list

# # List of genes that are in the current dataset but not in the IO profiles
# missing_io <- dplyr::setdiff(current_dataset_genes, io_genes_list)
# # Length of the list of genes that are in the current dataset but not in the IO profiles
# missing_io_length <- length(missing_io)

# BRAIN PROFILE -----------------------------------------------------------
# https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Cell-Profiles/refs/heads/main/Human/Brain/Brain.profiles.csv
brain_profiles <- read_csv(
  file = here("Analysis", "metadata", "NanoString.CosMx.Human.Brain.profiles.csv"),
  col_types = cols(
    `Gene` = col_character(),
    `Astrocyte A` = col_double(),
    `Astrocyte B` = col_double(),
    Endothelial = col_double(),
    `Inhibitory neuron A` = col_double(),
    `Inhibitory neuron B` = col_double(),
    `Inhibitory neuron C` = col_double(),
    `L2/3 neuron` = col_double(),
    `L4 neuron` = col_double(),
    `L6 neuron` = col_double(),
    `Microglia A` = col_double(),
    `Microglia B` = col_double(),
    `Oligodendrocyte` = col_double(),
    `Oligodendrocyte precursor cell` = col_double(),
  )
)
# Extract the list of genes
brain_genes_list = brain_profiles$Gene
# Length of the list of genes
brain_genes_list_length = length(brain_genes_list)
# Convert the profiles to a matrix with named rows
brain_profiles_matrix <- brain_profiles %>% dplyr::select(-Gene) %>% as.matrix
rownames(brain_profiles_matrix) = brain_genes_list

# # List of genes that are in the current dataset but not in the Brain profiles
# missing_brain <- dplyr::setdiff(current_dataset_genes, brain_genes_list)
# # Length of the list of genes that are in the current dataset but not in the Brain profiles
# missing_brain_length <- length(missing_brain)

# MERGED IO & BRAIN PROFILES ----------------------------------------------
# Find the shared genes between the two reference profiles
io_brain_shared_genes_list = dplyr::intersect(io_genes_list, brain_genes_list)
# Length of the list of shared genes
io_brain_shared_genes_list_length = length(io_brain_shared_genes_list)

# Align genes in the IO and Brain profiles
io_profiles_matrix_adj <- io_profiles_matrix[io_brain_shared_genes_list,]
brain_profiles_matrix_adj <- brain_profiles_matrix[io_brain_shared_genes_list,]

# Scale each dataset so that the 99th percentile value becomes 1000,  
# ensuring comparability while preserving relative differences.  
# The 99th percentile is the value below which 99% of the data falls,  
# helping to reduce the influence of extreme values.
io_profiles_matrix_adj <- io_profiles_matrix_adj / quantile(io_profiles_matrix_adj, 0.99) * 1000
brain_profiles_matrix_adj <- brain_profiles_matrix_adj / quantile(brain_profiles_matrix_adj, 0.99) * 1000

# Find shared cell types
shared_cell_types <- dplyr::intersect(colnames(io_profiles_matrix_adj), colnames(brain_profiles_matrix_adj))
# Rename columns
colnames(io_profiles_matrix_adj)[colnames(io_profiles_matrix_adj) == "Endothelial"] <- "Endothelial IO"
colnames(brain_profiles_matrix_adj)[colnames(brain_profiles_matrix_adj) == "Endothelial"] <- "Endothelial Brain"
# Select from which profile shared cell types should be removed
omit_from_brain <- c()
omit_from_io <- c()

# Merge the IO and Brain profiles using the common genes
io_brain_profiles_matrix <- cbind(
  io_profiles_matrix_adj[, dplyr::setdiff(colnames(io_profiles_matrix_adj), omit_from_io)], 
  brain_profiles_matrix_adj[, dplyr::setdiff(colnames(brain_profiles_matrix_adj), omit_from_brain)])
# Sort columns alphabetically
io_brain_profiles_matrix <- io_brain_profiles_matrix[,order(colnames(io_brain_profiles_matrix))]

# # List of genes that are in the current dataset but not in the IO and Brain profiles
# missing_io_brain <- dplyr::setdiff(current_dataset_genes, io_brain_common_genes_list)
# # Length of the list of genes that are in the current dataset but not in the IO and Brain profiles
# missing_io_brain_length <- length(missing_io_brain)
