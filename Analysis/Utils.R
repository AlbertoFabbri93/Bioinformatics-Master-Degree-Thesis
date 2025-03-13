### LIBRARIES
# Simplify saving the plots generated in base R
library("R.devices")

####### SAVE LIST OF ITEMS #######

# This function saves a list of plots and data frames to a specified folder
save_data <- function(data_list, folder_path, img_ext = gl_img_ext) {
  
  # Ensure the folder exists, if not, create it
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  # Flatten the plot/dataframe object list
  flat_items <- flatten_list(data_list)
  
  # Iterate over the object list and save each one
  for (item_name in names(flat_items)) {
    item <- flat_items[[item_name]]
    
    # Save data frames as CSV files
    if (is.data.frame(item)) {
      
      file_path <- file.path(folder_path, paste0(item_name, ".csv"))
      write.csv(item, file_path, row.names = FALSE)
      
    }
    # Handle plot objects (ggplot, trellis, or recordedplot)
    else if (inherits(item, "ggplot") || inherits(item, "trellis") || inherits(item, "recordedplot")) {
      
      # Check if custom dimensions are set as attributes on the plot object
      width <- attr(item, "width")
      height <- attr(item, "height")
      units <- attr(item, "units")
      dpi <- attr(item, "dpi")
      
      # Saving ggplot objects
      if (inherits(item, "ggplot")) {
        
        # Use ggsave() default values if no custom dimensions are provided
        width <- if (!is.null(width)) width else NA
        height <- if (!is.null(height)) height else NA
        units <- if (!is.null(units)) units else "px"
        dpi <- if (!is.null(dpi)) dpi else 300
        
        # Save the image in all the desired formats using a loop
        for (ext in img_ext) {
          file_path <- file.path(folder_path, paste0(item_name, ".", ext))
          ggsave(
            filename = file_path,
            plot = item,
            width = width,
            height = height,
            units = units,
            dpi = dpi)
        }
      }
      # Saving trellis and recordedplot objects
      else {
        
        # For trellis and recordedplot objects, use defaults if custom dimensions are not provided
        units <- if (!is.null(units)) units else "px"
        dpi <- if (!is.null(dpi)) dpi else 72
        width <- if (!is.null(width)) width else 480
        height <- if (!is.null(height)) height else 480
        
        # Save the image in all the desired formats using the R.devices package
        devEval(
          type = img_ext,
          path = folder_path,
          name = item_name,
          width = width,
          height = height,
          units = units,
          res = dpi,
          expr = { replayPlot(item) }
        )
      }
    }
  }
}

####### FLATTEN NESTED PLOT/DATAFRAME LIST #######

# Function to flatten a list of nested plots and data frames
# Images are also lists so simply unflattening the list does not work
flatten_list <- function(obj_list) {
  
  # Initialize an empty list to store the plot/dataframe objects
  flat_list <- list()
  
  # Iterate over the elements of the input list
  for (obj_name in names(obj_list)) {
    # "trellis" refers to a class of objects created by the lattice package
    # If the element is a ggplot, trellis, or recordedplot (base R plot) object, add it to the returned list
    if (inherits(obj_list[[obj_name]], "ggplot") || inherits(obj_list[[obj_name]], "trellis") || inherits(obj_list[[obj_name]], "recordedplot")) {
      flat_list[[obj_name]] <- obj_list[[obj_name]]
    }
    # If the element is a data frame, add it to the returned list
    else if (is.data.frame(obj_list[[obj_name]])) {
      flat_list[[obj_name]] <- obj_list[[obj_name]]
    }
    # If the element is a nested list, recursively flatten it
    else if (is.list(obj_list[[obj_name]])) {
      flat_list <- c(flat_list, flatten_list(obj_list[[obj_name]]))
    }
  }
  # Return a list of plot/dataframe objects
  return(flat_list)
}

####### PRINT A VECTOR COLOR PALETTE #######

print_palette <- function(colors, height = 3, width = 0.5) {
  # Check if colors is unnamed
  is_named <- !is.null(names(colors))
  
  # Get the color values
  color_values <- base::unname(colors)
  
  # Create a blank plot with appropriate dimensions
  n <- length(colors)
  plot(0, 0, type = "n", xlim = c(0, width * 2), ylim = c(0, n * height),
       xlab = "", ylab = "", axes = FALSE)
  
  # Plot each color rectangle and optionally its name
  for (i in 1:n) {
    y_bottom <- (n - i) * height
    y_top <- (n - i + 1) * height
    
    # Draw the colored rectangle
    rect(0, y_bottom, width, y_top, col = color_values[i], border = "black")
    
    # Add the color name only if the palette is named
    if (is_named) {
      text(width * 1.2, (y_bottom + y_top) / 2, names(colors)[i], 
           adj = 0, cex = 1)
    }
  }
}

####### GET PATIENT NUM #######

get_patient_num <- function(patient_data) {
  # The double [[ are necessary to make it return an unnamed number
  patient_num <- patient_data$Patient.ID[[1]]
  return(patient_num)
}

####### GET PATIENT DIRECTORIES #######

get_patient_dir_img <- function(patient_num) {
  patient_dir_img <- here(images_dir, paste0("Patient_", patient_num, "_plots/"))
  return(patient_dir_img)
}

get_patient_dir_rds <- function(patient_num) {
  patient_dir_rds <- here(objects_dir, paste0("Patient_", patient_num, "_data/"))
  return(patient_dir_rds)
}

###### MARKERS ######

# Returns a matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
# It uses a Wilcoxon Rank Sum Test and only returns positive markers
# By default, it returns the top 10 markers per cluster with a log2 fold change greater than 1
find_most_significant_markers <- function(
    patient_data,
    cluster_var,
    assay_name,
    log2fc_threshold = 1,
    genes_per_cluster = 10) {
  
  # Select the cluster as the identity
  Idents(patient_data) <- cluster_var
  
  # Find markers (differentially expressed genes) for each of the identity classes in the filtered dataset
  # If you have more than one assay it is necessary to specify the assay parameter
  markers_data <- FindAllMarkers(
    test.use = "wilcox",
    patient_data,
    assay = assay_name,
    only.pos = TRUE,
    random.seed = 5)
  
  # Filter markers to get the most significant ones per cluster
  most_significant_markers <- markers_data %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > log2fc_threshold) %>%
    slice_head(n = genes_per_cluster) %>%
    ungroup()
  
  return(most_significant_markers)
}

####### INSITUTYPE SEMISUPERVISED #######

run_ist_semisup_extract_data <- function(
    patient_data,
    ...) {
 
  # Get the data necessary to run insitutype
  patient_rna_data <- patient_data[["RNA"]]
  patient_cohort <- get_patient_cohort(patient_data)
  patient_rna_counts <- get_seurat_layer_data(patient_data, "RNA", "counts")
  patient_avg_neg_probes <- get_avg_neg_probes(patient_data[["NegativeProbes"]])
  
  # Call the insitutype semisupervised function
  patient_semisup <- run_Insitutype_semisupervised(
    patient_data = patient_rna_data,
    patient_cohort = patient_cohort,
    patient_rna_counts = patient_rna_counts,
    patient_avg_neg_probes = patient_avg_neg_probes,
    ...)
  
  # Return the Seurat object with the added clusters
  return(patient_semisup)
}

run_Insitutype_semisupervised <- function(
    patient_data,
    patient_cohort,
    patient_rna_counts,
    patient_avg_neg_probes,
    io_profiles = NULL,
    clusts_search_space = 1:8,
    phase_1_size = 20,
    phase_2_size = 50,
    phase_3_size = 200,
    iter_num = 1,
    max_iter = 5) {
  
  if (is.null(io_profiles)) {
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
      ))
    row_name = io_profiles$Gene
    io_profiles %<>% dplyr::select(-Gene) %>% as.matrix
    rownames(io_profiles) = row_name
  }
  
  set.seed(6)
  clusts_num <- chooseClusterNumber(
    counts = patient_rna_counts,
    neg = patient_avg_neg_probes,
    fixed_profiles = io_profiles,
    n_clusts = clusts_search_space)
  
  # Semi-supervised learning with insitutype and reference profiles
  # Insitutype needs integers, if given floating point numbers it fails with misleading errors
  patient_semisup <- insitutype(
    x = patient_rna_counts,
    neg = patient_avg_neg_probes,
    cohort = patient_cohort,
    reference_profiles = io_profiles,
    
    # Enter your own per-cell background estimates here if you
    # have them; otherwise insitutype will use the negprobes to
    # estimate background for you.
    bg = NULL,
    # Group the cells the do not correspond to any type in the reference matrix
    n_clusts = clusts_num$best_clust_number,
    # reference_profiles = updatedprofiles$updated_profiles,
    # Update the reference profile based on the current data
    update_reference_profiles = FALSE,
    # choosing inadvisably low numbers to speed the vignette; using the defaults
    # in recommended.
    # This is the number of cells used in each phase, because of random sampling
    n_phase1 = phase_1_size,
    n_phase2 = phase_2_size,
    n_phase3 = phase_3_size,
    n_starts = iter_num,
    max_iters = max_iter
  )
  
  return(patient_semisup)
}

####### EXTRACT PATIENT DATA #######

extract_patient_data <- function(all_patients_data, patient_num) {
  
  # Extract patient data from global Seurat object
  patient_data <- subset(x = all_patients_data, subset = Patient.ID == patient_num)
  
  return(patient_data)
}
  
get_avg_neg_probes <- function(patient_data) {
  
  # Calculate the average negative probes per cell
  patient_avg_neg_probes <- LayerData(patient_data, layer = "counts") %>%
  as.matrix() %>%
  t() %>%
  Matrix::rowMeans()

  return(patient_avg_neg_probes)
}

get_seurat_layer_data <- function(patient_data, assay_name, layer_name) {
  
  # Extract the data from the specified assay and layer
  patient_assay_layer_data <- GetAssayData(
    patient_data,
    assay = assay_name,
    layer = layer_name) %>%
  as.matrix() %>%
  t()
  # Return a matrix of the data
  return(patient_assay_layer_data)
}

get_patient_cohort <- function(patient_data) {
  
  # Features to be used for the cohorting
  features <- c("Mean.PanCK", "Mean.CD45", "Mean.CD68")
  # Cohort of the patient
  patient_immunofluorescence <- patient_data@meta.data %>% dplyr::select(all_of(features))
  # fastCohorting is stochastic, so set the seed for reproducibility
  set.seed(42);
  # "Gaussian_transform = TRUE" maps variables to gaussians in order to place dramatically different variables on the same scale
  patient_cohort <- fastCohorting(patient_immunofluorescence, gaussian_transform = TRUE, n_cohorts = 5)
  # check clusters and cohort numbers
  table(patient_cohort)
  
  return(patient_cohort)
}

####### PATIENT INFO #######

get_patient_info <- function(patient_data) {
  
  patient_num <- get_patient_num(patient_data)
  # Initialize list to store the cores, stamps, and fovs associated with the patient
  cores <- list()
  # Get cores for the current patient
  patient_cores <- sort(unique(patient_data@meta.data$core_serial))
  
  for (core in patient_cores) {
    stamps <- list()
    # Get stamps for the current core
    patient_core_stamps <- sort(unique(patient_data@meta.data$stamp[patient_data@meta.data$core_serial == core]))
    
    for (stamp in patient_core_stamps) {
      fovs <- table(patient_data@meta.data$fov
                    [patient_data@meta.data$core_serial == core & patient_data@meta.data$stamp == stamp])
      stamps[[stamp]] <- as.list(fovs)
    }
    cores[[core]] <- stamps
  }
  patient_info <- list(patient_num = patient_num, cores = cores)
  return(invisible(patient_info))
}

print_patient_info <- function(patient_data) {
    
  patient_info <- get_patient_info(patient_data)
  patient_num <- patient_info$patient_num
  # Print cores and fovs associated with the patient
  print(paste("FOVs and cell count associated with patient", patient_num))
  
  cores <- patient_info$cores
  for (core in names(cores)) {
    print(paste("----------", "CORE", core, "----------"))
    
    stamps <- cores[[core]]
    # Iterate over stamps using indices
    for (i in seq_along(stamps)) {
      print(paste("-----", "STAMP", i, "-----"))
      # Print FOVs from the stamps list
      fovs <- stamps[[i]]
      
      if (length(fovs) == 0) {
        cat("Missing data\n\n")
      } else {
        # Print FOVs in a table
        cat("FOV\tCell_count\n")
        for (fov_name in names(fovs)) {
          cat(fov_name, "\t", fovs[[fov_name]], "\n")
        }
        cat("\n")
      }
    }
  }
}

####### CLUSTER INFO #######

# Define function to extract cluster information
get_cluster_info <- function(seurat_obj, clustering_column) {
  
  # Check if clustering_column is valid
  if (!clustering_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", clustering_column, "not found in meta.data of Seurat object."))
  }
  
  seurat_obj@meta.data[[clustering_column]] <- as.character(seurat_obj@meta.data[[clustering_column]])
  
  # Get cluster information
  clusters <- unique(seurat_obj@meta.data[[clustering_column]])
  
  # Initialize list to store cluster names and cell counts
  cluster_info <- c()
  
  # Loop through clusters and get cell counts
  for (i in seq_along(clusters)) {
    cluster_name <- clusters[i]
    cell_count <- sum(seurat_obj@meta.data[[clustering_column]] == cluster_name)
    cluster_info[cluster_name] <- cell_count
  }
  
  # Sort cluster_info vector by cell_count in descending order
  cluster_info <- cluster_info[order(cluster_info, decreasing = TRUE)]
  
  return(cluster_info)
}

# Define function to create cluster summary dataframe with dynamic columns
create_cluster_summary_per_stamp <- function(patient_data, clustering_column) {
  
  # Check if clustering_column is valid
  if (!clustering_column %in% colnames(patient_data@meta.data)) {
    stop(paste("Column", clustering_column, "not found in meta.data of Seurat object."))
  }
  
  # Initialize a dataframe to store cluster summary
  cluster_summary <- data.frame(cluster_name = character(),
                                stringsAsFactors = FALSE)
  
  # Get unique clusters
  clusters <- as.vector(unique(patient_data[[clustering_column]][[clustering_column]]))
  
  # Combine two existing columns to create a new column
  dataa <- paste(patient_data$core_serial, patient_data$stamp, sep = "_")
  
  patient_data <- AddMetaData(patient_data, dataa, "combined_id")
  
  # Get unique combinations of core_serial and stamp in the current cluster
  unique_combinations <- sort(unique(dataa))
  
  
  # Loop through clusters
  for (cluster_name in clusters) {
    
    # Subset data for the current cluster
    cluster_data <- patient_data[ , patient_data[[clustering_column]] == cluster_name]
    
    # Get total cell count for the cluster
    cell_count <- nrow(cluster_data@meta.data)
    
    # Create a new row for the cluster in cluster_summary dataframe
    new_row <- data.frame(cluster_name = cluster_name,
                          total_cells = cell_count,
                          stringsAsFactors = FALSE)
    
    # Loop through unique combinations of core_serial and stamp
    for (core_stamp in unique_combinations) {
      
      # Count cells for the current core_stamp combination
      cell_count <- sum(cluster_data$combined_id == core_stamp)
      
      # Add a column for the current core_stamp combination
      new_row[[core_stamp]] <- cell_count
    }
    
    # Append new_row to cluster_summary dataframe
    cluster_summary <- rbind(cluster_summary, new_row)
  }
  
  # Reset row names of cluster_summary dataframe
  rownames(cluster_summary) <- NULL
  
  return(cluster_summary)
}

# Function to create cluster summary tibble with dynamic columns based on Patient.ID
create_cluster_summary_per_patient <- function(patient_data, clustering_column) {
  
  # Check if clustering_column is valid
  if (!clustering_column %in% colnames(patient_data@meta.data)) {
    stop(paste("Column", clustering_column, "not found in meta.data of Seurat object."))
  }
  
  # Extract metadata once to avoid repeated access
  meta_data <- patient_data@meta.data
  
  # Create a table that counts cells for each combination of cluster and Patient.ID
  count_table <- table(meta_data[[clustering_column]], meta_data$Patient.ID)
  
  # Convert the table to a tibble and reshape it
  cluster_summary <- as_tibble(as.data.frame(count_table)) %>%
    # Rename columns for clarity
    dplyr::rename(cluster_name = Var1, Patient.ID = Var2, cell_count = Freq) %>%
    # Spread the Patient.ID into columns
    pivot_wider(names_from = Patient.ID, values_from = cell_count, values_fill = 0) %>%
    # Add a total_cells column summing across patients
    mutate(total_cells = rowSums(dplyr::select(., -cluster_name))) %>%
    # Sort by cluster_name alphabetically
    arrange(cluster_name)
  
  return(cluster_summary)
}

####### ANALYZE RNA #######

normalize_cluster_data <- function(
    patient_data,
    assay = "RNA",
    clust_col_name = "rna_louvain_clusters",
    umap_name = "umap_rna",
    patient_dims = 1:25,
    patient_res = 0.8) {
  
  # Set the assay to the one containing the RNA data
  DefaultAssay(patient_data) <- assay
  # Normalize the count data present in a given assay
  patient_data <- Seurat::NormalizeData(
    object = patient_data,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE)
  # Scales and centers features in the dataset
  patient_data <- Seurat::ScaleData(
    object = patient_data,
    model.use = "linear",
    verbose = FALSE)
  # Detect highly variable genes for the pca
  # Identifies features that are outliers on a 'mean variability plot'
  patient_data <- Seurat::FindVariableFeatures(
    object = patient_data,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE)
  # Run a PCA dimensionality reduction
  patient_data <- Seurat::RunPCA(
    object = patient_data,
    reduction.name = "pca_rna",
    reduction.key = "PCRNA_",
    seed.use = 1,
    verbose = FALSE)
  # Computes the k.param nearest neighbors
  patient_data <- Seurat::FindNeighbors(
    object = patient_data,
    dims = patient_dims,
    reduction = "pca_rna",
    verbose = FALSE)
  # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm
  # Use the resolution parameter to fine tune the number of expected clusters
  patient_data <- Seurat::FindClusters(
    object = patient_data,
    resolution = patient_res,
    algorithm = 1, # 1 = Louvain algorithm
    cluster.name = clust_col_name,
    graph.name = paste0(assay, "_snn"),
    random.seed = 1,
    verbose = FALSE)
  # Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
  patient_data <- Seurat::RunUMAP(
    object = patient_data,
    n.neighbors = 30L,
    dims = patient_dims,
    repulsion.strength = 5,
    reduction = "pca_rna",
    reduction.name = umap_name,
    reduction.key = "UMAPRNA_",
    seed.use = 1,
    verbose = FALSE)
  
  return(patient_data)
}

####### COMPARE CLUSTERING METHODS #######

compare_clustering_methods <- function(patient_rna_data) {
  
  # Read patient number
  patient_num <- patient_rna_data$Patient.ID[1]
  
  # Create the contingency table (table used to study the correlation between the two variables)
  contingency_tab_clusters <- table(
    patient_rna_data$ist_semisup_clusters,
    patient_rna_data$rna_louvain_clusters,
    dnn = c("Insitutype", "Seurat"))
  
  # Convert the table to a data frame for ggplot2
  df_compare_clusters <- as.data.frame(as.table(contingency_tab_clusters))
  df_compare_clusters$Freq <- log10(df_compare_clusters$Freq + 10)
  
  # Plot using ggplot2 with geom_tile
  heatmap_seurat_vs_insitutype <- ggplot(df_compare_clusters, aes(Insitutype, Seurat, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c() +
    labs(fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste("Patient", patient_num),
      subtitle = "Seurat vs Insitutype clusters"
    )
  
  return(heatmap_seurat_vs_insitutype)
}

####### ANALYZE PROTEINS #######

analyze_proteins <- function(patient_data) {
  
  metadata_proteins_columns <- c(
    "Mean.PanCK",
    "Max.PanCK",
    "Mean.CD68",
    "Max.CD68",
    "Mean.Membrane",
    "Max.Membrane",
    "Mean.CD45",
    "Max.CD45",
    "Mean.DAPI",
    "Max.DAPI")
  
  for (col in metadata_proteins_columns) {
    patient_data[[col]] <- as.numeric(as.factor(patient_data@meta.data[[col]]))
  }
  
  # Extract proteins features and replace underscores in column names
  proteins_features <- patient_data@meta.data[, metadata_proteins_columns]
  colnames(proteins_features) <- gsub("_", "-", colnames(proteins_features))
  proteins_matrix <- as.matrix(proteins_features)
  
  # Ensure row names (cell names) match the Seurat object cell names
  rownames(proteins_matrix) <- rownames(patient_data@meta.data)
  
  # Create a new assay with the proteins data
  proteins_assay <- CreateAssay5Object(counts = t(proteins_matrix))
  
  # Add the new assay to the Seurat object with the name "proteins"
  patient_data[["proteins"]] <- proteins_assay
  DefaultAssay(patient_data) <- "proteins"
  
  # Scale the raw metadata values
  patient_data <- ScaleData(
    object = patient_data,
    layer = "counts",
    assay = "proteins",
    features = metadata_proteins_columns,
    do.center = TRUE,
    do.scale = TRUE)
  
  # Run PCA and specify the number of PCs to compute
  # Only 10 features available, use 10 PCs
  n_pcs <- 9
  patient_data <- RunPCA(
    object = patient_data,
    assay = "proteins",
    features = metadata_proteins_columns,
    npcs = n_pcs,
    reduction.name = "pca_proteins",
    reduction.key = "PCPR_",
    seed.use = 2)
  
  # Use the available number of PCs for FindNeighbors and clustering
  patient_data <- FindNeighbors(
    object = patient_data, 
    dims = 1:n_pcs, 
    assay = "proteins",
    reduction = "pca_proteins")
  patient_data <- FindClusters(
    patient_data,
    resolution = 0.5,
    assay = "proteins",
    cluster.name = "protein_clusters",
    graph.name = "proteins_snn",
    random.seed = 2)
  
  # Run UMAP for visualization
  patient_data <- RunUMAP(
    patient_data,
    assay = "proteins",
    dims = 1:n_pcs,
    reduction = "pca_proteins",
    reduction.name = "umap_proteins",
    reduction.key = "UMAPPR_",
    seed.use = 2)
  
  return(patient_data)
  
}

####### PRINT CLUSTERING PLOTS #######
  
generate_clustering_plots <- function(
    patient_data,
    cluster_var,
    cluster_assay,
    cluster_reduction = "umap",
    create_heatmap = FALSE,
    cluster_name = NULL,
    color_lookup_table = NULL
    ) {
  
  # If a human friendly name is not given, use the name of the column in the Seurat object
  if (is.null(cluster_name)) {
    cluster_name <- cluster_var
  }
  
  # If no color lookup table is given create one
  if (is.null(color_lookup_table)) {
    color_lookup_table <- generate_colors_lookup_table(patient_data, cluster_var)
  }

  # Get a pointer to the spatial representation of the data
  patient_image <- Images(patient_data)[1]
  
  # List to be returned with all the plots
  clustering_plots <- list()
  
  patient_num <- get_patient_num(patient_data)

  # Print some information about the clusters
  print(paste(cluster_name, "and number of cells in each of them associated with patient", patient_num))
  print(table(patient_data[[cluster_var]]))
  
  # Select the cluster as the identity
  Idents(patient_data) <- cluster_var
  # Plot the cells using their polygonal boundaries
  DefaultBoundary(patient_data[[patient_image]]) <- "segmentation"
  
  # Plot the umap of the clusters
  clustering_plots <- c(clustering_plots, generate_umap(
    patient_data,
    cluster_var, 
    cluster_reduction,
    cluster_name,
    color_lookup_table))
  
  # Plot cells in their spatial context
  stamps_list <- list()
  # Loop over every core associated with the patient
  for(curr_core in sort(unique(patient_data@meta.data$core_serial))) {
    # Loop over every stamp associated with the current core 
    # A core can be associated with multiple stamps (e.g. different areas of the same metastasis)
    for (curr_stamp in sort(unique(patient_data@meta.data$stamp[patient_data@meta.data$core_serial == curr_core]))) {
      
      # Subset data from current core and stamp
      core_stamp_subset <- subset(patient_data, subset = core_serial == curr_core & stamp == curr_stamp)
      
      # Plot the current core/stamp combination with all the cells in their spatial context and colored by cluster
      stamp_plot <- ImageDimPlot(
        core_stamp_subset,
        fov = patient_image,
        # Set border color to 'NA' as 'white' masks all cells when zoomed out
        border.color = NA,
        flip_xy = FALSE,
        cols = color_lookup_table) + theme(
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.5, 'lines'), # Adjust the size of the legend keys
          legend.spacing = unit(0.5, 'lines') # Adjust the spacing between legend items
        ) +
        labs(
          title = paste("Patient", patient_num, "Core", curr_core, ", Stamp", curr_stamp),
          subtitle = cluster_name
        )
      stamp_plot_name <- paste("Patient",  patient_num, cluster_var, "core", curr_core, "stamp", as.character(curr_stamp), sep = "_")
      clustering_plots[[stamp_plot_name]] <- stamp_plot
    }
  }
  
  # Create heatmap of most differentially expressed genes
  if (create_heatmap) {
    diff_expr_genes_heatmap <- generate_dyn_text_heatmap(
      patient_data,
      cluster_var,
      cluster_assay,
      color_lookup_table = color_lookup_table,
      cluster_name = cluster_name)
    # Save plot to list
    clustering_plots <- c(clustering_plots, diff_expr_genes_heatmap)
  }
  
  return(clustering_plots)
}

### EXTRACT/REMOVE CLUSTERS ###

# Return a new Seurat object with only the cells that belong to the specified clusters of the given metadata column
extract_clusters <- function(patient_data, cluster_col, cluster_vals) {
  
  # Identify cells that belong to the clusters specified in cluster_vals
  cells_to_keep <- patient_data[[cluster_col, drop=TRUE]] %in% cluster_vals
  
  # Subset the Seurat object to keep only the desired cells
  patient_data_clusters <- subset(patient_data, cells = Cells(patient_data)[cells_to_keep])
  
  return(patient_data_clusters)
}

# Remove a set of clusters from a Seurat object 
# given the name of the metadata column containing the clusters and the name of the clusters to be removed
remove_clusters <- function(patient_data, cluster_col, cluster_vals) {
  
  # Identify cells that do not belong to the clusters specified in cluster_vals
  cells_to_keep <- !patient_data[[cluster_col, drop=TRUE]] %in% cluster_vals
  
  # Subset the Seurat object to keep only the desired cells
  patient_data_clusters_removed <- subset(patient_data, cells = Cells(patient_data)[cells_to_keep])
  
  return(patient_data_clusters_removed)
}

####### GENERATE PLOTS #######

generate_proteins_plots <- function(patient_data, assay) {
  
  print("Generate proteins plots")
  
  patient_num <- get_patient_num(patient_data)
  
  # List with all the plots to be returned
  plot_list <- list()
  
  # Show the significance of every principal component of the PCA
  # It can be used to decide the number of dims of the FindNeighbors function
  elbow_plot_red = "pca_proteins"
  proteins_elbow_plot<- generate_elbow_plot(patient_data, elbow_plot_red, 9)
  plot_list <- c(plot_list, proteins_elbow_plot)
  
  proteins_features_plots <- generate_feature_plot(
    patient_data = patient_data,
    reduction = "umap_proteins",
    features = c("Mean.PanCK", "Mean.CD45", "Mean.CD68", "Mean.Membrane", "Mean.DAPI", "Area" ),
    max_cutoff = "q95")
  plot_list <- c(plot_list, proteins_features_plots)
  
  protein_cluster_var <- "protein_clusters"
  protein_color_lookup_table <- generate_colors_lookup_table(patient_data, protein_cluster_var, known_clusters_colors)
  protein_cluster <- generate_clustering_plots(
    patient_data,
    protein_cluster_var,
    cluster_assay = assay, 
    cluster_reduction = "umap_proteins",
    create_heatmap = FALSE,
    cluster_name = "Protein Clusters",
    color_lookup_table = protein_color_lookup_table)
  plot_list[[protein_cluster_var]] <- protein_cluster
  
  return(plot_list)
}

generate_rna_plots <- function(patient_data, assay, RNA_cluster_var) {
  
  print("Generate RNA plots")
  
  # List with all the plots to be returned
  plot_list <- list()
  
  patient_num <- get_patient_num(patient_data)
  
  # Show the significance of every principal component of the PCA
  # It can be used to decide the number of dims of the FindNeighbors function
  elbow_plot_red = "pca_rna"
  # By default the RunPCA function uses 50 dimensions, plot all of them
  RNA_elbow_plot <- generate_elbow_plot(patient_data, elbow_plot_red, 50)
  plot_list <- c(plot_list, RNA_elbow_plot)
  
  RNA_features_plots <- generate_feature_plot(
    patient_data = patient_data,
    reduction = "umap_rna",
    features = c("Mean.PanCK", "Mean.CD45", "Mean.CD68", "Mean.Membrane", "Mean.DAPI", "Area" ),
    max_cutoff = "q100")
  plot_list <- c(plot_list, RNA_features_plots)
  
  RNA_color_lookup_table <- generate_colors_lookup_table(patient_data, RNA_cluster_var, known_clusters_colors)
  RNA_cluster <- generate_clustering_plots(
    patient_data,
    RNA_cluster_var,
    cluster_assay = assay, 
    cluster_reduction = "umap_rna",
    create_heatmap = TRUE,
    cluster_name = NULL,
    color_lookup_table = RNA_color_lookup_table)
  plot_list[[RNA_cluster_var]] <- RNA_cluster
  
  return(plot_list)
}
 
generate_ist_plots <- function(patient_data, assay, ist_object) {
  
  print("Generate Insitutype plots")
  
  # List with all the plots to be returned
  Insitutype_clusters <- list()
  patient_num <- get_patient_num(patient_data)
  
  Insitutype_cluster_var <- "ist_semisup_clusters"
  Insitutype_color_lookup_table <- generate_colors_lookup_table(
    patient_data,
    Insitutype_cluster_var,
    known_clusters_colors)
  
  Insitutype_clusters <- c(
    Insitutype_clusters,
    generate_flightpath(ist_object, Insitutype_color_lookup_table, patient_num))
 
  Insitutype_clusters <- c(Insitutype_clusters, generate_clustering_plots(
    patient_data,
    Insitutype_cluster_var,
    cluster_assay = assay, 
    cluster_reduction = "umap_RNA",
    create_heatmap = TRUE,
    cluster_name = "Insitutype Semisupervised Clusters",
    color_lookup_table = Insitutype_color_lookup_table))
  
  return(Insitutype_clusters)
}

generate_comparison_plots <- function(patient_data) {
  
  print("Generate Comparison plots")
  
  # List with all the plots to be returned
  plot_list <- list()
  
  patient_num <- get_patient_num(patient_data)

  seurat_vs_insitutype_plot_name <- paste0("Patient_",  patient_num, "_heatmap_seurat_vs_insitutype")
  seurat_vs_insitutype_plot <- compare_clustering_methods(patient_data)
  plot_list[[seurat_vs_insitutype_plot_name]] <- seurat_vs_insitutype_plot
  
  return(plot_list)
}

####### STATISTICAL FUNCTIONS #######

# Calculate the trimean of vector of values
trimean <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  result <- (q[1] + 2 * q[2] + q[3]) / 4
  return(base::unname(result))
}
