###### COLORS ######

# Know clusters that should have consistent colors
known_clusters_colors <- c(
  "B cell" = "#5A5156",
  "Dendritic cell" = "#E4E1E3",
  "Endothelial" = "#F6222E",
  "Fibroblast" = "#FE00FA",
  "Macrophage" = "#16FF32",
  "Mast cell" = "#3283FE",
  "Monocyte" = "#FEAF16",
  "Neutrophil" = "#B00068",
  "NK cell" = "#1CFFCE",
  "Plasma" = "#90AD1C",
  "Plasmablast" = "#2ED9FF",
  "Plasmacytoid dendritic cell" = "#DEA0FD",
  "T cell CD4" = "#AA0DFE",
  "T cell CD8" = "#F8A19F",
  "T cell regulatory" = "#325A9B",
  "Tumor" = "#C4451C",
  "T cell" = "#1C8356",
  "Unknown" = "#f5ee59"
)

# Create a vector with a color for every cluster
# It can be used to create a palette which includes known clusters (e.g. B cell, T cell, etc.) and unknown clusters (e.g. 0, 1, ..., a, b, etc.)
generate_colors_lookup_table <- function(
    data,
    cluster_column_name,
    known_clusters_colors = NULL,
    color_palette = DiscretePalette(36, palette = "polychrome")) {
  
  # Remove the colors that are used by the known clusters
  usable_color_palette <- dplyr::setdiff(color_palette, known_clusters_colors)
  
  # Random colors for the unknown clusters
  # The cluster 0 in Seurat for example is different for every patient
  usable_color_palette <- sample(usable_color_palette)
  
  # "data" must be either a data frame or a tibble
  # If the object is of type Seurat, extract the meta data
  if ("Seurat" %in% class(data)) {
    data <- data@meta.data
  }
  
  # Get the unique clusters in the data (works with both base R data frames and tibbles)
  clusters <- unique(dplyr::pull(data, cluster_column_name))
  
  # Remove the known clusters from the list of clusters as they already have a color assigned
  unknown_clusters <- dplyr::setdiff(clusters, names(known_clusters_colors))
  
  # Assign a random colors to the unknown clusters
  unknown_clusters_colors <- setNames(usable_color_palette[1:length(unknown_clusters)], unknown_clusters)
  
  # Return the colors for the known and unknown clusters
  return(c(known_clusters_colors, unknown_clusters_colors))
}

###### GENERATE SINGLE PLOTS ######

# Use this plot to print the level of expression of a protein or rna
generate_feature_plot <- function(patient_data, reduction, features, max_cutoff = NA) {
  
  print(paste("Generate FeaturePlots from", reduction, "reduction of features", features))
  
  patient_num <- get_patient_num(patient_data)
  
  features_plots <- FeaturePlot(
    object = patient_data,
    features = features,
    reduction = reduction,
    max.cutoff = max_cutoff) +
    plot_annotation(
      title = 'Patient 1',
      subtitle = reduction,
    )  & NoLegend() & NoAxes()
  
  features_plot_name <- paste0("Patient_",  patient_num, "_featureplots_", reduction)
  return(setNames(list(features_plots), features_plot_name))
}

# Show the significance of every principal component of the PCA
# It can be used to decide the number of dims of the FindNeighbors function
# By default the RunPCA function uses 50 dimensions, plot all of them
generate_variance_elbow_plots <- function(patient_data, reduction, dims = 50) {
  
  print(paste("Generate ElbowPlot from", reduction))
  
  # Get patient number
  patient_num <- get_patient_num(patient_data)
  
  # Variance explained calculation
  variance_explained <- patient_data@reductions[[reduction]]@stdev^2
  cumulative_variance <- cumsum(variance_explained) / sum(variance_explained)
  
  # Dynamically find the PCs that explain 70%, 80%, and 90% variance (within the specified dims)
  pc_70 <- which(cumulative_variance >= 0.7 & seq_along(cumulative_variance) <= dims)[1]
  pc_80 <- which(cumulative_variance >= 0.8 & seq_along(cumulative_variance) <= dims)[1]
  pc_90 <- which(cumulative_variance >= 0.9 & seq_along(cumulative_variance) <= dims)[1]
  
  # Define dynamic breaks based on number of PCs
  if (dims <= 20) {
    x_breaks <- seq(1, dims, 1)  # Show every PC if <= 20 PCs
  } else if (dims <= 50) {
    x_breaks <- seq(1, dims, 2)  # Show every 2nd PC if <= 50 PCs
  } else {
    x_breaks <- seq(1, dims, 5)  # Show every 5th PC if > 50 PCs
  }
  
  # Common theme for both plots
  common_theme <- theme_minimal(base_size = 14)
  
  # Plot 1: Cumulative variance explained plot (limited to the specified number of PCs)
  cumulative_plot <- ggplot(data.frame(PC = 1:dims, CumulativeVariance = cumulative_variance[1:dims]), aes(x = PC, y = CumulativeVariance)) +
    geom_line(color = "black") +
    geom_point(color = "black") +
    labs(subtitle = "Cumulative Variance Explained", x = "Principal Components", y = "Variance Explained") +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "green") +
    scale_x_continuous(breaks = x_breaks) +  # Dynamic x-axis breaks
    common_theme
  
  # Add vertical lines and annotations only if the PCs exist within the specified range
  if (!is.na(pc_70)) {
    cumulative_plot <- cumulative_plot +
      geom_vline(xintercept = pc_70, linetype = "dotted", color = "purple") +
      annotate("text", x = pc_70, y = 0.7, label = paste("PC", pc_70), color = "purple", hjust = -0.2)
  }
  if (!is.na(pc_80)) {
    cumulative_plot <- cumulative_plot +
      geom_vline(xintercept = pc_80, linetype = "dotted", color = "purple") +
      annotate("text", x = pc_80, y = 0.8, label = paste("PC", pc_80), color = "purple", hjust = -0.2)
  }
  if (!is.na(pc_90)) {
    cumulative_plot <- cumulative_plot +
      geom_vline(xintercept = pc_90, linetype = "dotted", color = "purple") +
      annotate("text", x = pc_90, y = 0.9, label = paste("PC", pc_90), color = "purple", hjust = -0.2)
  }
  
  # Plot 2: Standard elbow plot (limited to the specified number of PCs)
  elbow_plot <- ElbowPlot(patient_data, reduction = reduction, ndims = dims) +
    labs(subtitle = "Elbow Plot", x = "Principal Components", y = "Standard Deviation") +
    scale_x_continuous(breaks = x_breaks) +  # Dynamic x-axis breaks
    common_theme
  
  # Combine the plots vertically using patchwork and add shared title and subtitle
  combined_plot <- (cumulative_plot / elbow_plot) +
    plot_annotation(title = paste("Patient", patient_num),
                    subtitle = paste("Reduction:", reduction))
  
  # Set name to be used when generating an image from the plot
  combined_plot_name <- paste("Patient", patient_num, "combined_elbow_variance_plots", reduction, sep = "_")
  
  return(setNames(list(combined_plot), combined_plot_name))
}

generate_umap <- function(
    patient_data,
    cluster_var,
    cluster_reduction,
    cluster_name = NULL,
    color_lookup_table = NULL) {
  
  if (is.null(cluster_name)) {
    cluster_name <- cluster_var
  }
  
  patient_num <- get_patient_num(patient_data)
  
  # Graphs the output of a dimensional reduction technique on a 2D scatter plot
  # Each point is a cell and it's positioned based on the cell embeddings determined by the reduction technique
  umap_clusters <- Seurat::DimPlot(
    object = patient_data,
    reduction = cluster_reduction,,
    group.by = cluster_var,
    label=TRUE,
    label.box=TRUE,
    repel=TRUE,
    cols = color_lookup_table) +
    labs(
      title = paste("Patient", patient_num),
      subtitle = cluster_name) +
    NoLegend()
  
  # Return the plot inside a list with a name
  umap_clusters_plot_name <- paste("Patient",  patient_num, cluster_var, "umap", sep = "_")
  return(setNames(list(umap_clusters), umap_clusters_plot_name))
}

# The flightpath plot places cells according to their posterior probabilities of belonging to each cell type.
# It conveys the tendency of different cell types to be confused with each other.
generate_flightpath <- function(
    ist_object,
    patient_num,
    color_lookup_table = NULL) {
  
  set.seed(6)
  fp <- flightpath_plot(flightpath_result = NULL,
                        insitutype_result = ist_object,
                        col = color_lookup_table[ist_object$clust],
                        showclusterconfidence = TRUE) +
    labs(title = paste("Patient", patient_num),
         subtitle = "Insitutype Semisupervised Clustering")
  
  fp_name <- paste("Patient",  patient_num, "Insitutype_semisup_clusters_flightpath_plot", sep = "_")
  return(setNames(list(fp), fp_name))
}

####### GENERATE HEATMAP #######

# Heatmap plot, define a function to dynamically calculate label size based on the number of features
calculate_label_size <- function(num_features) {
  
  base_size <- 400
  scaled_size <- base_size / (num_features)
  return(scaled_size)
}

# Define the size parameter (size of the cluster names on top) of the heatmap
# proportionally to the size of the smallest cluster
calculate_clusters_names_size <- function(smallest_cluster_proportion) {
  
  base_size <- 9  # Base size for the heatmap
  # scaling_factor <- 0.1
  return(base_size * exp(log(smallest_cluster_proportion)/3))
}

get_clusters_proportions <- function(patient_data) {
  
  # Calculate the size of each cluster
  cluster_sizes <- table(Idents(patient_data))
  total_cells <- sum(cluster_sizes)
  cluster_proportions <- cluster_sizes / total_cells
  
  return(cluster_proportions)
}

generate_dyn_text_heatmap <- function(
    patient_data,
    cluster_var,
    assay_name,
    cluster_name = NULL,
    color_lookup_table = NULL) {
  
  # If a human friendly name is not given, use the name of the column in the Seurat object
  if (is.null(cluster_name)) {
    cluster_name <- cluster_var
  }
  
  # Use a better palette by default
  if (is.null(color_lookup_table)) {
    color_lookup_table <- generate_colors_lookup_table(patient_data, cluster_var)
  }
  
  # Select the cluster as the identity
  Idents(patient_data) <- cluster_var
  
  most_significant_markers <- find_most_significant_markers(
    patient_data,
    cluster_var,
    assay_name)
  
  # Size of smallest cluster
  smallest_cluster_proportion <- min(get_clusters_proportions(patient_data))
  # Size of cluster label
  cluster_label_size <- calculate_clusters_names_size(smallest_cluster_proportion)
  
  # Calculate the number of features (rows)
  num_features <- length(most_significant_markers$gene)
  # Calculate the label size for the current figure
  features_label_size <- calculate_label_size(num_features)
  
  # Read patient number
  patient_num <- get_patient_num(patient_data)
  
  # Create the heatmap with the scaled text
  diff_expr_genes_heatmap <- DoHeatmap(
    patient_data,
    features = most_significant_markers$gene,
    assay = assay_name,
    label = TRUE,
    size = 2,
    group.colors = color_lookup_table,
  ) + theme(
    axis.text.y = element_text(size = features_label_size),
  ) + labs(
    title = paste("Patient", patient_num, cluster_name),
    subtitle = "Top 10 Differentially Expressed Genes per Cluster"
  )
  
  result_list <- list(
    diff_expr_genes_heatmap,
    most_significant_markers
  )
  
  # Dynamically assign the desired names
  names(result_list) <- c(
    paste("Patient", patient_num, cluster_var, "diff_expr_genes_heatmap", sep = "_"),
    paste("Patient", patient_num, cluster_var, "diff_expr_genes_heatmap_list", sep = "_")
  )
  
  # Return the named list
  return(result_list)
}

###### GENERATE STAMP PLOTS ######

# Print all the stamps associated with a patient using the cluster information to color the cells
generate_spatial_plots <- function(
    patient_data,
    cluster_var,
    cluster_name = NULL,
    color_lookup_table = NULL,
    legend_fill = "ident"
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
        ) +
        guides(
          fill = guide_legend(legend_fill)
        )
      stamp_plot_name <- paste("Patient",  patient_num, cluster_var, "core", curr_core, "stamp", as.character(curr_stamp), sep = "_")
      clustering_plots[[stamp_plot_name]] <- stamp_plot
    }
  }
  
  return(clustering_plots)
}