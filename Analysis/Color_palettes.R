## COLOR PALETTES
# Color vector for patients
patient_colors <- c(
  "1" = "#D5809B",
  "2" = "#C23621",
  "3" = "#5670B5",
  "4" = "#53274E",
  "5" = "#A5BBD5",
  "6" = "#315A2F",
  "7" = "#63AC74")

# Color vector for cell types, similar cell types have similar colors
cell_type_colors <- c(
  "Tumor" = "#C0C0C0",
  # Lymphocites
  # A type of immune cells that are made in the bone marrow and are found in the blood and in lymph tissue
  "B cell" = "#D4819C",
  "NK cell" = "#7790BD",
  "Plasma" = "#955C97",
  "Plasmablast" = "#B3709B",
  "T cell CD4" = "#315A2F",
  "T cell CD8" = "#64AB73",
  "T cell regulatory" = "#526F55",
  # Myeloid cells
  # Family of immune cells comprising monocytes, macrophages, myeloid dendritic cells (mDCs),
  # granulocytes, and mast cells that originate from a common myeloid progenitor in the bone marrow
  "Dendritic cell" = "#5671B3",
  "Macrophage" = "#55274E",
  "Mast cell" = "#C8689F",
  "Monocyte" = "#6E3987",
  "Neutrophil" = "#9C6DA8",
  "Plasmacytoid dendritic cell" = "#212222",
  # Stromal cells
  # A type of cell that makes up certain types of connective tissue (supporting tissue that surrounds other tissues and organs)
  "Endothelial" = "#EDBD39",
  "Fibroblast" = "#C43522",
  # Others
  "Unknown" = "#f5ee59")
