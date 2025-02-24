## COLOR PALETTES

# Color vector for patients
palette_patients <- c(
  
  "1" = "#ef476f",
  "2" = "#ffd166",
  "3" = "#06d6a0",
  "4" = "#118ab2",
  "5" = "#073b4c",
  "6" = "#84a59d",
  "7" = "#8187dc")


# Color vector for tumor types
palette_tumor_types <- c(
  "ER+" = "#63AC74",
  "HER2+" = "#D5809B",
  "TN" = "#5670B5",
  "ER-" = "#9A4D32"
)

# Color vector for cell types, similar cell types have similar colors
palette_cell_types <- c(
  
  # A type of immune cells that are made in the bone marrow and are found in the blood and in lymph tissue
  "Lymphocytes" = "#30D5C8",
    "B cell" = "#D4819C",
    "NK cell" = "#7790BD",
    "Plasma" = "#955C97",
    "Plasmablast" = "#B3709B",
    "T cell CD4" = "#315A2F",
    "T cell CD8" = "#64AB73",
    "T cell regulatory" = "#526F55",
  
  # Family of immune cells comprising monocytes, macrophages, myeloid dendritic cells (mDCs),
  # granulocytes, and mast cells that originate from a common myeloid progenitor in the bone marrow
  "Myeloids" = "#BF489F",
    "Dendritic cell" = "#5671B3",
    "Macrophage" = "#55274E",
    "Mast cell" = "#C8689F",
    "Monocyte" = "#6E3987",
    "Neutrophil" = "#9C6DA8",
    "Plasmacytoid dendritic cell" = "#212222",
  
  # A type of cell that makes up certain types of connective tissue (supporting tissue that surrounds other tissues and organs)
  "Stromal cells" = "#D48542",
    "Endothelial" = "#EDBD39",
    "Fibroblast" = "#C43522",
    
  "Brain cells" = "#111111",
    # Inhibitory neurons
    "Inhibitory neuron" = "#E3D26F",
    "Inhibitory neuron A" = "#E3D27F",
    "Inhibitory neuron B" = "#E3D28F",
    "Inhibitory neuron C" = "#E3D29F",
    # Excitatory neurons
    "L2/3 neuron" = "#89808C",
    "L4 neuron" = "#89808C",
    "L6 neuron" = "#89808C",
    # Glial cells
    "Oligodendrocyte-like" = "#CEBACF",
    "Oligodendrocyte precursor cell" = "#CEBACF",
    "Oligodendrocyte" = "#CEBACF",
    # Microglia cells
    "Microglia" = "#C08552",
    "Microglia A" = "#C08562",
    "Microglia B" = "#C08572",
    # Astrocyte cells
    "Astrocyte" = "#F9E0D9",
    "Astrocyte A" = "#F9E0E9",
    "Astrocyte B" = "#F9E0F9",
  
  # Other cells
    "Tumor" = "#C0C0C0",
    "Unknown" = "#f5ee59",
    "Other" = "#B0B0B0"
)
