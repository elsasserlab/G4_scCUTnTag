# NOTE: This script generates scRNA_Seq-Zeisel_et_al-neuron.Rds from the
# level-2 neuron loom files (l2_neurons_*.agg.loom) of the Zeisel et al. (2018)
# mouse nervous system atlas (Cell 174:999-1014, DOI:
# 10.1016/j.cell.2018.06.021, http://mousebrain.org). Those loom files are no
# longer publicly hosted (mousebrain.org is offline), so the script cannot run
# as-is. It is kept here to document exactly how the Rds was produced. To run
# it, place the loom files in this same directory; the Rds and norm data are
# written here as well.
print("Load R packages")
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("devtools",
               "loomR",
               "Seurat",
               "Signac",
               "tidyverse",
               "RColorBrewer")

# script directory: loom inputs are read from here and outputs are written here
args = commandArgs(trailingOnly = FALSE)
file_arg = sub("^--file=", "", args[grep("--file=", args)])
SCRIPT_DIR = if (length(file_arg) && nzchar(file_arg) && file.exists(file_arg)) {
  dirname(normalizePath(file_arg))
} else {
  getwd()
}

# result folder
result_folder = paste0(SCRIPT_DIR, "/")
dir.create(result_folder, showWarnings = FALSE, recursive = TRUE)

# helper function for creating Seurat objects
make_seurat_object = function(label, lfile) {
  lfile = connect(filename = lfile,
                  mode = "r+",
                  skip.validate = TRUE)
  counts = lfile[["matrix"]][, ]
  colnames(counts) = lfile[["row_attrs/Gene"]][]
  
  
  row_names = character()
  for (i in seq(1, length(lfile[["col_attrs/Class"]][]))) {
    row_names = c(row_names, paste0(label, "_", as.character(i)))
  }
  
  rownames(counts) = row_names
  counts = t(counts)
  counts = counts[!duplicated(rownames(counts)), ]
  
  seurat = CreateSeuratObject(
    counts = counts,
    project = label,
    min.cells = 3,
    min.features = 200
  )
  seurat@meta.data = seurat@meta.data %>% mutate(orig.ident = label)
  
  return(seurat)
}

# loom files (level 2) from http://mousebrain.org/ (Linnarsson Lab, Zeisel et al. 2018)
loom_path = function(f) file.path(SCRIPT_DIR, f)
lfiles = list.files(
  SCRIPT_DIR,
  full.names = TRUE,
  pattern = "*.loom"
)
amygdala = make_seurat_object(label = "neuron_amygdala", lfile = loom_path("l2_neurons_amygdala.agg.loom"))
cerebellum = make_seurat_object(label = "neuron_cerebellum", lfile = loom_path("l2_neurons_cerebellum.agg.loom"))
cortex1 = make_seurat_object(label = "neuron_cortex1", lfile = loom_path("l2_neurons_cortex1.agg.loom"))
cortex2 = make_seurat_object(label = "neuron_cortex2", lfile = loom_path("l2_neurons_cortex2.agg.loom"))
cortex3 = make_seurat_object(label = "neuron_cortex3", lfile = loom_path("l2_neurons_cortex3.agg.loom"))
drg = make_seurat_object(label = "neuron_drg", lfile = loom_path("l2_neurons_drg.agg.loom"))
enteric = make_seurat_object(label = "neuron_enteric", lfile = loom_path("l2_neurons_enteric.agg.loom"))
hippocampus = make_seurat_object(label = "neuron_hippocampus", lfile = loom_path("l2_neurons_hippocampus.agg.loom"))
hypothalamus = make_seurat_object(label = "neuron_hypothalamus", lfile = loom_path("l2_neurons_hypothalamus.agg.loom"))
medulla = make_seurat_object(label = "neuron_medulla", lfile = loom_path("l2_neurons_medulla.agg.loom"))
midbraindorsal = make_seurat_object(label = "neuron_midbraindorsal", lfile = loom_path("l2_neurons_midbraindorsal.agg.loom"))
midbrainventral = make_seurat_object(label = "neuron_midbrainventral", lfile = loom_path("l2_neurons_midbrainventral.agg.loom"))
striatumdorsal = make_seurat_object(label = "neuron_striatumdorsal", lfile = loom_path("l2_neurons_striatumdorsal.agg.loom"))
olfactory = make_seurat_object(label = "neuron_olfactory", lfile = loom_path("l2_neurons_olfactory.agg.loom"))
striatumventral = make_seurat_object(label = "neuron_striatumventral", lfile = loom_path("l2_neurons_striatumventral.agg.loom"))
pons = make_seurat_object(label = "neuron_pons", lfile = loom_path("l2_neurons_pons.agg.loom"))
sympathetic = make_seurat_object(label = "neuron_sympathetic", lfile = loom_path("l2_neurons_sympathetic.agg.loom"))
thalamus = make_seurat_object(label = "neuron_thalamus", lfile = loom_path("l2_neurons_thalamus.agg.loom"))

# Seurat workflow
seurat_neurons = merge(
  x = amygdala,
  y = list(
    cerebellum,
    cortex1,
    cortex2,
    cortex3,
    drg,
    enteric,
    hippocampus,
    hypothalamus,
    medulla,
    midbraindorsal,
    midbrainventral,
    striatumdorsal,
    olfactory,
    striatumventral,
    pons,
    sympathetic,
    thalamus
  )
)
all.genes = rownames(seurat_neurons)
neuron_rna = NormalizeData(seurat_neurons,
                           normalization.method = "LogNormalize",
                           scale.factor = 10000)
neuron_rna = FindVariableFeatures(neuron_rna,
                                  selection.method = "vst",
                                  nfeatures = 2000)
neuron_rna = ScaleData(neuron_rna, features = all.genes)
neuron_rna = RunPCA(neuron_rna, features = VariableFeatures(object = neuron_rna))
neuron_rna = RunUMAP(neuron_rna, dims = 1:20, return.model = TRUE)

# export normalized data
data.neuron_amygdala = neuron_rna@assays$RNA@layers$data.neuron_amygdala
data.neuron_cerebellum = neuron_rna@assays$RNA@layers$data.neuron_cerebellum
data.neuron_cortex1 = neuron_rna@assays$RNA@layers$data.neuron_cortex1
data.neuron_cortex2 = neuron_rna@assays$RNA@layers$data.neuron_cortex2
data.neuron_cortex3 = neuron_rna@assays$RNA@layers$data.neuron_cortex3
data.neuron_drg = neuron_rna@assays$RNA@layers$data.neuron_drg
data.neuron_enteric = neuron_rna@assays$RNA@layers$data.neuron_enteric
data.neuron_hippocampus = neuron_rna@assays$RNA@layers$data.neuron_hippocampus
data.neuron_hypothalamus = neuron_rna@assays$RNA@layers$data.neuron_hypothalamus
data.neuron_medulla = neuron_rna@assays$RNA@layers$data.neuron_medulla
data.neuron_midbraindorsal = neuron_rna@assays$RNA@layers$data.neuron_midbraindorsal
data.neuron_midbrainventral = neuron_rna@assays$RNA@layers$data.neuron_midbrainventral
data.neuron_striatumdorsal = neuron_rna@assays$RNA@layers$data.neuron_striatumdorsal
data.neuron_olfactory = neuron_rna@assays$RNA@layers$data.neuron_olfactory
data.neuron_striatumventral = neuron_rna@assays$RNA@layers$data.neuron_striatumventral
data.neuron_pons = neuron_rna@assays$RNA@layers$data.neuron_pons
data.neuron_sympathetic = neuron_rna@assays$RNA@layers$data.neuron_sympathetic
data.neuron_thalamus = neuron_rna@assays$RNA@layers$data.neuron_thalamus

layers = c(
  "data.neuron_amygdala",
  "data.neuron_cerebellum",
  "data.neuron_cortex1",
  "data.neuron_cortex2",
  "data.neuron_cortex3",
  "data.neuron_drg",
  "data.neuron_enteric",
  "data.neuron_hippocampus",
  "data.neuron_hypothalamus",
  "data.neuron_medulla",
  "data.neuron_midbraindorsal",
  "data.neuron_midbrainventral",
  "data.neuron_olfactory",
  "data.neuron_striatumventral",
  "data.neuron_pons",
  "data.neuron_sympathetic",
  "data.neuron_thalamus",
  "data.neuron_striatumdorsal"
)

norms = lapply(layers, function(x) {
  GetAssayData(object = neuron_rna, layer = x)
})

inner_join_by_rowname <- function(dfs) {
  # Start with the first data frame
  result <- dfs[[1]]
  
  for (i in 2:length(dfs)) {
    # Match row names (gene symbols)
    result <- merge(result, dfs[[i]], by = "row.names", all = FALSE)
    
    # Restore row names and remove the temporary Row.names column
    rownames(result) <- result$Row.names
    result <- result[, -which(names(result) == "Row.names")]
  }
  return(result)
}

norm_data = inner_join_by_rowname(norms)
norm_data = norm_data %>% mutate(gene_symbol = rownames(norm_data)) %>%
  dplyr::select(gene_symbol, everything())

# UMAP
colors = c(
  "#c994c7",
  "#feb24c",
  "#deebf7",
  "#9ecae1",
  "#3182bd",
  "#636363",
  "#c51b8a",
  "#e5f5e0",
  "#fc9272",
  "#31a354",
  "#a1d99b",
  "#ffffb2",
  "#bdbdbd",
  "#8c2d04",
  "#d94801",
  "#016450",
  "#8c6bb1",
  "#e31a1c"
)

umap = DimPlot(
  object = neuron_rna,
  group.by = "orig.ident",
  reduction = "umap",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = FALSE
) +
  xlim(-20, 20) +
  ylim(-20, 20) +
  scale_color_manual(values = colors) +
  ggtitle("Zeisel et al. - neuron scRNA-Seq") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
umap

ggsave(
  glue("{result_folder}Zeisel_et_al-neuron_scRNA_Seq.pdf"),
  plot = umap,
  width = 10,
  height = 10,
  device = "pdf"
)

ggsave(
  glue("{result_folder}Zeisel_et_al-neuron_scRNA_Seq.png"),
  plot = umap,
  width = 10,
  height = 10,
  dpi = 300,
)

# export Seurat object
neuron_rna@meta.data = neuron_rna@meta.data %>% rename(cell_type = orig.ident)
saveRDS(neuron_rna,
        file.path(SCRIPT_DIR, "scRNA_Seq-Zeisel_et_al-neuron.Rds"))
write_tsv(norm_data,
          file.path(SCRIPT_DIR, "scRNA_Seq-Zeisel_et_al-neuron-norm_data.tsv"))
