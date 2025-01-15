library(Seurat)
library(ggplot2)
library(monocle3)
#library(dplyr)

#source('/home/6j9/projects/atopic_dermatitis/scripts/trajectories/fixed_plot-cells.R')


obj_path <- '/home/6j9/projects/atopic_dermatitis/data/seurat_objects/ko/imputed.rds'
cell_meta_path <- '/home/6j9/projects/atopic_dermatitis/data/seurat_objects/ko/metadata/cell_meta_imputed.tsv'
clusters_path <- '/home/6j9/projects/atopic_dermatitis/data/seurat_objects/ko/metadata/clusters_obj-type-preprocessed_resolution-0.3_num-pcs-50.tsv'
celltypes_path <- '/home/6j9/projects/atopic_dermatitis/data/seurat_objects/ko/metadata/ko_pcs-50_res-0.3_min-dist-0.3_nn-30_celltypes.tsv'
#celltypes_path <- '/home/6j9/projects/atopic_dermatitis/data/seurat_objects/ko/metadata/annotated_clusters_obj-type-preprocessed_resolution-0.3_num-pcs-50.tsv'



# read in seurat object
obj <- readRDS(obj_path)

# format counts for monocle
mtx <- obj[['RNA']]$data
mtx <- as.matrix(mtx)
rownames(mtx) <- rownames(obj)
colnames(mtx) <- colnames(obj)

# format gene meta for monocle
gene_meta <- as.data.frame(rownames(obj))
rownames(gene_meta) <- rownames(obj)
colnames(gene_meta) <- c('gene_short_name')

# format cell meta for monocle
cell_meta <- read.table(cell_meta_path, sep='\t', header=TRUE)

# create cell_data_set (cds) - monocle3 equivalent of seurat obj
cds <- new_cell_data_set(mtx, cell_meta, gene_meta)

# add clusters to cell metadata
clusters <- read.table(clusters_path, sep='\t', header=TRUE, row.names=1)
colData(cds) <- cbind(colData(cds), clusters)

# add cluster celltype labels to cell metadata
# by using cluster number to celltype label map
ct_map <- read.table(celltypes_path, sep='\t', header=TRUE)
cell_meta <- colData(cds)
cell_meta$ct_labels <- ct_map$celltype[match(cell_meta$seurat_clusters, ct_map$cluster)]
cell_meta$ct_labels <- as.factor(cell_meta$ct_labels)
colData(cds) <- cell_meta

# run UMAP and trajectory analysis
cds <- preprocess_cds(cds, method='PCA', num_dim=50)
cds <- reduce_dimension(cds, reduction_method='UMAP', preprocess_method='PCA')
cds <- cluster_cells(cds, reduction_method='UMAP', k=50)
cds <- learn_graph(cds)

# plot regular UMAP for sanity check of correct celltypes
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/ct_umap.png'
out_plot <- plot_cells(cds, reduction_method='UMAP', color_cells_by='ct_labels',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       #group_cells_by='ct_labels',
                       label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# plot by condition because useful
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/condition.png'
out_plot <- plot_cells(cds, reduction_method='UMAP', color_cells_by='condition',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       label_groups_by_cluster=FALSE, label_cell_groups=FALSE)
ggsave(out_path, out_plot)





# plot trajectory analysis
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/trajectories_pc-pts.png'
out_plot <- plot_cells(cds, reduction_method='UMAP',
                       label_cell_groups=FALSE,
                       group_label_size=3, show_trajectory_graph=TRUE,
                       label_principal_points=TRUE,
                       #group_cells_by='ct_labels',
                       label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/trajectories_less.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=TRUE,
                       #group_cells_by='ct_labels',
                       label_leaves=FALSE,
                       label_branch_points=FALSE,
                       label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# plotting certain genes to figure out trajectories
genes <- c('CCL11',
           'CEBPB',
           'CCL2',
           'COL6A5')
           #'MDK',
           #'LEPR',
           #'CXCL1',
           #'HAS2',
           #'TNFAIP6')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/genes_1.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       #group_cells_by='ct_labels',
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# order cells by principal nodes and plot by pseudotime
cds <- order_cells(cds, root_pr_nodes=c('Y_58', 'Y_97'))
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/pseudotime_col6a5-leukos.png'
out_plot <- plot_cells(cds, color_cells_by='pseudotime', reduction_method='UMAP',
                       label_cell_groups=FALSE, label_leaves=FALSE,
                       label_branch_points=FALSE)
ggsave(out_path, out_plot)

cds <- order_cells(cds, root_pr_nodes=c('Y_58'))
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/pseudotime_col6a5.png'
out_plot <- plot_cells(cds, color_cells_by='pseudotime', reduction_method='UMAP',
                       label_cell_groups=FALSE, label_leaves=FALSE,
                       label_branch_points=FALSE)
ggsave(out_path, out_plot)












