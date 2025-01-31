library(Seurat)
library(ggplot2)
library(monocle3)


obj_path <- '/home/6j9/projects/atopic_dermatitis/data/seurat_objects/ko/imputed.rds'
cell_meta_path <- '/home/6j9/projects/atopic_dermatitis/data/seurat_objects/ko/metadata/cell_meta_imputed.tsv'
clusters_path <- '/home/6j9/projects/atopic_dermatitis/data/seurat_objects/ko/metadata/clusters_obj-type-preprocessed_resolution-0.3_num-pcs-50.tsv'
celltypes_path <- '/home/6j9/projects/atopic_dermatitis/data/seurat_objects/ko/metadata/ko_pcs-50_res-0.3_min-dist-0.3_nn-30_celltypes.tsv'



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
           'PRRX1',
           'COL6A5')
           #'MDK',
           #'LEPR',
           #'CXCL1',
           #'HAS2',
           #'TNFAIP6')
genes <- c('COL6A5')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/col6a5-traj.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=TRUE,
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


### fibroblast subtyping ###

# pan-fibroblasts #
genes <- c('PDGFRA',
           'DPT',
           'COL1A2',
           'COL3A1',
           'TWIST2',
           'VIM')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/pan-fibro.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# putative common dermal fibro-adipogenic progenitor #
genes <- c('EN1',
           'DLK1',
           'HIC1')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/fibro-adipo-prog.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# papillary dermal fibroblasts #
genes <- c('DPP4',
           'DLK1',
           'ATXN1', #SCA1
           'LY6A',
           'PRDM1',
           'EPHB2',
           'LRIG1',
           'TRPS1')
genes <- c('OL6A5',
           'COL23A1',
           'HSPB3', #SCA1
           'APCDD1',
           'WIF1',
           'NTN1',
           'PDPN')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/papillary-fibro.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# reticular dermal fibroblasts #
genes <- c('DPP4',
           'DLK1',
           'ATXN1',
           'LY6A')
genes <- c('THY1',
           'FMO1',
           'MYOC',
           'LSP1',
           'MGP',
           'ACTA2',
           'PPARG',
           'CD36')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/reticular-fibro.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# dermal papilla fibroblasts #
genes <- c('SOX2',
           'LEF1',
           'CRABP1',
           'RSPO3',
           'CORIN',
           'ALPL',
           'VCAN')
genes <- c('APCDD1',
           'AXIN2',
           'COLEC12',
           'PTGDS',
           'COL18A1',
           'SFRP2',
           'SOX18')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/dermal-papilla-fibro.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# myofibroblasts #
genes <- c('ACTA1',
           'TAGLN',
           'FN1')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/myofibroblasts.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# mesenchymal fibroblasts #
genes <- c('ASPN',
           'POSTN',
           'GPC3',
           'TNN',
           'SFRP1')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/mesenchymal-fibro.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# secretory reticular #
genes <- c('WISP2',
           'SLPI',
           'CTHRC1',
           'MFAP5',
           'TSPAN8')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/secretory-reticular.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# secretory-papillary #
genes <- c('APCDD1',
           'ID1',
           'WIF1',
           'COL18A1',
           'PTGDS')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/secretory-papillary.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# dermal fibroblasts #
genes <- c('COL11A1',
           'MYL4',
           'HES1',
           'CTNNB1',
           'SFRP2',
           'CRABP1',
           'FMO2',
           'PRG4')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/dermal-fibroblasts.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# epithelial vs endothelial

genes <- c('KRT8',
           'KRT18',
           'EPCAM',
           'CDH1',
           'LGR5',
           'PECAM1',
           'CDH5',
           'VWF',
           'NOS3')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/epi-endo.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# neutrophils

genes <- c('FCER1G',
           'S100A9',
           'S100A8',
           'LYZ2')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/neutrophils.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# inflammatory fibroblasts

genes <- c('COL6A5',
           'POSTN',
           'FN1',
           'IL4RA')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/inf-fibro.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# macrophages

genes <- c('CD68',
           'ADGRE1',
           'LY6C1',
           'MARCO',
           'CSF1R',
           'CCR2',
           'CX3CR1',
           'FCER2',
           'MRC1',
           'CXCL9',
           'ITGAM')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/macrophages.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# m1 macrophages #
genes <- c('NOS2',
           'TNF',
           'IL1B',
           'IL6',
           'IL12A')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/m1-macrophages.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)

# m2 macrophages #
genes <- c('MRC1',
           'CD163',
           'ARG1',
           'CCL17',
           'TGFB1')
out_path <- '/home/6j9/projects/atopic_dermatitis/plots/trajectories/fibroblast_subtypes/m2-macrophages.png'
out_plot <- plot_cells(cds, color_cells_by='ct_labels', reduction_method='UMAP',
                       group_label_size=3, show_trajectory_graph=FALSE,
                       genes=genes,
                       label_cell_groups=FALSE, label_groups_by_cluster=FALSE)
ggsave(out_path, out_plot)
