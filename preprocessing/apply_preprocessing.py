import pandas as pd
import numpy as np
import scipy.io
import os, re, time

# settings
lower_umi = 500
upper_umi = 25000
mouse_mito = 0.05
human_mito = 0.10
species_map = {'he': 'human',
               'hughes': 'human',
               'ko': 'mouse',
               'leyva-castillo': 'mouse',
               'reynolds': 'human'}

def apply_preprocessing(mtx, species):
    umi_counts = np.sum(mtx.values, axis=0)
    umi_filter = umi_counts > lower_umi
    mtx = mtx.loc[:, umi_filter]
    umi_counts = np.sum(mtx.values, axis=0)
    umi_filter = umi_counts < upper_umi
    mtx = mtx.loc[:, umi_filter]
    umi_counts = np.sum(mtx.values, axis=0)
    mito_idx = np.array([True if re.match('^MT-', gene) else False for gene in  mtx.index])
    mito_counts = np.sum(mtx.values[mito_idx, :], axis=0)
    mito_proportion = mito_counts / umi_counts
    if species == 'mouse':
        mito_filter = mito_proportion < mouse_mito
    elif species == 'human':
        mito_filter = mito_proportion < human_mito
    mtx = mtx.loc[:, mito_filter]
    return mtx

def read_mtx(sample_dir, sample_name):
    mtx = scipy.io.mmread(os.path.join(sample_dir, 'matrix.mtx'))
    bcs = np.squeeze(pd.read_csv(os.path.join(sample_dir, 'barcodes.tsv'), header=None, index_col=None).values)
    bcs = np.array([bc + '_' + sample_name for bc in bcs])
    genes = np.squeeze(pd.read_csv(os.path.join(sample_dir, 'features.tsv'), sep = '\t', header=None, index_col=None).iloc[:, 1].values)
    genes = [gene.upper() for gene in genes]
    uniq_genes, counts = np.unique(genes, return_counts=True)
    dup_genes = uniq_genes[counts > 1]
    for gene in dup_genes:
        gene_idxs = np.arange(len(genes))[np.array([True if x == gene else False for x in genes])]
        for i, idx in enumerate(gene_idxs):
            genes[idx] = genes[idx] + f'_{i+1}'
    mtx = pd.DataFrame.sparse.from_spmatrix(mtx)
    mtx.index = genes
    mtx.columns = bcs
    return mtx

def write_mtx(out_dir, mtx, sample_name):
    sample_out_dir = os.path.join(out_dir, sample_name)
    if not os.path.isdir(sample_out_dir):
        os.mkdir(sample_out_dir)
    np.savetxt(os.path.join(sample_out_dir, 'barcodes.tsv'), mtx.columns.values, fmt='%s')
    features = np.stack((mtx.index, mtx.index, mtx.index)).T
    np.savetxt(os.path.join(sample_out_dir, 'features.tsv'), features, delimiter='\t', fmt='%s')
    scipy.io.mmwrite(os.path.join(sample_out_dir, 'matrix.mtx'), scipy.sparse.csr_matrix(mtx.values))

# paths
root_in_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/human_mouse_skin/data/bam'
root_out_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/human_mouse_skin/data/count_matrices'

author_dirs = os.listdir(root_in_dir)
author_dirs = [os.path.join(root_in_dir, path) for path in author_dirs if 'archive' not in path]

for in_dir in author_dirs:
    _, name = os.path.split(in_dir)
    print(name)
    out_dir = os.path.join(root_out_dir, name)
    data_dirs = []
    for r, d, f in os.walk(in_dir):
        for path in f:
            if 'matrix.mtx' in path and 'filtered' in r:
                data_dirs.append(r)
    data_dirs.sort()
    # apply preprocessing per-author and per-sample
    for sample_dir in data_dirs:
        sample_name = re.search('\/([a-zA-Z0-9\_]+)_Solo.out', sample_dir).groups(1)[0]
        mtx = read_mtx(sample_dir, sample_name)
        mtx = apply_preprocessing(mtx, 'mouse')
        write_mtx(out_dir, mtx, sample_name)
