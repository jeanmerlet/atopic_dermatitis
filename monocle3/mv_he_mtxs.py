import subprocess
import os


data_dir = '/home/6j9/projects/atopic_dermatitis/data/bam/he'
data_dirs = [os.path.join(data_dir, data_subdir) for data_subdir in os.listdir(data_dir)]
out_dir = '/home/6j9/projects/atopic_dermatitis/data/count_matrices/he'

for data_dir in data_dirs:
    _, dirname = os.path.split(data_dir)
    dirname = dirname.split('_')[0]
    out_mtx_dir = os.path.join(out_dir, dirname)
    os.makedirs(out_mtx_dir, exist_ok=True)
    subprocess.run(f'mv {data_dir}/Gene/filtered/* {out_mtx_dir}', shell=True)
