from spring_helper import *
from doublet_detector import *
from collections import defaultdict
# Adapted from SPRING example at https://github.com/AllonKleinLab/SPRING_dev/blob/master/data_prep/spring_example_HPCs.ipynb

#setting parameters for presentation items
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42









#Load in Data (Human)
sample_name = ['GSE102580_filtered_normalized_counts_human']

min_tot = [1000 for s in sample_name] # initial guess for total transcript counts threshold
nSamp = len(sample_name)
input_path = 'filtered_normalized_counts/'


# D stores all the data; one entry per library
D = {}

for j,s in enumerate(sample_name):
    D[s] = {}
    D[s]['meta'] = {'min_tot': min_tot[j]}

# load counts matrices -- try to load from npz file (fast),
# otherwise load from text file (slow) and build npz file for next time

for s in sample_name:
    print '_________________', s

    if os.path.isfile(input_path + s + '.raw_counts.unfiltered.npz'):
        print 'Loading from npz file'
        D[s]['E'] = scipy.sparse.load_npz(input_path + s + '.raw_counts.unfiltered.npz')
    else:
        print 'Loading from text file'
        D[s]['E'] = load_text(file_opener(input_path + s + '.tsv.gz'), delim = '\t')
        scipy.sparse.save_npz(input_path + s + '.filtered_normalized_counts.npz', D[s]['E'], compressed = True)
    print D[s]['E'].shape

gene_list = np.array(load_genes(input_path + 'genes.txt'))









# create master dataset (all SPRING subsets will refer back to this)

# build sample index
samp_id_flat = np.array([],dtype=str)
for s in sample_name:
    samp_id_flat = np.append(samp_id_flat, [s] * D[s]['E'].shape[0])

# merge quality metrics
total_counts = np.zeros(len(samp_id_flat), dtype=int)
mito_frac = np.zeros(len(samp_id_flat), dtype=float)
for s in sample_name:
    total_counts[samp_id_flat == s] = D[s]['total_counts']
    mito_frac[samp_id_flat == s] = D[s]['mito_frac']

# merge counts matrices
E = scipy.sparse.vstack([D[s]['E'] for s in sample_name]).tocsc()









##Save base directory files
# Set path for saving data -- you'll have to change this for your own setup.
# This path should be a subdirectory of your local copy of SPRING,
# specifically, {path_to_SPRING}/datasets/{main_dataset_name}.
# See example below, where springViewer_1_6_dev.html is located in ../

main_spring_dir = '../datasets/Plasschaert_human_homeostasis/'

if not os.path.exists(main_spring_dir):
    os.makedirs(main_spring_dir)

np.savetxt(main_spring_dir + 'genes.txt', gene_list, fmt='%s')
np.savetxt(main_spring_dir + 'total_counts.txt', total_counts)

# save master expression matrices

print 'Saving hdf5 file for fast gene loading...'
save_hdf5_genes(E, gene_list, main_spring_dir + 'counts_norm_sparse_genes.hdf5')

##############
print 'Saving hdf5 file for fast cell loading...'
save_hdf5_cells(E, main_spring_dir + 'counts_norm_sparse_cells.hdf5')

##############
save_sparse_npz(E, main_spring_dir + 'counts_norm.npz', compressed = False)









##
