import matplotlib
matplotlib.use('Agg')
from spring_helper import *
from doublet_detector import *
from collections import defaultdict
import csv

# Adapted from SPRING example at https://github.com/AllonKleinLab/SPRING_dev/blob/master/data_prep/spring_example_HPCs.ipynb



#Load in Data (Human)
sample_name = ['GSE102580_filtered_normalized_counts_human_transposed']

min_tot = [1 for s in sample_name] # initial guess for total transcript counts threshold
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
    print D[s]['E'].shape

gene_list = np.array(load_genes(input_path + 'genes_human.txt'))







#####filter by total counts
# adjust total counts thresholds
D['GSE102580_filtered_normalized_counts_human_transposed']['meta']['min_tot'] = 1


for s in sample_name:
    D[s]['total_counts'] = np.sum(D[s]['E'], axis=1).A[:,0]

    ix = D[s]['total_counts'] >= D[s]['meta']['min_tot']
    print s, np.sum(ix), '/', D[s]['E'].shape[0], np.median(D[s]['total_counts'][ix]), np.mean(D[s]['total_counts'][ix])

# Actually filter out low-count barcodes

for s in sample_name:
    print '---  %s ---' %s
    print 'Pre-filter: %i barcodes' %D[s]['E'].shape[0]
    D[s]['cell_index'] = np.arange(D[s]['E'].shape[0])
    tmpfilt = np.nonzero(D[s]['total_counts'] >= D[s]['meta']['min_tot'])[0]
    D[s] = filter_dict(D[s], tmpfilt)
    print 'Post-filter: %i barcodes' %D[s]['E'].shape[0]

del tmpfilt








# create master dataset (all SPRING subsets will refer back to this)

# build sample index
samp_id_flat = np.array([],dtype=str)
for s in sample_name:
    samp_id_flat = np.append(samp_id_flat, [s] * D[s]['E'].shape[0])

# merge quality metrics
total_counts = np.zeros(len(samp_id_flat), dtype=int)
for s in sample_name:
    total_counts[samp_id_flat == s] = D[s]['total_counts']

# merge counts matrices
E = scipy.sparse.vstack([D[s]['E'] for s in sample_name]).tocsc()









##Save base directory files
# Set path for saving data -- you'll have to change this for your own setup.
# This path should be a subdirectory of your local copy of SPRING,
# specifically, {path_to_SPRING}/datasets/{main_dataset_name}.
# See example below, where springViewer_1_6_dev.html is located in ../

main_spring_dir = '../datasets/Plasschaert_human_homeostasis/'
print main_spring_dir

if not os.path.exists(main_spring_dir):
    os.makedirs(main_spring_dir)

np.savetxt(main_spring_dir + 'genes.txt', gene_list, fmt='%s')
np.savetxt(main_spring_dir + 'total_counts.txt', total_counts)







##Save SPRING files
save_path = main_spring_dir + 'full_k=4'
print save_path
spec_clusters = []
with open(save_path + '/SpectralClusters.csv') as inputfile:
    for row in csv.reader(inputfile):
        spec_clusters.append(row[0])
cell_groupings= {'cluster':spec_clusters}
dict.items(cell_groupings)
out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 1, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 4,
                    cell_groupings = cell_groupings,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')

##Save SPRING files
save_path = main_spring_dir + 'full_k=5'
print save_path
spec_clusters = []
with open(save_path + '/SpectralClusters.csv') as inputfile:
    for row in csv.reader(inputfile):
        spec_clusters.append(row[0])
cell_groupings= {'cluster':spec_clusters}
dict.items(cell_groupings)
out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 1, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 5,
                    cell_groupings = cell_groupings,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')

##Save SPRING files
save_path = main_spring_dir + 'full_k=6'
print save_path
spec_clusters = []
with open(save_path + '/SpectralClusters.csv') as inputfile:
    for row in csv.reader(inputfile):
        spec_clusters.append(row[0])
cell_groupings= {'cluster':spec_clusters}
dict.items(cell_groupings)
out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 1, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 6,
                    cell_groupings = cell_groupings,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')

##Save SPRING files
save_path = main_spring_dir + 'full_k=7'
print save_path
spec_clusters = []
with open(save_path + '/SpectralClusters.csv') as inputfile:
    for row in csv.reader(inputfile):
        spec_clusters.append(row[0])
cell_groupings= {'cluster':spec_clusters}
dict.items(cell_groupings)
out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 1, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 7,
                    cell_groupings = cell_groupings,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')

##Save SPRING files
save_path = main_spring_dir + 'full_k=8'
print save_path
spec_clusters = []
with open(save_path + '/SpectralClusters.csv') as inputfile:
    for row in csv.reader(inputfile):
        spec_clusters.append(row[0])
cell_groupings= {'cluster':spec_clusters}
dict.items(cell_groupings)
out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 1, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 8,
                    cell_groupings = cell_groupings,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')

##Save SPRING files
save_path = main_spring_dir + 'full_k=9'
print save_path
spec_clusters = []
with open(save_path + '/SpectralClusters.csv') as inputfile:
    for row in csv.reader(inputfile):
        spec_clusters.append(row[0])
cell_groupings= {'cluster':spec_clusters}
dict.items(cell_groupings)
out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 1, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 9,
                    cell_groupings = cell_groupings,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')

##Save SPRING files
save_path = main_spring_dir + 'full_k=10'
print save_path
spec_clusters = []
with open(save_path + '/SpectralClusters.csv') as inputfile:
    for row in csv.reader(inputfile):
        spec_clusters.append(row[0])
cell_groupings= {'cluster':spec_clusters}
dict.items(cell_groupings)
out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 1, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 10,
                    cell_groupings = cell_groupings,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')

##Save SPRING files
save_path = main_spring_dir + 'full_k=11'
print save_path
spec_clusters = []
with open(save_path + '/SpectralClusters.csv') as inputfile:
    for row in csv.reader(inputfile):
        spec_clusters.append(row[0])
cell_groupings= {'cluster':spec_clusters}
dict.items(cell_groupings)
out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 1, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 11,
                    cell_groupings = cell_groupings,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')


##Save SPRING files
save_path = main_spring_dir + 'full_k=12'
print save_path
spec_clusters = []
with open(save_path + '/SpectralClusters.csv') as inputfile:
    for row in csv.reader(inputfile):
        spec_clusters.append(row[0])
cell_groupings= {'cluster':spec_clusters}
dict.items(cell_groupings)
out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 1, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 12,
                    cell_groupings = cell_groupings,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')
