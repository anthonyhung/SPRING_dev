from spring_helper import *
from doublet_detector import *
from collections import defaultdict



sample_name = ['P9A', 'P11A', 'P11B', 'P12A']

min_tot = [1000 for s in sample_name] # initial guess for total transcript counts threshold
nSamp = len(sample_name)
input_path = 'raw_counts/'


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
        D[s]['E'] = load_text(file_opener(input_path + s + '.counts.tsv.gz'), delim = '\t')
        scipy.sparse.save_npz(input_path + s + '.raw_counts.unfiltered.npz', D[s]['E'], compressed = True)
    print D[s]['E'].shape


gene_list = np.array(load_genes(input_path + 'genes.txt'))


# plot total counts histograms - don't actually filter out any barcodes yet


# adjust total counts thresholds
D['P9A']['meta']['min_tot'] = 900
D['P11A']['meta']['min_tot'] = 700
D['P11B']['meta']['min_tot'] = 800
D['P12A']['meta']['min_tot'] = 800

for s in sample_name:
    D[s]['total_counts'] = np.sum(D[s]['E'], axis=1).A[:,0]
    ix = D[s]['total_counts'] >= D[s]['meta']['min_tot']



# Actually filter out low-count barcodes

for s in sample_name:
    print '---  %s ---' %s
    print 'Pre-filter: %i barcodes' %D[s]['E'].shape[0]
    D[s]['cell_index'] = np.arange(D[s]['E'].shape[0])
    tmpfilt = np.nonzero(D[s]['total_counts'] >= D[s]['meta']['min_tot'])[0]
    D[s] = filter_dict(D[s], tmpfilt)
    print 'Post-filter: %i barcodes' %D[s]['E'].shape[0]

del tmpfilt





# get mitochondrial genes

mt_ix = [i for i,g in enumerate(gene_list) if g.startswith('mt-')]
print [gene_list[i] for i in mt_ix]



# plot mito-gene frac histograms - don't actually filter out any cells yet

# set mito-gene frac threshold
for s in sample_name:
    D[s]['meta']['max_mt'] = 0.15

    D[s]['mito_frac'] = np.sum(D[s]['E'][:,mt_ix], axis=1).A[:,0] / np.sum(D[s]['E'], axis=1,dtype=float).A[:,0]

# Actually filter out mito-high cells

for s in sample_name:
    print '---  %s ---' %s
    print 'Pre-filter: %i barcodes' %D[s]['E'].shape[0]
    tmpfilt = np.nonzero(D[s]['mito_frac'] <= D[s]['meta']['max_mt'])[0]
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
mito_frac = np.zeros(len(samp_id_flat), dtype=float)
for s in sample_name:
    total_counts[samp_id_flat == s] = D[s]['total_counts']
    mito_frac[samp_id_flat == s] = D[s]['mito_frac']

# merge counts matrices
E = scipy.sparse.vstack([D[s]['E'] for s in sample_name]).tocsc()






# normalize by total counts
E = tot_counts_norm(E)[0]





# Set path for saving data -- you'll have to change this for your own setup.
# This path should be a subdirectory of your local copy of SPRING,
# specifically, {path_to_SPRING}/datasets/{main_dataset_name}.
# See example below, where springViewer_1_6_dev.html is located in ../

main_spring_dir = '../datasets/hpc/'

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






cellCycleGenes_G2M = ['Ube2c','Hmgb2','Hmgn2','Tuba1b','Mki67','Ccnb1','Tubb','Top2a','Tubb4b']



for s in sample_name:
    t0 = time.time()
    print '________________', s

    cell_ix = samp_id_flat == s
    save_path = main_spring_dir + s

    out = make_spring_subplot(E[cell_ix,:], gene_list, save_path,
                        normalize = False, tot_counts_final = total_counts[cell_ix],
                        min_counts = 3, min_cells = 3, min_vscore_pctl = 85,show_vscore_plot = True, 
                        num_pc = 30,
                        k_neigh = 4,
                        num_force_iter = 200)

    np.save(save_path + '/cell_filter.npy', np.nonzero(cell_ix)[0])
    np.savetxt(save_path + '/cell_filter.txt', np.nonzero(cell_ix)[0], fmt='%i')

    print time.time() - t0
