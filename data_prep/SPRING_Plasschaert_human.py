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
