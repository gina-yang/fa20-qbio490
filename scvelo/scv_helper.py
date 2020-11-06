import scvelo as scv

# Perform setup and data loading
# Returns object that stores data matrix with observations and variables
def setup(loomfile):
	scv.settings.verbocity = 3
	scv.settings.presenter_view = True
	# scv.set_figure_params('scvelo') # For beautified visualization

	crcdata = scv.read(loomfile, cache=True)
	crcdata.var_names_make_unique()
	return(crcdata)



# Preprocessing
	# Filter genes based on number of counts; extract highly variable genes; normalize each 
	# cell by total counts over all genes; logarithmize the data matrix
		# min_shared_counts: minimum number of counts (both spliced and unspliced) required for a gene
		# n_top_genes: number of highly variable genes to keep
def preproc(data_obj):
	scv.pp.filter_and_normalize(data_obj, min_shared_counts=20, n_top_genes=2000)
	# Computes moments for velocity estimation with pca and computing neighbors
	scv.pp.moments(data_obj, n_pcs=30, n_neighbors=30)
	return(data_obj)


# Estimates RNA velocity
def compute_velocity(data_obj, louvain_res):
	scv.tl.velocity(data_obj)
	scv.tl.velocity_graph(data_obj)
	scv.tl.louvain(data_obj, resolution=louvain_res) # resolution default = 1, try 0.5
	# Set to UMAP embedding
	scv.tl.umap(data_obj)
	

