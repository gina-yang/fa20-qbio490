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

### ANNDATA UTILITIES ### https://anndata.readthedocs.io/en/latest/api.html
# adata = adata1.concatenate(adata2) # to concatenate two data objects (adata = adata1 + adata2)
# adata = tumorborder.concatenate(tumorcore, batch_key='batch',batch_categories=['tumor border', 'tumor core']) 
# adata.obs # shows observations
# adata.var # shows cell info
# adata.obs[name] # show obs called 'name'

### WRITING DATA MATRIX TO FILE ###
# adata.X # this is the expression matrix of n_obs x n_vars. it is a scipy/numpy sparse matrix
# kul01.to_df().to_csv('myfile.csv') # to save expression matrix to csv file (reference https://github.com/theislab/scanpy/issues/262)


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
	
# Plot velocity as streamlines with cells colored by cluster (found by louvain)
scv.pl.velocity_embedding_stream(data_obj, color='louvain')
# This time colored by batch
scv.pl.velocity_embedding_stream(adata, color='batch', legend_loc='lower left') 


# Plot phase portraits of marker genes
scv.pl.velocity(data_obj, ['CD3D', 'CD68', 'DCN', 'EPCAM', 'KIT', 'CD79A'], ncols=2)


# def get_top_genes(data_obj, csvfile):
# 	# sc.tl.rank_genes_groups() to obtain cluster markers
# 	scv.tl.rank_velocity_genes(data_obj, groupby='louvain', min_corr=.3)
# 	df = scv.DataFrame(data_obj.uns['rank_velocity_genes']['names']) 
# 	#  df.head() # See top 5 genes for each cluster
# 	df.to_csv(csvfile, index=False)

