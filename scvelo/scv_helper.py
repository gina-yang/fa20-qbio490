import scvelo as scv

# Perform setup and data loading
# Returns object that stores data matrix with observations and variables
def setup(loomfile):
	scv.settings.verbocity = 3
	scv.settings.presenter_view = True
	scv.set_figure_params('scvelo') # For beautified visualization

	crcdata = scv.read(loomfile, cache=True)
	crcdata.var_names_make_unique()
	return(crcdata)


# crcdata # shows dimensions: n_obs x n_vars (n_obs=number of cells, n_vars=cell info)
# adata = adata1.concatenate(adata2) # to concatenate two data objects (adata = adata1 + adata2)


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

# Estimates RNA velocity and projects and visualizes them
def plot_velocity(data_obj, louvain_res):
	scv.tl.velocity(data_obj)
	scv.tl.velocity_graph(data_obj)
	scv.tl.louvain(data_obj, resolution=louvain_res) # resolution default = 1, try 0.5
	# Set to UMAP embedding
	scv.tl.umap(data_obj)
	# Plot velocity as streamlines
	scv.pl.velocity_embedding_stream(data_obj, basis='umap')

# Pass marker genes in list to plot phase portraits
# marker_genes = ['CD3D', 'CD68', 'DCN', 'EPCAM', 'KIT', 'CD79A']

# Plot phase portraits of marker genes
def plot_genes(data_obj, gene_list):
	scv.pl.velocity(data_obj, gene_list, ncols=2)


# Output high expressed genes for each louvain cluster to a csv file
def get_top_genes(data_obj, csvfile):
	scv.tl.rank_velocity_genes(data_obj, groupby='louvain', min_corr=.3)
	df = scv.DataFrame(data_obj.uns['rank_velocity_genes']['names']) 
	#  df.head() # See top 5 genes for each cluster
	df.to_csv(csvfile, index=False)

