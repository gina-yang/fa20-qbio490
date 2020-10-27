import scvelo as scv

def setup(loomfile):
	scv.settings.verbocity = 3
	scv.settings.presenter_view = True

	crcdata = scv.read(loomfile, cache=True)
	crcdata.var_names_make_unique()
	crcdata # shows dimensions: n_obs x n_vars (n_obs=number of cells, n_vars=cell info)

	# scv.pl.proportions(adata) # View proportion spliced/unspliced counts (in pie chart)

	# Preprocessing
	# Filter genes based on number of counts; extract highly variable genes; normalize each 
	# cell by total counts over all genes; logarithmize the data matrix
		# min_shared_counts: minimum number of counts (both spliced and unspliced) required for a gene
		# n_top_genes: number of highly variable genes to okeep
	scv.pp.filter_and_normalize(crcdata, min_shared_counts=20, n_top_genes=2000)
	# Computes moments for velocity estimation with pca and computing neighbors
	scv.pp.moments(crcdata, n_pcs=30, n_neighbors=30)

	scv.tl.velocity(crcdata)
	scv.tl.velocity_graph(crcdata)
	return(crcdata)


def plot_velocity(data_obj, louvain_res):
	scv.tl.louvain(data_obj, resolution=louvain_res) # resolution default = 1, try 0.5
	# Set to UMAP embedding
	scv.tl.umap(data_obj)
	# Plot velocity as streamlines
	scv.pl.velocity_embedding_stream(data_obj, basis='umap')
	# Pass marker genes in list to plot expression/velocity

# marker_genes = ['CD3D',  'CD68', 'DCN', 'EPCAM', 'KIT', 'CD79A']

def plot_genes(data_obj, gene_list):
	scv.pl.velocity(data_obj, gene_list, ncols=2)
