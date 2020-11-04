import scvelo as scv

scv.settings.verbocity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo') # For beautified visualization

crcdata = scv.read('kul19_normal.loom', cache=True)
crcdata.var_names_make_unique()
crcdata # shows dimensions: n_obs x n_vars (n_obs=number of cells, n_vars=cell info)

### ANNDATA UTILITIES ### https://anndata.readthedocs.io/en/latest/api.html
# adata = adata1.concatenate(adata2) # to concatenate two data objects (adata = adata1 + adata2)
# adata = tumorborder.concatenate(tumorcore, batch_key='batch',batch_categories=['tumor border', 'tumor core']) 
# adata.obs # shows observations
# adata.var # shows cell info
# adata.obs[name] # show obs called 'name'

### WRITING DATA MATRIX TO FILE ###
# adata.X # this is the expression matrix of n_obs x n_vars. it is a scipy/numpy sparse matrix
# kul01.to_df().to_csv('myfile.csv') # to save expression matrix to csv file (reference https://github.com/theislab/scanpy/issues/262)

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


scv.tl.louvain(crcdata, resolution=0.5) # resolution default=1, try 0.5

# Set to UMAP embedding
scv.tl.umap(crcdata)
# Plot velocity as streamlines
scv.pl.velocity_embedding_stream(crcdata, basis='umap')


# Plot velocity as streamlines with cells colored by cluster (found by louvain)
scv.pl.velocity_embedding_stream(crcdata, color='louvain')
# This time colored by batch
scv.pl.velocity_embedding_stream(crcdata, color='batch', legend_loc='lower left') 


# Plot phase portraits of marker genes
scv.pl.velocity(crcdata, ['CD3D', 'CD68', 'DCN', 'EPCAM', 'KIT', 'CD79A'], ncols=2)


# def get_top_genes(data_obj, csvfile):
# 	# sc.tl.rank_genes_groups() to obtain cluster markers
# 	scv.tl.rank_velocity_genes(data_obj, groupby='louvain', min_corr=.3)
# 	df = scv.DataFrame(data_obj.uns['rank_velocity_genes']['names']) 
# 	#  df.head() # See top 5 genes for each cluster
# 	df.to_csv(csvfile, index=False)