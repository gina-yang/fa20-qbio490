import scvelo as scv

scv.settings.verbocity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo') # For beautified visualization

crcdata = scv.read('kul19_normal.loom', cache=True)
crcdata.var_names_make_unique()
crcdata # shows dimensions: n_obs x n_vars (n_obs=number of cells, n_vars=cell info)

# To concatenate data objects (adata3 = adata2 + adata1) : adata3 = adata2.concatenate(adata1)

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



# Plot phase portraits for marker genes
scv.pl.velocity(crcdata, ['CD3D',  'CD68', 'DCN', 'EPCAM', 'KIT', 'CD79A'], ncols=2)

# Identify hihgly expressed genes in each cluster
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)
