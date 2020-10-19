import scvelo as scv

scv.settings.verbocity = 3
scv.settings.presenter_view = True

crcdata = scv.read('kul01.loom', cache=True)
crcdata.var_names_make_unique()
crcdata

# scv.pl.proportions(adata) # View proportion spliced to unspliced (pie chart)

# ??
scv.pp.filter_and_normalize(crcdata, min_shared_counts=20, n_top_genes=2000)
# 
scv.pp.moments(crcdata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(crcdata)
scv.tl.velocity_graph(crcdata)

# Set to UMAP embedding
scv.tl.umap(crcdata)
# Plot velocity
scv.pl.velocity_embedding_stream(crcdata, basis='umap')
# Pass marker genes in list to plot expression/velocity
scv.pl.velocity(crcdata, ['CD3D',  'CD68', 'DCN', 'EPCAM', 'KIT', 'CD79A'], ncols=2)