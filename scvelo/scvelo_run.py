import scvelo as scv

scv.settings.verbocity = 3
scv.settings.presenter_view = True

crcdata = scv.read('crctest.loom', cache=True)
crcdata.var_names_make_unique()
crcdata

# scv.pl.proportions(adata) # view proportion spliced to unspliced

scv.pp.filter_and_normalize(crcdata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(crcdata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(crcdata)
scv.tl.velocity_graph(crcdata)

scv.tl.umap(crcdata) # set to umap
scv.pl.velocity_embedding_stream(adata, basis='umap') # plot velocity 
scv.pl.velocity(crcdata, ['CD3D',  'CD68', 'DCN', 'EPCAM', 'KIT', 'CD79A'], ncols=2) # view marker genes