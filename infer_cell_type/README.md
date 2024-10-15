# infer_cell_type

This subfolders job is to contain the code necessary for 1. inferring cell type from gene expression matrix 2. Determining whether a genes are encoding transcription factors or cell surface receptors to send to cell2stl to render models of those specific gene types 3. (optional) some sensible optimizations using caching etc. so a user could do minor updates to the gene expression matrix and see if that changes the inferred cell type or encoded TFs/cs receptors

# TO-DOs
1. Infer cell type from a gene expression matrix (just use python for this)
2. Determine whether genes are encoding TFs or cell surface receptors (this might be best to use an SQL table lookup)
3. Need some sensible interface for cell2stl to use. I'm thinking users can have an SQL table with info and cell2stl can do a lookup for the transcription factors and cellsurface receptors it needs, + inferred cell type