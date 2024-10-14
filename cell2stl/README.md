# cell2stl

This subfolder should contain a C++ or python program to take 1. a cell type and produce a crude 3d model of what that cell type looks like for rendering on the webstie and 2. some number of cell surface proteins and transcription factors and create a 3d model of what the cell surface would look like, with weighted 3d arrows that point to and from any relevant cs proteins and tf proteins. 

# Order of To-dos (hopefully broken down so you can downsize the project if necessary)
1. Simply cell type -> .stl 3d model of that cell type. This should be more simple, we can start by only producing the 10 most populous cell types from tabula sapiens.
2. Use PDB to take known cell surface protein and transcription factors and get the .stl 3d model for those proteins. Then, try to render it such that those proteins are embedded in a reasonable approximation of what a cell surface and internal part of a cell would look like. Finally, for extra super duper bonus points, make zooming out dynamically switch rendering view from a zoomed in high resolution of those proteins to an approximation of the proportion of cell surface and internal part made up by those proteins, with some color coding. Start with just 5 or 10 proteins each and go more complex from there.
3. Weighted 3d arrows while zoomed in. 
