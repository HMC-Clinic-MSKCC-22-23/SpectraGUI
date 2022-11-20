# SpectraGUI
A GUI to sit over top SPECTRA, a factor analysis tool built by MSKCC.

## How to download
- packages to install?

## How to Run

### Running the script
Run the python script SpectraGUI.py the way you normally run scripts

### When the window opens
1. AnnData file path - enter the file path to your annData set
2. AnnData key for cell type - enter the cell type key that's specified in your AnnData set. 
3. Pathway Annotations - if you click edit, it will open a new window. You can now either manually edit annotations or upload a csv file. If you upload a file, the first row must contain a column for "Cell Type", "Pathway Names" and "Genes". Double clicking on a row in the top box will fill the bottom box with the encoded dictionary of genes. Clicking cancel will cancel any changes, you must click "Save" to save changes.
4. Lambda Value - insert the lambda value
5. Advanced options - provides options for "Use highly variable," rho, delta, kappa, whether or not to use weights, and the top genes factor.
6. Click "Run" when you're ready!

### Output Screen
1. UMAP (changing factors, recoloring, etc)
2. Heatmap
3. gene-gene graph
4. save options
5. re-running the model
