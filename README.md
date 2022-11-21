# SpectraGUI
A GUI to sit over top SPECTRA, a factor analysis tool built by MSKCC.

## How to download
Assuming you already have SPECTRA, download this and place in the same folder as SPECTRA.

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
1. UMAP - can select factors to recolor the map, to see changes, click the recolor button.
2. Heatmap - can select aspects of the adata to highlight, choose from the dropdown.
3. gene-gene graph - click the button and a separate window will appear that you can resize in order to see the whole graph better.
4. Save options - check the things you want to save and click save, and it will download in the same folder as the GUI
5. Re-running the model - click re-run to go back to the input screen. It will keep your parameters, but will lose your trained model. That functionality hasn't been implemented yet.
