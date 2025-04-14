# Overview
Visualisation codes for Xenium datasets

## Visualisation of single cells
This code allows visualisation of *N* random cells which will pass some criteria. It can be cell area, number of transcripts, cell type etc. Basically anything, what sits (or can be added) to `adata.obs`. Each cell image contains crop with DAPI image, cell label polygon and transcripts (please note at the moment I only select top 5 gene types in all regions, or list of!) For an example check **visualise_cells_example.ipynb**
**plot_cells_by_par** function has next inputs (*mandatory inputs):
 - sdata*: spatial data object
 - n_hor*: number of columns
 - n_vert*: number of rows. In total number of cells *N = n_hor x n_vert*
 - col_name*: name of the column in 'adata.obs' to perform filtering on
 - min_value*: either minimum value for float columns, or string for categorical columns
 - max_value*: maqximum value for float columns
 - transcripts_all*: dataframe with transcripts positions and types (TODO: it can be read directly from sdata!)
 - pixelsize*: pixel size, can be found in image metadata (see Jupyter Notebook)
 - border_size_px: padding size around borders of the cell segmented polygon
 - label_col: color of the cell labelled polygon
 - label_width: width of the cell labelled polygon
 - save: path to the image be saved
 - spot_size: size of the spot transcripts
 - text_size: size of the legend text for the gene names
 - display_genes: list of teh gene names to be shown (for transcripts). If the list is empty (default option) shows top5 genes present in all selected FOVs

## Radial distribution of genes/cell types
This code can visualise a distribution as a function of angle and distance from the reference point. This can be either gene distribution or distribtuion of any other (categorical) entity from 'adata.obs'. For an example check **example_spatial_distribution.ipynb**
**plot_spatial_distribution** function has next inputs (*mandatory inputs):
 - adata*: AnnData object
 - ref_point*: either list or numpy array with 2 numbers - x,y position of the reference point 
 - sample_name: sample name (for the plot title)
 - save_folder: folder to save output images
 - angle_bin_size: bin size for the angular distribution
 - distance_bin_size: bin size for the radial distribution
 - distance_max: maximum distance for the radial distribution
 - gene_list: list of teh genes which distributions you want to plot
 - obs_sets: sets of columns and values to be plotted in the format: `obs_sets = [['column1', ['col1_val1', 'col1_val2']], ['column2', ['col2_val1', 'col2_val2']]]`
 - colormap: colormap for the heatmap histogram distribution
