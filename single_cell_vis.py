import numpy as np
import matplotlib.pyplot as plt
import spatialdata as sd
from spatialdata_io import xenium
from tifffile import TiffFile
import pandas as pd


def plot_one_img(ax, full_img, polygon, pixelsize, text, value, df_transcripts, top5_genes, border_size, label_col, label_width, spot_size, text_size, add_leg, display_genes):
    x, y, x_0, x_1, y_0, y_1 = define_box_boundaries_and_label_pts(polygon, pixelsize, border_size)
    #print(x_0); print(x_1); print(y_0); print(y_1)
    sub_img = full_img[0,y_0:y_1, x_0:x_1].compute().to_numpy()
    #print(sub_img)
    ax.imshow(sub_img, cmap = 'gray', vmin = np.min(sub_img), vmax = np.max(sub_img), origin = 'lower')
    ax.plot(x-x_0,y-y_0, color = label_col, linewidth = label_width)
    plot_transcripts_only(ax, df_transcripts, x_0, x_1, y_0, y_1, top5_genes, pixelsize, spot_size, text_size, add_leg)
    if type(value)==str:
        text_str = str(text) + '=' + value
    else:
        text_str = str(text) + '=' + "{0:,.2f}".format(value)
    ax.set_title(text_str)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim([0, x_1-x_0])
    ax.set_ylim([0, y_1-y_0])

def define_box_boundaries_and_label_pts(polygon, pixelsize, border_size):
    xx,yy = polygon.exterior.xy; 
    x = (np.array(xx))/pixelsize;y = (np.array(yy))/pixelsize
    x_0 = int(np.min(x)-border_size); x_1 = int(np.max(x)+border_size);
    y_0 = int(np.min(y)-border_size); y_1 = int(np.max(y)+border_size);
    return x, y, x_0, x_1, y_0, y_1

def find_pixelsize_from_image(folder_path):
    try:
        with TiffFile(folder_path + '/morphology_mip.ome.tif') as tif:
            meta = tif.ome_metadata
            index = meta.find('PhysicalSizeX')
            pixelsize = float(meta[index+15:index+21])
    except:
        with TiffFile(folder_path + '/morphology_focus.ome.tif') as tif:
            meta = tif.ome_metadata
            index = meta.find('PhysicalSizeX')
            pixelsize = float(meta[index+15:index+21])
    return pixelsize

def get_only_transcripts_around_cells(transcripts_all, polygons_df, pixelsize, border_size):
    df_transcripts = pd.DataFrame()
    for polygon in polygons_df['geometry']:
        x, y, x_0, x_1, y_0, y_1 = define_box_boundaries_and_label_pts(polygon, pixelsize, border_size)
        subtr = transcripts_all[(transcripts_all['x']/pixelsize>x_0) & (transcripts_all['x']/pixelsize<x_1) & (transcripts_all['y']/pixelsize>y_0) & (transcripts_all['y']/pixelsize<y_1)]
        df_transcripts = pd.concat([df_transcripts, subtr], ignore_index=True)
    return df_transcripts
        
def determine_top_5_genes(df_transcripts):
    counts = df_transcripts.groupby('feature_name').size()
    top5_genes = list(counts.sort_values(ascending = False)[:5].index.astype('str'))
    return top5_genes

def plot_transcripts_only(ax, df_transcripts, x_0, x_1, y_0, y_1, top5_genes, pixelsize, spot_size, size_text, add_legend = True):
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231']
    df_transcripts_sub = df_transcripts[(df_transcripts['x']/pixelsize>x_0) & (df_transcripts['x']/pixelsize<x_1) & (df_transcripts['y']/pixelsize>y_0) & (df_transcripts['y']/pixelsize<y_1)]
    #print(df_transcripts_sub)
    for gene,col in zip(top5_genes,colors):
        df_transcripts_sub_gene = df_transcripts_sub[df_transcripts_sub['feature_name']==gene]
        #print(df_transcripts_sub_gene)
        if df_transcripts_sub_gene.shape[0]>0:
            ax.scatter(df_transcripts_sub_gene['x']/pixelsize-x_0, df_transcripts_sub_gene['y']/pixelsize-y_0, c = col, s = spot_size)
            
    #add "legend"
    if add_legend:
        i=1;
        for gene,col in zip(top5_genes,colors):
            ax.text(0, (y_1-y_0)-i*(y_1-y_0)/20, gene, fontsize = size_text, c = col,  weight='bold')
            i+=1

def get_N_random_polygons(sdata_sub_obs, N, polygons_df, col_name):
    sdata_sub_obs_rand = sdata_sub_obs.sample(N)
    #pol_df = polygons_df(polygons_df.index == sdata_sub_obs_rand['cell_id'])
    pol_df = polygons_df.loc[sdata_sub_obs_rand['cell_id']]
    values_df = sdata_sub_obs_rand[col_name]
    return pol_df, values_df


def plot_cells_by_par(sdata, n_hor, n_vert, col_name, min_value, max_value, transcripts_all, pixelsize, border_size_px = 10, label_col = 'magenta', label_width = 2, save = 'None', spot_size = 5, text_size = 10, display_genes = []):
    print('Preparation and filtering')
    if type(min_value)==str:
        sdata_sub_obs = sdata.tables['table'].obs[sdata.tables['table'].obs[col_name] == min_value]
    else:
        sdata_sub_obs = sdata.tables['table'].obs[(sdata.tables['table'].obs[col_name]>=min_value) & (sdata.tables['table'].obs[col_name]<max_value)]
    
    polygons_df, values_df = get_N_random_polygons(sdata_sub_obs, n_hor*n_vert, sdata.shapes['cell_boundaries'], col_name)
    
    #transcripts
    df_transcripts = get_only_transcripts_around_cells(transcripts_all, polygons_df, pixelsize, border_size_px)
    #print(df_transcripts)
    if len(display_genes)==0:
        display_genes = determine_top_5_genes(df_transcripts)
        
    fig, axs = plt.subplots(n_vert, n_hor, figsize = (n_hor*10, n_vert*10))
    #print(polygons_df); 
    i=0; add_legend = True
    for val,polygon in zip(values_df, polygons_df['geometry']):
        #print(i//n_hor)
        #print(i%n_hor)
        print('Plotting subimage ' + str(i), end = '\r')
        #print(polygon)
        if n_hor!=1 and n_vert!=1:
            ax = axs[i//n_hor,i%n_hor]
        else:
            ax = axs[i]
        i+=1
        try:
            plot_one_img(ax, sdata.images['morphology_mip'].scale0.image, polygon, pixelsize, col_name, val, df_transcripts, display_genes, border_size_px, label_col, label_width, spot_size, text_size, add_legend, display_genes)
        except:
            plot_one_img(ax, sdata.images['morphology_focus'].scale0.image, polygon, pixelsize, col_name, val, df_transcripts, display_genes, border_size_px, label_col, label_width, spot_size, text_size, add_legend, display_genes)
        add_legend = False
    if save!='None':
        plt.savefig(save, dpi = 300)
