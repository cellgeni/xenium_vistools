import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines
import os
from IPython.display import display, clear_output
import pandas as pd


def pick_line_and_reference(adata, save=None, callback=None):
    """
    Jupyter notebook-friendly point picker with callback.
    
    Args:
        adata: AnnData object
        save: optional file path to save figure
        callback: function to call with the (p1, p2, p_proj) result
    """
    spatial_key = 'spatial' if 'spatial' in adata.obsm else 'X_spatial'
    spatial = adata.obsm[spatial_key]

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(spatial[:, 0], spatial[:, 1], s=1)
    ax.set_title('Double-click 2 points for line + 1 for reference')
    ax.axis('equal')
    coords = []

    def onclick(event):
        if event.dblclick:
            coords.append((event.xdata, event.ydata))
            ax.plot(event.xdata, event.ydata, 'ro')
            fig.canvas.draw()

            if len(coords) == 2:
                line = mlines.Line2D(
                    [coords[0][0], coords[1][0]],
                    [coords[0][1], coords[1][1]],
                    color='red'
                )
                ax.add_line(line)
                fig.canvas.draw()

            if len(coords) == 3:
                fig.canvas.mpl_disconnect(cid)
                clear_output(wait=True)
                if save:
                    plt.savefig(save, dpi=300)
                #plt.close(fig)
                result = process_coords(coords)
                if callback:
                    callback(result)

    def process_coords(coords):
        p1 = np.array(coords[0])
        p2 = np.array(coords[1])
        p_ref = np.array(coords[2])
        line_vec = p2 - p1
        line_vec_norm = line_vec / np.linalg.norm(line_vec)
        proj_length = np.dot(p_ref - p1, line_vec_norm)
        p_proj = p1 + proj_length * line_vec_norm
        print(f"P1: {p1}\nP2: {p2}\nReference (corrected): {p_proj}")
        return p1, p2, p_proj

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    #display(fig)

def plot_group_distribution_along_line(adata, picked_points, xboundaries = (None, None), selected_groups = None, bins=30, column_name='group', color = None, units = 'pixels', grid = True, save_folder = None):
    """
    Plot distribution of group categories along the picked line.

    Args:
        adata (AnnData): AnnData object with spatial and group info.
        picked_points (tuple): (p1, p2, reference_point) from pick_line_and_reference().
        bins (int): Number of bins for histogram.
    """
    # Unpack points
    p1, p2, ref_point = picked_points

    # Pick spatial coordinates
    spatial_key = 'spatial' if 'spatial' in adata.obsm else 'X_spatial'
    spatial = adata.obsm[spatial_key]
    
    # Compute line unit vector
    line_vec = p2 - p1
    line_vec_norm = line_vec / np.linalg.norm(line_vec)

    # Project all points onto the line
    projections = np.dot(spatial - ref_point, line_vec_norm)

    # Prepare dataframe for plotting
    df = pd.DataFrame({
        'distance_along_line': projections,
        column_name: adata.obs[column_name].values
    })

    # Plot
    
    
    # Plot histogram for each group separately
    if selected_groups is not None:
        #plt.figure(figsize=(5*len(selected_groups), 6))
        #f, axs = plt.subplots(len(selected_groups),1, figsize=(5*len(selected_groups), 9))
        for selected_group in selected_groups:
            f, ax = plt.subplots(figsize=(8, 8))
            subset = df[df[column_name] == selected_group]
            if xboundaries == (None, None):
                ax.hist(
                    subset['distance_along_line'],
                    bins=bins,
                    alpha=0.8,
                    label=str(selected_group),
                    histtype='stepfilled',
                    color = color,
                )
            else:
                ax.hist(
                    subset['distance_along_line'],
                    bins=bins,
                    alpha=0.8,
                    label=str(selected_group),
                    histtype='stepfilled',
                    color = color,
                    range = xboundaries,
                )
            ax.set_title("Distribution of " + str(selected_group) + " along the selected line")
            ax.axvline(0, color='black', linestyle='--', label='Reference Point')
            ax.set_xlabel('Distance along the line, ' + str(units))
            ax.set_ylabel('Number of cells')
            ax.legend()
            
            if grid:
                ax.grid(True)

            if save_folder:
                path = os.path.join(save_folder, f"{selected_group}_distribution.png")
                plt.savefig(path, dpi=300)
    else:
        plt.figure(figsize=(12, 12))
        groups = df[column_name].unique()
        for group in groups:
            subset = df[df[column_name] == group]
            plt.hist(
                subset['distance_along_line'],
                bins=bins,
                alpha=0.6,
                label=str(group),
                histtype='stepfilled'
            )
        
        #plt.title("Distribution of " + str(selected_group) + " along the selected line")
        plt.axvline(0, color='black', linestyle='--', label='Reference Point')
        plt.xlabel('Distance along the line, ' + str(units))
        plt.ylabel('Number of cells')
        plt.legend()
        if grid:
             plt.grid(True)

    
    
    plt.show()


def plot_gene_expression_along_line(adata, picked_points, xboundaries = (None, None), selected_genes = [None], bins=30, column_name='group', color = None, units = 'pixels', grid = True, save_folder = None):
    """
    Plot distribution of group categories along the picked line.

    Args:
        adata (AnnData): AnnData object with spatial and group info.
        picked_points (tuple): (p1, p2, reference_point) from pick_line_and_reference().
        bins (int): Number of bins for histogram.
    """
    # Unpack points
    p1, p2, ref_point = picked_points

    # Pick spatial coordinates
    spatial_key = 'spatial' if 'spatial' in adata.obsm else 'X_spatial'
    spatial = adata.obsm[spatial_key]
    
    # Compute line unit vector
    line_vec = p2 - p1
    line_vec_norm = line_vec / np.linalg.norm(line_vec)

    # Project all points onto the line
    projections = np.dot(spatial - ref_point, line_vec_norm)

    # Prepare dataframe for plotting
    df = pd.DataFrame({
        'distance_along_line': projections,
        column_name: adata.obs[column_name].values
    })

    # Plot
    
    
    # Plot histogram for each group separately
    if selected_genes is not None:
        #plt.figure(figsize=(5*len(selected_groups), 6))
        #f, axs = plt.subplots(len(selected_groups),1, figsize=(5*len(selected_groups), 9))
        for gene in selected_genes:
            f, ax = plt.subplots(figsize=(8, 8))
            gene_index = adata.var_names.get_loc(gene)
            gene_expression = adata.X[:, gene_index].flatten() if issubclass(type(adata.X), np.ndarray) else adata.X[:, gene_index].toarray().flatten()
            if xboundaries==(None, None):
                ax.hist(
                    df['distance_along_line'],
                    bins=bins,
                    alpha=0.8,
                    label=str(gene),
                    histtype='stepfilled',
                    color = color,
                    weights=gene_expression
                )
            else:
                ax.hist(
                    df['distance_along_line'],
                    bins=bins,
                    alpha=0.8,
                    label=str(gene),
                    histtype='stepfilled',
                    color = color,
                    range = xboundaries,
                    weights=gene_expression
                )
            ax.set_title("Distribution of " + str(gene) + " expression along the selected line")
            ax.axvline(0, color='black', linestyle='--', label='Reference Point')
            ax.set_xlabel('Distance along the line, ' + str(units))
            ax.set_ylabel('Gene expression')
            ax.legend()
            
            if grid:
                ax.grid(True)

            if save_folder:
                path = os.path.join(save_folder, f"{gene}_distribution.png")
                plt.savefig(path, dpi=300)
    else:
        print('No genes has been chosen to display!')
        
        
        plt.axvline(0, color='black', linestyle='--', label='Reference Point')
        plt.xlabel('Distance along the line, ' + str(units))
        plt.ylabel('Number of cells')
        plt.legend()
        if grid:
             plt.grid(True)

    if xboundaries:
        plt.xlim(xboundaries)
    
    plt.show()

