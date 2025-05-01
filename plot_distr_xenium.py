import numpy as np
import pandas as pd
import anndata
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

def polar_density_plot_hist(hist_list, labels, angle_bins, distance_bins,
                                 colormap='Set1', save = None, distance_max = None):
    """
    Plot multiple polar histograms of density along angle and distance.

    Args:
        data_angles_list: list of np.ndarray, angles (radians) for each histogram
        data_distances_list: list of np.ndarray, distances for each histogram
        labels: list of strings, labels for each histogram (legend)
        angle_bin_size: size of angle bins (radians)
        distance_bin_size: size of distance bins (same units as input)
        cmap: matplotlib colormap name for distinct colors
    """
    if not isinstance(hist_list, list):
        hist_list = [hist_list]
        labels = [labels]

    angle_bin_centers = (angle_bins[:-1] + angle_bins[1:]) / 2
    angle_bin_size = angle_bins[1] - angle_bins[0]

    # Define bins
    angle_bins = np.arange(0, 2*np.pi + angle_bin_size, angle_bin_size)
    n_sets = len(hist_list)
    

    # Prepare plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 8))
    if n_sets == 1:
        hist = hist_list[0]
        angle_bin_centers = (angle_bins[:-1] + angle_bins[1:]) / 2
        distance_bin_centers = (distance_bins[:-1] + distance_bins[1:]) / 2
        cmap = plt.get_cmap(colormap)
        norm = mpl.colors.Normalize(vmin=hist.min(), vmax=hist.max())
        for i in range(len(distance_bin_centers)):
            radii = hist[:, i]
            colors = cmap(norm(radii))  # Map the radii to colors
            bars = ax.bar(angle_bin_centers, radii, width=(angle_bins[1] - angle_bins[0]), bottom=i*(distance_bins[1] - distance_bins[0]), color=colors, alpha=0.75)
        ax.set_xticks(np.arange(0, 2.0 * np.pi, angle_bin_size))

        ax.set_title(labels[0])
        # Create a colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.05, pad=0.05)
        cbar.set_label('Density')
    else:
        print('multiple histograms')
        colors = plt.get_cmap(colormap).colors  
        for i in range(n_sets):
            
            label = labels[i] if labels else f"Hist {i+1}"
            hist = hist_list[i]
            distance_bin_size = distance_bins[1] - distance_bins[0]
            # Plot each distance bin as a ring of bars
            for j in range(hist.shape[1]):
                radii = hist[:, j]
                theta_centers = angle_bins[:-1] + angle_bin_size / 2
                bar_width = angle_bin_size
                bottom = j * distance_bin_size
                ax.bar(theta_centers, radii, width=bar_width, bottom=bottom,
                    color=colors[i % len(colors)], alpha=1/n_sets, label=label if j == 0 else None)

        #ax.set_theta_zero_location('E')
        #ax.set_theta_direction(1)
        ax.set_title("Multiple gene expression distribution")
        ax.legend(loc='upper right')
        plt.show()
    
    if save:
        plt.savefig(save, dpi = 300)


def polar_density_plot(data_angles_list, data_distances_list, labels=None,
                                 angle_bin_size=np.pi/12, distance_bin_size=1000,
                                 colormap='Set1', save = None, distance_max = None):
    """
    Plot multiple polar histograms of density along angle and distance.

    Args:
        data_angles_list: list of np.ndarray, angles (radians) for each histogram
        data_distances_list: list of np.ndarray, distances for each histogram
        labels: list of strings, labels for each histogram (legend)
        angle_bin_size: size of angle bins (radians)
        distance_bin_size: size of distance bins (same units as input)
        cmap: matplotlib colormap name for distinct colors
    """
    if not isinstance(data_angles_list, list):
        data_angles_list = [data_angles_list]
        data_distances_list = [data_distances_list]
        labels = [labels]

    assert len(data_angles_list) == len(data_distances_list), "Angle and distance lists must match"
    n_sets = len(data_angles_list)

    # Define bins
    angle_bins = np.arange(0, 2*np.pi + angle_bin_size, angle_bin_size)
    if distance_max:
        distance_max = distance_max
    else:
        distance_max = max([np.max(d) for d in data_distances_list])
    
    distance_bins = np.arange(0, distance_max + distance_bin_size, distance_bin_size)

    # Prepare plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 8))
    if n_sets == 1:
        hist, _, _ = np.histogram2d(data_angles_list[0], data_distances_list[0], bins=(angle_bins, distance_bins))
        distance_bins = np.arange(0, distance_max + distance_bin_size, distance_bin_size)
        angle_bin_centers = (angle_bins[:-1] + angle_bins[1:]) / 2
        distance_bin_centers = (distance_bins[:-1] + distance_bins[1:]) / 2
        cmap = plt.get_cmap(colormap)
        norm = mpl.colors.Normalize(vmin=hist.min(), vmax=hist.max())
        for i in range(len(distance_bin_centers)):
            radii = hist[:, i]
            colors = cmap(norm(radii))  # Map the radii to colors
            bars = ax.bar(angle_bin_centers, radii, width=(angle_bins[1] - angle_bins[0]), bottom=i*(distance_bins[1] - distance_bins[0]), color=colors, alpha=0.75)
        ax.set_xticks(np.arange(0, 2.0 * np.pi, angle_bin_size))

        ax.set_title(labels[0])
        # Create a colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.05, pad=0.05)
        cbar.set_label('Density')
    else:
        colors = plt.get_cmap(colormap).colors  
        for i in range(n_sets):
            angles = data_angles_list[i]
            distances = data_distances_list[i]
            label = labels[i] if labels else f"Hist {i+1}"

            hist, _, _ = np.histogram2d(angles, distances, bins=[angle_bins, distance_bins])

            # Plot each distance bin as a ring of bars
            for j in range(hist.shape[1]):
                radii = hist[:, j]
                theta_centers = angle_bins[:-1] + angle_bin_size / 2
                bar_width = angle_bin_size
                bottom = j * distance_bin_size
                ax.bar(theta_centers, radii, width=bar_width, bottom=bottom,
                    color=colors[i % len(colors)], alpha=1/n_sets, label=label if j == 0 else None)

        #ax.set_theta_zero_location('E')
        #ax.set_theta_direction(1)
        ax.set_title("Polar Density Plot (Multiple Histograms)")
        ax.legend(loc='upper right')
        
    plt.show()
    if save:
        plt.savefig(save, dpi = 300)

def polar_density_plot_hist_old(hist, angle_bins, distance_bins, title, save = None, colormap = 'Purples'):
    """
    Create a polar histogram plot with color scale using precomputed histogram data.
    
    Args:
    hist (np.ndarray): A 2D array with the histogram data.
    angle_bins (np.array): An array of bins for the angles (radians).
    distance_bins (np.array): An array of bins for the distances.
    title (str): Title of the plot.
    colormap (str): Name of the colormap to use.
    """
    # Ensure angle bins cover full circle if they do not
    
    angle_bin_centers = (angle_bins[:-1] + angle_bins[1:]) / 2
    angle_bin_size = angle_bins[1] - angle_bins[0]
    # Plotting
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    norm = mpl.colors.Normalize(vmin=hist.min(), vmax=hist.max())
    cmap = plt.get_cmap(colormap)
    
    for i in range(len(distance_bins) - 1):
        radii = hist[:, i]
        colors = cmap(norm(radii))  # Map the radii to colors
        bars = ax.bar(angle_bin_centers, radii, width=angle_bin_size, bottom=i*(distance_bins[1] - distance_bins[0]), color=colors, alpha=0.75)
    ax.set_xticks(np.arange(0, 2.0 * np.pi, angle_bin_size))
    '''
    for i in range(3,5):
        radii = hist[:, i]
        colors = cmap(norm(radii))  # Map the radii to colors
        print(angle_bin_centers)
        #print(radii)
        bars = ax.bar(angle_bin_centers, radii, width=(angle_bins[1] - angle_bins[0]), bottom=i*(distance_bins[1] - distance_bins[0]), color=colors, alpha=0.75)
    '''
    # Create a colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.05, pad=0.05)
    cbar.set_label('Density')

    ax.set_title(title)
    plt.show()
    if save:
        plt.savefig(save, dpi = 300)

def calculate_distances_and_angles(ref_point, points):
    """
    Calculate the distances and angles from a reference point to each point in a given matrix.
    
    Args:
    ref_point (tuple): A tuple (x, y) representing the reference point coordinates.
    points (np.array): A 2D numpy array where each row is (x, y) coordinates of a point.
    
    Returns:
    tuple: Two numpy arrays, first with angles in radians, second with distances.
    """
    # Extract coordinates from the input
    ref_x, ref_y = ref_point
    points_x = points[:, 0]
    points_y = points[:, 1]
    
    # Calculate differences
    dx = points_x - ref_x
    dy = points_y - ref_y
    
    # Calculate distances using the Euclidean distance formula
    distances = np.sqrt(dx**2 + dy**2)
    
    # Calculate angles using arctan2 to handle different quadrants correctly
    angles = np.arctan2(dy, dx)
    angles+= np.pi
    return angles, distances

def gene_expression_distribution(adata, gene_name, ref_point, angle_bin_size, distance_bin_size, distance_max = None):
    """
    Calculate the distribution of a specified gene's expression as a function of angle and distance from a reference point.

    Args:
    adata (anndata.AnnData): The annotated data matrix.
    gene_name (str): The name of the gene to analyze.
    ref_point (tuple): A tuple (x, y) representing the reference point coordinates.
    angle_bin_size (float): The size of each angle bin in radians.
    distance_bin_size (float): The size of each distance bin.

    Returns:
    np.ndarray: A 2D array with the distribution of gene expression across angle and distance bins.
    """
    # Extract spatial coordinates and gene expression
    coords = adata.obsm['spatial']  # adjust column names based on your data
    gene_index = adata.var_names.get_loc(gene_name)
    gene_expression = adata.X[:, gene_index].flatten() if issubclass(type(adata.X), np.ndarray) else adata.X[:, gene_index].toarray().flatten()

    # Calculate distances and angles
    dx = coords[:, 0] - ref_point[0]
    dy = coords[:, 1] - ref_point[1]
    distances = np.sqrt(dx**2 + dy**2)
    angles = np.arctan2(dy, dx)
    angles+= np.pi
    # Define bins
    if distance_max:
        max_distance = distance_max
    else:
        max_distance = np.max(distances)
    
    distance_bins = np.arange(0, max_distance + distance_bin_size, distance_bin_size)
    angle_bins = np.arange(0, 2*np.pi+angle_bin_size, angle_bin_size)

    # Histogram the data into bins
    hist, _, _ = np.histogram2d(angles, distances, bins=(angle_bins, distance_bins), weights=gene_expression)

    return hist, angle_bins, distance_bins

def get_positions_from_obs_column(adata, column_name, value):
    """
    Filter an AnnData object to get spatial coordinates for observations that match a specific value in a column.

    Args:
    adata (anndata.AnnData): The annotated data matrix.
    column_name (str): The column in adata.obs to filter by.
    value: The value to match in the column.

    Returns:
    np.ndarray: Array of spatial coordinates for matching observations.
    """
    if column_name not in adata.obs.columns:
        raise ValueError(f"Column '{column_name}' does not exist in adata.obs.")
    
    if 'spatial' not in adata.obsm:
        raise ValueError("Spatial data is not available in adata.obsm['spatial'].")

    # Filter the observations where the column matches the specified value
    matching_indices = adata.obs[column_name] == value
    
    # Retrieve the spatial data for these indices
    spatial_data = adata.obsm['spatial'][matching_indices]

    return spatial_data

def plot_spatial_distribution(adata, ref_point, sample_name = None, save_folder = None, angle_bin_size = np.pi/4, distance_bin_size = 5000, 
                              distance_max = None, gene_list = [], obs_sets = [], colormap = 'Purples', multiple_histograms = False,
                              group_together = False):
    ## structure of obs sets: [['column1', ['col1_val1', 'col1_val2']], ['column2', ['col2_val1', 'col2_val2']]]

    #make sure we have adata.obsm['spatial']
    if 'spatial' not in adata.obsm:
        adata.obsm['spatial'] = adata.obsm['X_spatial']
    
    ## lets plot genes if there are any in the list
    if gene_list:
        if not multiple_histograms:
        #one plot per gene
            for gene in gene_list:
                print(gene)
                hist, angle_bins, distance_bins = gene_expression_distribution(adata, gene, ref_point,  angle_bin_size, distance_bin_size, distance_max)
                if sample_name: 
                    plot_title = sample_name + ' -- ' + str(gene)
                else:
                    plot_title = str(gene)
                if save_folder:
                    save_path = os.path.join(save_folder, plot_title + '.png')
                    polar_density_plot_hist(hist, angle_bins=angle_bins, distance_bins=distance_bins, labels=plot_title, save = save_path, colormap = colormap)
                else:
                    polar_density_plot_hist(hist, angle_bins=angle_bins, distance_bins=distance_bins, labels=plot_title, colormap = colormap)
        else:
            #plot all genes histograms in one plot
            hist_list = [];
            for gene in gene_list:
                    print(gene)
                    hist, angle_bins, distance_bins = gene_expression_distribution(adata, gene, ref_point,  angle_bin_size, distance_bin_size, distance_max)
                    hist_list.append(hist)

            if not group_together:
            #we plot them all with different colors
                if save_folder:
                    save_path = os.path.join(save_folder, 'genes.png')
                    polar_density_plot_hist(hist_list, angle_bins=angle_bins, distance_bins=distance_bins, labels=gene_list, save = save_path)
                else:
                    polar_density_plot_hist(hist_list, angle_bins=angle_bins, distance_bins=distance_bins, labels=gene_list)
            else:
            #we group them together
                hist_sum = np.mean(hist_list, axis=0)
                plot_title = "grouped genes " + "_".join(gene_list)
                if save_folder:
                    save_path = os.path.join(save_folder, plot_title + '.png')
                    polar_density_plot_hist(hist_sum, angle_bins=angle_bins, distance_bins=distance_bins, labels=plot_title, save = save_path, colormap = colormap)
                else:
                    polar_density_plot_hist(hist_sum, angle_bins=angle_bins, distance_bins=distance_bins, labels=plot_title, colormap = colormap)

            
            

    ## lets plot any adata.obs columns
    if obs_sets:
        for columns in obs_sets:
            column_name = columns[0]
            if not multiple_histograms:
                #one plot per one value
                for value in columns[1]:
                    cell_positions = get_positions_from_obs_column(adata, column_name, value)
                    angles, distances = calculate_distances_and_angles(ref_point, cell_positions)
                    if sample_name: 
                        plot_title = sample_name + ' -- ' + str(value)
                    else:
                        plot_title = str(value)
                    
                    if save_folder:
                        save_path = os.path.join(save_folder, plot_title + '.png')
                        polar_density_plot(angles, distances, angle_bin_size=angle_bin_size, distance_bin_size=distance_bin_size, labels=plot_title, save = save_path, colormap = colormap, distance_max = distance_max)
                    else:
                        polar_density_plot(angles, distances, angle_bin_size=angle_bin_size, distance_bin_size=distance_bin_size, labels=plot_title, colormap = colormap, distance_max = distance_max)
            else:
                #plot all values histograms in one plot
                angles_list = []; distances_list = []; values_list = []
                for value in columns[1]:
                    cell_positions = get_positions_from_obs_column(adata, column_name, value)
                    angles, distances = calculate_distances_and_angles(ref_point, cell_positions)
                    angles_list.append(angles)
                    distances_list.append(distances)
                    values_list.append(value)
                if not group_together:
                #we plot them all with different colors
                    if save_folder:
                        save_path = os.path.join(save_folder, column_name + '.png')
                        polar_density_plot(angles_list, distances_list, angle_bin_size=angle_bin_size, distance_bin_size=distance_bin_size, labels=values_list, save = save_path, distance_max = distance_max)
                    else:S
                        polar_density_plot(angles_list, distances_list, angle_bin_size=angle_bin_size, distance_bin_size=distance_bin_size, labels=values_list, distance_max = distance_max)
                else:
                #we group them together
                    angles = np.concatenate(angles_list)
                    distances = np.concatenate(distances_list)
                    plot_title = "grouped " + column_name + "-" + "_".join(values_list)
                    if save_folder:
                        save_path = os.path.join(save_folder, plot_title + '.png')
                        polar_density_plot(angles, distances, angle_bin_size=angle_bin_size, distance_bin_size=distance_bin_size, labels=plot_title, save = save_path, distance_max = distance_max, colormap = colormap)
                    else:
                        polar_density_plot(angles, distances, angle_bin_size=angle_bin_size, distance_bin_size=distance_bin_size, labels=plot_title, distance_max = distance_max, colormap = colormap)
class InteractivePlot:
    def __init__(self, points, size = 1):
        self.fig, self.ax = plt.subplots()
        self.ax.scatter(points[:, 0], points[:, 1], s = size)
        self.clicked_position = None
        plt.axis('equal')
        # Text annotation for displaying click position
        self.text = self.ax.text(0.5, 0.5, '', transform=self.ax.transAxes)
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)

    def onclick(self, event):
        # Update clicked position
        self.clicked_position = (event.xdata, event.ydata)
        self.text.set_text(f'Clicked at: {self.clicked_position[0]:.2f}, {self.clicked_position[1]:.2f}')
        self.fig.canvas.draw()

    def get_click_position(self):
        return self.clicked_position
    
