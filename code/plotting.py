import numpy as np
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from helper import find_geo_coords

def global_hatched_map(ax, data, sign_agreement, lsm, time, levels, cmap, title ):
    norm = BoundaryNorm(levels, ncolors=256, clip=False, extend='both')  # Create a BoundaryNorm

    lon_name, lat_name = find_geo_coords(data)
    # filled contour map
    im = ax.contourf(data[lon_name], data[lat_name], data, transform=ccrs.PlateCarree(), levels=levels, norm=norm, cmap=cmap, extend='both')

    # add hatching where models agree on sign
    if ( sign_agreement is not None):
        ax.contourf(data[lon_name], data[lat_name], sign_agreement.astype(int), 3, hatches=['', '////'], colors='none', transform=ccrs.PlateCarree())
        plt.rcParams['hatch.color'] = (0.0, 0.0, 0.0, 0.3)

    # add coastlines
    if (time == 'modern'):
        ax.contour(lsm.lon, lsm.lat, lsm, transform=ccrs.PlateCarree(), levels=[0.5], colors=['black'], linewidths=1)
        # Using features to add modern coastlines
        ax.add_feature(cfeature.COASTLINE, linestyle='--', edgecolor='black', linewidth=1)
    elif (time == 'eocene'):
        ax.contour(lsm.lon, lsm.lat, lsm, transform=ccrs.PlateCarree(), levels=[0.5], colors=['black'], linewidths=1)
    elif (time == 'modern_only'):
        ax.add_feature(cfeature.COASTLINE, linestyle='--', edgecolor='black', linewidth=1)

        
    # Add gridlines and labels
    gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.0, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}
    # ax.set_title(f'{exp_labels[col]} ({labels[row]})')
    ax.set_title(title, weight='bold', fontsize=12)

    # Generate the plot to enable access to the labels' attributes
    plt.draw()

    # Fix the right latitude labels
    # https://stackoverflow.com/questions/75597673/hide-right-side-axis-latitude-labels-for-high-latitude-non-rectangular-project/75600501#75600501
    # Iterate for the y-labels
    # The right labels have x coordinates > 0
    # The left labels < 0
    for ea in gl.ylabel_artists:
        right_label = ea.get_position()[0] > 0
        # print(ea, ea.get_position()[0], ea.get_visible())
        if right_label:
            ea.set_visible(False)

    return im

def common_colorbar(fig, im, levels, label, x=0.15, width=0.7):
    cbar_ax = fig.add_axes([x, 0.05, width, 0.03])  # Adjust the axes dimensions as necessary
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal', extend='both')
    cbar.set_label(label)

    # Set custom ticks and labels on the colorbar
    cbar.set_ticks(levels)
    cbar.set_ticklabels([str(level) for level in levels])

def custom_colormap(original_cmap='RdBu_r'):
    # Create a custom colormap with white around zero
    original_cmap = plt.colormaps[original_cmap]  # Get the original RdBu_r colormap
    newcolors = original_cmap(np.linspace(0, 1, 256))  # Get the color array
    midpoint = 0.5  # Find the relative position of zero
    whitening_intensity = 20  # Define range of colors to make white around zero
    start = int(midpoint * 256 - whitening_intensity)
    end = int(midpoint * 256 + whitening_intensity)
    newcolors[start:end, :] = [1, 1, 1, 1]  # Set the colors to white in the middle of the colormap
    custom_cmap = LinearSegmentedColormap.from_list('custom_RdBu', newcolors)

    return custom_cmap