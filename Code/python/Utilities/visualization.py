import matplotlib
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from seaborn.palettes import _ColorPalette
import numpy as np


def adjust_spines(ax, spines, offset=3, smart_bounds=False):

    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', offset))
            spine.set_smart_bounds(smart_bounds)
        else:
            spine.set_color('None')  # don't draw spine

    # Turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # No y-axis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # No x-axis ticks
        ax.xaxis.set_ticks([])


# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


def center_diverging_colormap(vmin, vmax, cmap_name='seismic'):
    # Center colormap around 0
    # Copied from seaborn.matrix._determine_cmap_params() v0.9.0
    cmap = matplotlib.cm.get_cmap(cmap_name)
    center = 0
    values_range = max(vmax - center, center - vmin)
    normalizer = matplotlib.colors.Normalize(center - values_range, center + values_range)
    cmin, cmax = normalizer([vmin, vmax])
    cc = np.linspace(cmin, cmax, 256)
    cmap = matplotlib.colors.ListedColormap(cmap(cc))

    return cmap


def colorbar(mappable):
    # Get axes and figure of the selected image
    ax = mappable.axes
    fig = ax.figure
    # Make an axis on the right of the image
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.2)
    # Add colorbar
    cbar = fig.colorbar(mappable, cax=cax)
    # Set properties of ticks
    cbar.ax.tick_params(axis='y', which='major', direction='in', width=.5, length=10, labelsize=8)
    return cbar


def make_seaborn_palette(list_colors):
    list_colors = np.array(list_colors, dtype=float) / 255
    palette = map(colors.colorConverter.to_rgb, list_colors)
    palette = _ColorPalette(palette)
    return palette
