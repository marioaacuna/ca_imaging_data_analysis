from __future__ import print_function, division, absolute_import

import itertools

import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.tight_layout import get_renderer


def _process_data(data, sort_by, sort_sets_by):
    # check all indices are vertical
    assert all(set([True, False]) >= set(level) for level in data.index.levels)
    if not data.index.is_unique:
        data = (data
                .groupby(level=list(range(data.index.nlevels)))
                .sum())

    data_tall = data.reset_index()
    data_tall_columns = [i for i in data_tall.columns if i != 'value']
    totals = [data_tall.loc[np.where(data_tall[i])[0], 'value'].sum() for i in data_tall_columns]
    totals = pd.Series(totals, index=data.index.names)
    if sort_sets_by == 'cardinality':
        totals.sort_values(ascending=False, inplace=True)
    elif sort_sets_by is not None:
        raise ValueError('Unknown sort_sets_by: %r' % sort_sets_by)
    data = data.reorder_levels(totals.index.values)

    if sort_by == 'cardinality':
        data = data.sort_values(ascending=False)
    elif sort_by == 'degree':
        comb = itertools.combinations
        o = pd.DataFrame([{name: True for name in names}
                          for i in range(data.index.nlevels + 1)
                          for names in comb(data.index.names, i)],
                         columns=data.index.names)
        o.fillna(False, inplace=True)
        o = o.astype(bool)
        o.set_index(data.index.names, inplace=True)
        # FIXME: should use reindex(index=...) ??
        data = data.loc[o.index]
    else:
        raise ValueError('Unknown sort_by: %r' % sort_by)

    min_value = 0
    max_value = np.inf
    data = data[np.logical_and(data >= min_value, data <= max_value)]

    return data, totals


class _Transposed:
    """Wrap an object in order to transpose some plotting operations

    Attributes of obj will be mapped.
    Keyword arguments when calling obj will be mapped.

    The mapping is not recursive: callable attributes need to be _Transposed
    again.
    """

    def __init__(self, obj):
        self.__obj = obj

    def __getattr__(self, key):
        return getattr(self.__obj, self._NAME_TRANSPOSE.get(key, key))

    def __call__(self, *args, **kwargs):
        return self.__obj(*args, **{self._NAME_TRANSPOSE.get(k, k): v
                                    for k, v in kwargs.items()})

    _NAME_TRANSPOSE = {
        'width': 'height',
        'height': 'width',
        'hspace': 'wspace',
        'wspace': 'hspace',
        'hlines': 'vlines',
        'vlines': 'hlines',
        'bar': 'barh',
        'barh': 'bar',
        'xaxis': 'yaxis',
        'yaxis': 'xaxis',
        'left': 'bottom',
        'right': 'top',
        'top': 'right',
        'bottom': 'left',
        'sharex': 'sharey',
        'sharey': 'sharex',
        'get_figwidth': 'get_figheight',
        'get_figheight': 'get_figwidth',
        'set_figwidth': 'set_figheight',
        'set_figheight': 'set_figwidth',
    }


def _transpose(obj):
    if isinstance(obj, str):
        return _Transposed._NAME_TRANSPOSE.get(obj, obj)
    return _Transposed(obj)


def _identity(obj):
    return obj


class UpSet:
    """Manage the data and drawing for a basic UpSet plot

    Primary public method is :meth:`plot`.

    Parameters
    ----------
    data : pandas.Series
        Values for each set to plot.
        Should have multi-index where each level is binary,
        corresponding to set membership.
    orientation : {'horizontal' (default), 'vertical'}
        If horizontal, intersections are listed from left to right.
    sort_by : {'cardinality', 'degree'}
        If 'cardinality', set intersections are listed from largest to
        smallest value.
        If 'degree', they are listed in order of the number of sets
        intersected.
    sort_sets_by : {'cardinality', None}
        Whether to sort the overall sets by total cardinality, or leave them
        in the provided order.
    facecolor : {str, list}
        Color for bar charts and dots.
    with_lines : bool
        Whether to show lines joining dots in the matrix, to mark multiple sets
        being intersected.
    element_size : float or None
        Side length in pt. If None, size is estimated to fit figure
    intersection_plot_elements : int
        The intersections plot should be large enough to fit this many matrix
        elements.
    totals_plot_elements : int
        The totals plot should be large enough to fit this many matrix
        elements.
    """

    def __init__(self, data, orientation='horizontal', sort_by='degree', sort_sets_by=None,
                 barcolor='darkgrey', facecolor='black', linecolor='darkgrey',
                 with_lines=True, element_size=32, dimmed_element_size=20,
                 intersection_labels=True, total_labels=True, label_format='%.1f',
                 intersection_plot_elements=10, totals_plot_elements=5):

        self._horizontal = orientation == 'horizontal'
        self._reorient = _identity if self._horizontal else _transpose
        self._barcolor = barcolor
        self._facecolor = facecolor
        self._linecolor = linecolor
        self._with_lines = with_lines
        self._element_size = element_size
        self._dimmed_element_size = dimmed_element_size
        self._intersection_labels = intersection_labels
        self._total_labels = total_labels
        self._label_format = label_format
        self._totals_plot_elements = totals_plot_elements
        self._intersection_plot_elements = intersection_plot_elements

        self.intersections, self.totals = _process_data(data, sort_by=sort_by, sort_sets_by=sort_sets_by)
        if not self._horizontal:
            self.intersections = self.intersections[::-1]

        self._n_sets = self.intersections.index.nlevels

    def _swapaxes(self, x, y):
        if self._horizontal:
            return x, y
        return y, x

    def make_grid(self, fig=None):
        """Get a SubplotSpec for each Axes, accounting for label text width
        """
        n_cats = len(self.totals)
        n_inters = len(self.intersections)

        if fig is None:
            fig = plt.gcf()

        # Determine text size to determine figure size / spacing
        r = get_renderer(fig)
        t = fig.text(0, 0, '\n'.join(self.totals.index.values), fontdict={'fontsize': 14})
        textw = t.get_window_extent(renderer=r).width
        t.remove()

        MAGIC_MARGIN = 10  # FIXME
        figw = self._reorient(fig.get_window_extent(renderer=r)).width
        if self._element_size is None:
            colw = (figw - textw - MAGIC_MARGIN) / (len(self.intersections) + self._totals_plot_elements)
        else:
            fig = self._reorient(fig)
            render_ratio = figw / fig.get_figwidth()
            colw = self._element_size / 72 * render_ratio
            figw = (colw * (len(self.intersections) + self._totals_plot_elements) + MAGIC_MARGIN + textw)
            fig.set_figwidth(figw / render_ratio)
            fig.set_figheight((colw * (n_cats + self._intersection_plot_elements)) / render_ratio)

        text_nelems = int(np.ceil(figw / colw - (len(self.intersections) + self._totals_plot_elements)))

        GS = self._reorient(matplotlib.gridspec.GridSpec)
        gridspec = GS(*self._swapaxes(n_cats + self._intersection_plot_elements, n_inters + text_nelems + self._totals_plot_elements), hspace=1)
        if self._horizontal:
            return {'intersections': gridspec[:-n_cats, -n_inters:],
                    'matrix': gridspec[-n_cats:, -n_inters:],
                    'shading': gridspec[-n_cats:, :],
                    'totals': gridspec[-n_cats:, :self._totals_plot_elements],
                    'gs': gridspec}
        else:
            return {'intersections': gridspec[-n_inters:, n_cats:],
                    'matrix': gridspec[-n_inters:, :n_cats],
                    'shading': gridspec[:, :n_cats],
                    'totals': gridspec[:self._totals_plot_elements, :n_cats],
                    'gs': gridspec}

    def plot_matrix(self, ax):
        """Plot the matrix of intersection indicators onto ax
        """
        ax = self._reorient(ax)
        data = self.intersections

        if not isinstance(self._facecolor, list):
            facecolor = [self._facecolor] * self._n_sets
        else:
            facecolor = list(self._facecolor)
        facecolor = np.array(facecolor, dtype=object)

        values = data.index.to_frame()[data.index.names].values
        col_idx = np.where(values)[1]
        idx = np.flatnonzero(values)
        non_idx = np.setdiff1d(np.arange(values.ravel().shape[0]), idx)
        c = np.array(['lightgrey'] * data.shape[0] * self._n_sets, dtype=object)
        c[idx] = facecolor[col_idx]
        x = np.repeat(np.arange(len(data)), self._n_sets)
        y = np.tile(np.arange(self._n_sets), len(data))
        if self._element_size is not None:
            s = (self._element_size * .35) ** 2
            s_dimmed = (self._dimmed_element_size * .35) ** 2
        else:
            # TODO: make s relative to colw
            s = 200
        x, y = self._swapaxes(x, y)
        # Plot non-highlighted points
        ax.scatter(x[non_idx], y[non_idx], c=c[non_idx], linewidth=0, s=s_dimmed, zorder=1)
        ax.scatter(x[idx], y[idx], c=c[idx], linewidth=0, s=s, zorder=3)

        if self._with_lines:
            line_data = (pd.Series(y[idx], index=x[idx])
                         .groupby(level=0)
                         .aggregate(['min', 'max']))
            ax.vlines(line_data.index.values,
                      line_data['min'], line_data['max'],
                      lw=2, colors=self._linecolor, zorder=2)

        tick_axis = ax.yaxis
        tick_axis.set_ticks(np.arange(self._n_sets))
        ax.set_yticklabels(data.index.names, rotation=0 if self._horizontal else -90, fontdict={'fontsize': 14})
        ax.xaxis.set_visible(False)
        ax.tick_params(axis='both', which='both', length=0)
        if not self._horizontal:
            ax.yaxis.set_ticks_position('top')
        ax.set_frame_on(False)

        # Invert axis so first group is on top
        ax.invert_yaxis()

    def plot_intersections(self, ax):
        """Plot bars indicating intersection size
        """
        ax = self._reorient(ax)
        # Plot bars corresponding to intersections
        values = self.intersections.index.to_frame()[self.intersections.index.names].values
        n_values = np.sum(values.astype(int), axis=1)
        idx_not1 = np.where(n_values != 1)[0]
        ax.bar(idx_not1, self.intersections.iloc[idx_not1], .5, color=self._barcolor, zorder=10, align='center')
        # Plot bars corresponding to unique sets
        idx1 = np.where(n_values == 1)[0]
        color1 = np.where(values[idx1, :])[1]
        if not isinstance(self._facecolor, list):
            facecolor = [self._facecolor] * self._n_sets
        else:
            facecolor = list(self._facecolor)
        facecolor = np.array(facecolor, dtype=object)
        ax.bar(idx1, self.intersections.iloc[idx1], .5, color=facecolor[color1], zorder=10, align='center')
        # ax.bar(np.arange(len(self.intersections)), self.intersections, .5, color=self._barcolor, zorder=10, align='center')
        ax.xaxis.set_visible(False)
        for x in ['top', 'bottom', 'right']:
            ax.spines[self._reorient(x)].set_visible(False)

        tick_axis = ax.yaxis
        tick_axis.grid(True)
        tick_axis.set_tick_params(length=0)

        # Add labels
        if self._intersection_labels:
            ylims = ax.get_ylim()
            yrange = abs(ylims[1] - ylims[0])
            ypadding = yrange * .05
            # Convert numbers to strings
            values = self.intersections.values
            values_str = [self._label_format % i for i in values]
            for ii in range(len(values)):
                ax.text(ii, values[ii] + ypadding, values_str[ii], horizontalalignment='center', verticalalignment='top')

    def plot_totals(self, ax):
        """Plot bars indicating total set size
        """
        orig_ax = ax
        ax = self._reorient(ax)
        ax.barh(np.arange(len(self.totals.index.values)), self.totals, .5, color=self._facecolor, align='center')
        max_total = self.totals.max()
        if self._horizontal:
            orig_ax.set_xlim(max_total, 0)
        for x in ['top', 'left', 'right']:
            ax.spines[self._reorient(x)].set_visible(False)
        ax.yaxis.set_visible(False)
        ax.xaxis.set_tick_params(length=0)
        ax.xaxis.grid(True)
        ax.set_axisbelow(True)
        ax.patch.set_visible(False)

        # Add labels
        if self._total_labels:
            xlims = ax.get_xlim()
            xrange = abs(xlims[1] - xlims[0])
            xpadding = xrange * .05
            # Convert numbers to strings
            values = self.totals.values
            values_str = [self._label_format % i for i in values]
            for ii in range(len(values)):
                ax.text(values[ii] + xpadding, ii, values_str[ii], horizontalalignment='right', verticalalignment='center')

    def plot_shading(self, ax):
        # alternating row shading
        for i in range(0, len(self.totals), 2):
            rect = plt.Rectangle(self._swapaxes(0, i - .4),
                                 *self._swapaxes(*(1, .8)),
                                 facecolor='#f5f5f5', lw=0, zorder=0)
            ax.add_patch(rect)
        ax.set_frame_on(False)
        ax.tick_params(
            axis='both',
            which='both',
            left=False,
            right=False,
            bottom=False,
            top=False,
            labelbottom=False,
            labelleft=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    def plot(self, fig=None):
        """Draw all parts of the plot onto fig or a new figure

        Parameters
        ----------
        fig : matplotlib.figure.Figure, optional
            Defaults to a new figure.

        Returns
        -------
        subplots : dict of matplotlib.axes.Axes
            Keys are 'matrix', 'intersections', 'totals', 'shading'
        """
        if fig is None:
            fig = plt.figure(figsize=(10, 6))
        specs = self.make_grid(fig)
        shading_ax = fig.add_subplot(specs['shading'])
        self.plot_shading(shading_ax)
        matrix_ax = self._reorient(fig.add_subplot)(specs['matrix'], sharey=shading_ax)
        self.plot_matrix(matrix_ax)
        inters_ax = self._reorient(fig.add_subplot)(specs['intersections'], sharex=matrix_ax)
        self.plot_intersections(inters_ax)
        totals_ax = self._reorient(fig.add_subplot)(specs['totals'], sharey=matrix_ax)
        self.plot_totals(totals_ax)
        return {'fig': fig,
                'matrix': matrix_ax,
                'intersections': inters_ax,
                'shading': shading_ax,
                'totals': totals_ax}


def plot(data, fig=None, **kwargs):
    """Make an UpSet plot of data on fig

    Parameters
    ----------
    data : pandas.Series
        Values for each set to plot.
        Should have multi-index where each level is binary,
        corresponding to set membership.
    fig : matplotlib.figure.Figure, optional
        Defaults to a new figure.
    kwargs
        Other arguments for :class:`UpSet`

    Returns
    -------
    subplots : dict of matplotlib.axes.Axes
        Keys are 'matrix', 'intersections', 'totals', 'shading'
    """
    return UpSet(data, **kwargs).plot(fig)
