# System packages
import os, warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

# Numerical packages
import numpy as np
from sklearn.neighbors.kde import KernelDensity
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

# Visualization packages
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Parallel computation
from joblib import Parallel, delayed

# Local repository
from Utilities.visualization import adjust_spines, MidpointNormalize
from Utilities.IO_operations import log


def compute_boosted_KDE(data, k_iterations, PARAMETERS_output_figures=None, class_name='', overwrite_figures=False, debug=False):
    # Initialize variables
    weights = None

    # Normalize data to unit standard-deviation
    data_norm = StandardScaler(with_mean=True, with_std=True, copy=True).fit_transform(data)

    # Get data dimensionality
    m, d = data_norm.shape

    # Allocate weights (initialized as uniform distribution), and arrays to keep
    # track of PDF and leave-one-out PDF at each boosting iteration
    if k_iterations == 0:
        return np.ones((m, ), dtype=np.float32) / m
    # Continue normally
    w = np.ones((m, k_iterations), dtype=np.float32) / m
    pdf_kde = np.zeros((m, k_iterations), dtype=np.float32) * np.nan
    loo_kde = np.zeros((m, k_iterations), dtype=np.float32) * np.nan

    # Compute bandwidth according to Silverman's rule
    bw = (m * (d + 2) / 4.) ** (-1. / (d + 4))

    # Boost KDE with selected bandwidth
    for k in range(k_iterations):
        # Adjust weights in all iterations after the first one
        if k > 0:
            # In the original paper, the log of the ratio was considered, but
            # here both PDFs are already logs. Therefore, apply transformation
            # according to which log(x/y) = log(x) - log(y)
            w[:, k] = w[:, k-1] + (pdf_kde[:, k - 1] - loo_kde[:, k - 1])
            # Make sure weights sum up to 1
            w[:, k] /= np.sum(w[:, k])
        # Reshape weights
        weights = np.atleast_2d(w[:, k]).transpose()

        # Compute PDF
        pdf_kde[:, k] = compute_kde(data_norm, bw, weights)
        # Compute leave-one-out PDF
        loo_kde[:, k] = compute_loo_kde(data_norm, bw, weights).ravel()

        # Do final adjustment of weights
        if k == k_iterations - 1:
            weights = w[:, -1] + (pdf_kde[:, -1] - loo_kde[:, -1])
            weights /= np.sum(weights)

    # Store weights of last iteration
    final_weights = weights.copy()
    # Add final weights to main variable w
    w = np.hstack((w, np.atleast_2d(final_weights).T))

    ### OLD ###
    # Invert and normalize the final probability weights (this creates
    # asymmetries in the final distribution)
    # boosted_weights = 1. / final_weights

    ### NEW ###
    boosted_weights = 2 * (1 / m) - final_weights
    # Flip weights around value of uniform distribution
    # Cap negative weights to 0
    neg_idx = np.where(boosted_weights < 0)[0]
    boosted_weights[neg_idx] = 1. / final_weights[neg_idx] / np.sum(1. / final_weights)
    boosted_weights[boosted_weights < 0] = 0
    # Make sure weights sum up to 1
    boosted_weights /= np.sum(boosted_weights)

    # Make diagnostic plots on boosting
    if PARAMETERS_output_figures['outlier_detection'] is not None:
        # Sort data according to sum of all neurons
        data_sum = np.sum(data_norm, axis=1)
        idx = np.argsort(data_sum)
        x = data_sum[idx]

        # Make folder if it doesn't exist
        if not os.path.exists(PARAMETERS_output_figures['folder']):
            os.mkdir(PARAMETERS_output_figures['folder'])

        if 'KDE_weights' in PARAMETERS_output_figures['outlier_detection']:
            output_filename = os.path.join(PARAMETERS_output_figures['folder'], PARAMETERS_output_figures['base_filename'] + class_name + '_KDE_weights_normalization.pdf')
            if not os.path.exists(output_filename) or overwrite_figures:
                if debug:
                    log('Printing plot of KDE weights normalization')
                fig = plt.figure(figsize=(10, 10))
                plt.clf()
                ax = fig.add_subplot(1, 1, 1)
                ax.plot(x, w[idx, -1], 'o', label='original')
                ax.plot(x, 1 / w[idx, -1] / np.sum(1 / w[idx, -1]), 'o', label='1/n')
                ax.plot(x, boosted_weights[idx], 'o', label='inverted')
                ax.axhline(w[0, 0], color='r', label='uniform')
                ax.legend()
                plt.grid(True, axis='both')
                y_lims = [0, np.quantile(w[:, -1], .975)]
                y_lim_range_padding = (y_lims[1] - y_lims[0]) * .05
                ax.set_ylim(y_lims[0] - y_lim_range_padding, y_lims[1] + y_lim_range_padding)
                ax.set_title('Probability weights')
                ax.set_xlabel('Sum of activity')
                ax.set_ylabel('Probability weight')
                adjust_spines(ax, ['bottom', 'left'], smart_bounds=True)
                plt.tight_layout()
                # Save figure
                plt.savefig(fname=output_filename, bbox_inches='tight', format='pdf')
                plt.close(fig)

            output_filename = os.path.join(PARAMETERS_output_figures['folder'], PARAMETERS_output_figures['base_filename'] + class_name + '_KDE_weights_boosting.pdf')
            if not os.path.exists(output_filename) or overwrite_figures:
                if debug:
                    log('Printing plot of KDE weights boosting')
                # Make another figure that shows how weights were normalized at each
                # iteration and shows a raster plot of what weights will influence.
                data_to_show = np.hstack((np.atleast_2d(np.sum(data_norm, axis=1)).T, data_norm))
                idx = np.argsort(data_to_show[:, 0])[::-1]
                # Show adjustment of weights
                w_to_show = np.log(np.hstack((w, np.atleast_2d(boosted_weights).T)))
                w_to_show[np.logical_not(np.isfinite(w_to_show))] = w_to_show[np.isfinite(w_to_show)].min()
                fig = plt.figure(figsize=(25, 11))
                plt.clf()
                ax1 = fig.add_subplot(1, 2, 1)
                im = ax1.imshow(w_to_show[idx, :], aspect='auto', cmap='RdBu_r', norm=MidpointNormalize(vmin=w_to_show.min(), vmax=w_to_show.max(), midpoint=w_to_show[0, 0]))
                cbar = fig.colorbar(im, ax=ax1)
                cbar.set_label('Log probability weight')
                ax1.set_xticks(range(k_iterations+2)); ax1.set_xticklabels(['uniform'] + list(range(1, k_iterations+1)) + ['final'])
                ax1.set_xlabel('Iteration #')
                ax1.set_ylabel('Sorted bins')
                ax1.set_title('Weight adjustment')
                # Show corresponding neural activity
                ax2 = fig.add_subplot(1, 2, 2)
                im = plt.imshow(data_to_show[idx, :], aspect='auto', cmap='RdBu_r', norm=MidpointNormalize(vmin=data_norm.min(), vmax=data_norm.max(), midpoint=0))
                cbar = fig.colorbar(im, ax=ax2)
                cbar.set_label('Event amplitude (z-score)')
                ax2.set_xticks(np.arange(data.shape[1]+1))
                xticks = np.arange(0, data.shape[1]+1)[::5][1:]
                xticklabels = np.empty(shape=data.shape[1], dtype=object)
                xticklabels[xticks - 1] = np.array(['%i' % i for i in xticks])
                xticklabels = [i if i is not None else '' for i in xticklabels]
                xticklabels.insert(0, 'sum')
                ax2.set_xticklabels(xticklabels)
                ax2.set_xlabel('Neuron #')
                ax2.get_yaxis().set_visible(False)
                ax2.set_title('Neural activity')
                # Adjust layout
                plt.tight_layout()
                plt.savefig(fname=output_filename, bbox_inches='tight', format='pdf')
                plt.close(fig)

        if 'PDF_all_neurons' in PARAMETERS_output_figures['outlier_detection']:
            # Create pdf file
            output_filename = os.path.join(PARAMETERS_output_figures['folder'], PARAMETERS_output_figures['base_filename'] + class_name + '_KDE_weighted_PDFs.pdf')
            if not os.path.exists(output_filename) or overwrite_figures:
                if debug:
                    log('Printing probability distribution function (PDF) of all neurons')
                # Open pdf file
                PDF_file = PdfPages(output_filename)
                for neuron_idx in range(d):
                    # Get data
                    data_to_show_original = data[:, neuron_idx].reshape(-1, 1)
                    data_bootstrapped = data[np.random.choice(m, m, p=boosted_weights, replace=True), neuron_idx].reshape(-1, 1)
                    pdf_limits = np.hstack((data_to_show_original.ravel(), data_bootstrapped.ravel()))
                    pdf_limits = [pdf_limits.min(), pdf_limits.max()]
                    pdf_grid = np.linspace(pdf_limits[0], pdf_limits[1], 100).reshape(-1, 1)
                    # Estimate unweighted PDF and weighted PDF
                    pdf = np.exp(compute_kde(data_to_show_original, bw, np.ones((m, k_iterations), dtype=np.float32) / m, return_pdf_at=pdf_grid))
                    pdf_bs = np.exp(compute_kde(data_to_show_original, bw, boosted_weights, return_pdf_at=pdf_grid))

                    # Plot PDFs
                    fig = plt.figure(figsize=(10, 5))
                    plt.clf()
                    ax = fig.add_subplot(1, 2, 1)
                    ax.plot(pdf_grid, pdf, '-k', label='original')
                    ax.plot(pdf_grid, pdf_bs, '-r', label='bootstrap')
                    ax.legend()
                    ax2 = fig.add_subplot(1, 2, 2)
                    ax2.plot(pdf_grid, pdf, '-k', label='original')
                    ax2.plot(pdf_grid, pdf_bs, '-r', label='bootstrap')
                    ax.set_xlabel(r'Event amplitude ($\it{a.u.}$)')
                    ax2.set_xlabel('Event amplitude ($\it{a.u.}$)')
                    ax.set_title('Probability distribution function')
                    ax2.set_title('Log probability distribution function')
                    # Adjust layout
                    ax.set_xlim(pdf_grid.min(), pdf_grid.max())
                    ax2.set_xlim(pdf_grid.min(), pdf_grid.max())
                    ax2.set_yscale('log')
                    adjust_spines(ax, ['bottom', 'left'], smart_bounds=True)
                    adjust_spines(ax2, ['bottom', 'left'], smart_bounds=True)
                    fig.text(0.01, 0.99, 'Cell %i' % (neuron_idx + 1), horizontalalignment='left', verticalalignment='top', fontsize=18)
                    plt.tight_layout()
                    # Print figure to file
                    PDF_file.savefig(fig, bbox_inches='tight')
                    plt.close(fig)

                # Flush images to disk
                PDF_file.close()

    # End function by returning weights
    return boosted_weights, final_weights


def compute_loo_kde(data, bw, weights):
    # Compute leave-one-out KDE
    return np.array(Parallel(n_jobs=-1)(
            delayed(compute_kde)(data[train_index, :], bw, weights[train_index, :], data[test_index, :])
            for train_index, test_index in LeaveOneOut().split(data)))


def compute_kde(data, bw, weights, return_pdf_at=None):
    # Compute KDE
    # Fit KDE on training data
    kde = KernelDensity(bandwidth=bw, kernel='gaussian', metric='euclidean', algorithm='ball_tree', breadth_first=True, leaf_size=40).fit(data, sample_weight=weights.ravel())
    # Set data points at which to return the PDF
    if return_pdf_at is None:
        return_pdf_at = data
    # Return PDF at test points
    return kde.score_samples(return_pdf_at)

