import warnings
warnings.simplefilter('ignore')
import os
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.collections import PathCollection
from Utilities.visualization import make_seaborn_palette
from scipy.stats import ttest_ind, ranksums


# analysis_type = 'all'
# analysis_type = 'selective'
# analysis_type = 'stable'
# analysis_type = 'ensemble'
analysis_type = 'specific'


do_print_figures = False
do_statistic_tests = 'bootstrap'  # None, 'ttest', 'wilcoxon' and 'bootstrap'
is_bootstrap_test_parametric = True
n_bootstrap_iterations = 1000
correct_p_values_sampling_error = True
show_stars = True


# Load data
filename = r'D:\_MATLAB_CaImaging\%s_cells.csv' % analysis_type
data = pd.read_csv(filename)
np.random.seed(17)

# Set stimulus order and color palette
stimulus_order = ('HPS', 'FPS', 'temp_48', 'temp_43', 'temp_38', 'pinprick', 'puff', 'touch', 'sound', 'odor_plus')
palette = make_seaborn_palette([(255, 0, 4), (128, 128, 128)])
jitter = .3
mean_range = jitter / 2


def make_graph(data, title_str, pdf_filename):
    """

    :param data:
    :param title_str:
    :param pdf_filename:
    :return:
    """
    p_values = None

    # Test sham and CCI
    if do_statistic_tests is not None:
        # Define bootstrap test function
        def bootstrap_test(x, y):
            n_x = x.shape[0]
            n_y = y.shape[0]
            if n_x < 2 or n_y < 2:
                return None, 1

            # Compute empirical differences
            empirical_difference = bootstrap_diff_function(x) - bootstrap_diff_function(y)

            # Set H0 as true (i.e. mean is 0)
            xy = np.hstack((x, y))
            gm = bootstrap_diff_function(xy)
            xy = np.hstack((x - bootstrap_diff_function(x), y - bootstrap_diff_function(y))) + gm

            n_samples = xy.shape[0]
            bootstrap_difference = np.zeros((n_bootstrap_iterations, ), dtype=float)
            for i_iter in range(n_bootstrap_iterations):
                idx_x = np.random.randint(0, n_samples - 1, size=(n_x, ))
                idx_y = np.random.randint(0, n_samples - 1, size=(n_y, ))
                bootstrap_difference[i_iter] = bootstrap_diff_function(xy[idx_x]) - bootstrap_diff_function(xy[idx_y])

            # Compute p-value
            p = (np.sum(np.abs(bootstrap_difference) >= np.abs(empirical_difference)) + 1) / (n_bootstrap_iterations + 1)
            # p has to be the second output
            return None, p

        if do_statistic_tests == 'ttest':
            test_function = ttest_ind
        elif do_statistic_tests == 'wilcoxon':
            test_function = ranksums
        elif do_statistic_tests == 'bootstrap':
            test_function = bootstrap_test
            if is_bootstrap_test_parametric:
                bootstrap_diff_function = np.nanmean
            else:
                bootstrap_diff_function = np.nanmedian

        # Define function to select performance in each group so that groups can be tested against each other
        def select(df, group):
            idx = np.where(df['group'] == group)[0]
            performance = df.reset_index(drop=True).loc[idx, 'performance']
            return performance

        p_values = this_data.groupby(by=['response_stimulus', 'session_idx']).apply(lambda df: test_function(select(df, 'sham'), select(df, 'CCI'))[1])
        p_values = pd.DataFrame(p_values, columns=['p'])
        p_values['sig'] = p_values['p'] < 0.05
        p_values.reset_index(inplace=True)

        # If there is a difference at baseline, use that p-value as the new alpha
        if correct_p_values_sampling_error:
            def set_new_alpha(df):
                # Take p-values in sessions before surgery
                if analysis_type != 'stable':
                    p = df.iloc[np.where(df['session_idx'] < 0)[0]]['p'].values.min()
                else:
                    p = df.iloc[np.where(np.in1d(df['session_idx'], 'before'))[0]]['p'].values.min()
                if p < 0.05:
                    df['sig'] = df['p'] < p
                return df

            # Apply correction
            p_values = p_values.groupby(by='response_stimulus').apply(lambda df: set_new_alpha(df))

    # Draw scatterplots
    g = sns.catplot(kind='strip', data=data, x='session_idx', y='performance',
                    hue='group', col='response_stimulus',
                    hue_order=['CCI', 'sham'], col_order=stimulus_order,
                    col_wrap=4, jitter=jitter, dodge=True, orient='v',
                    edgecolor='w', linewidth=1,
                    height=2.5, aspect=.8, palette=palette)

    for ax in g.axes:
        ch = ax.get_children()
        pc = [i for i in ch if isinstance(i, PathCollection)]
        if len(pc) == 0:
            continue
        for pc_i in pc:
            # Increase scatter size
            pc_i.set_sizes([50])

            # Draw mean
            offsets = pc_i.get_offsets()
            y = np.mean(offsets[:, 1])
            x = np.median(offsets[:, 0])
            ax.plot([x - mean_range, x + mean_range], [y, y], color='k', lw=2, zorder=100)

    # Draw line at 0.5
    [i.axhline(.5, color=(.7, .7, .7), lw=.5) for i in g.axes]

    # Fix axes titles
    for ax in g.axes:
        old_title = ax.title.get_text()
        new_title = old_title.split('=')[1].strip()
        ax.set_title(new_title, fontsize=14, fontweight='bold')

    # Set figure title and subplot spacing
    g.fig.suptitle(title_str, fontsize=18)
    g.fig.set_size_inches(15, 10, forward=True)
    g.fig.subplots_adjust(left=.05, top=.9, right=.92, wspace=.1)

    # Add stars for significance
    if show_stars and do_statistic_tests is not None:
        significant_differences = np.where(p_values['sig'])[0]
        n_stars = significant_differences.shape[0]
        for i_star in range(n_stars):
            info = p_values.loc[significant_differences[i_star]]
            ax_idx = stimulus_order.index(info['response_stimulus'])
            session_idx = info['session_idx']
            x = pd.unique(data['session_idx']).tolist().index(session_idx)
            y = max(g.axes[ax_idx].get_ylim())
            g.axes[ax_idx].text(x, y, '*', ha='center', va='top', fontsize=18, fontweight='bold')

    # Print figure
    if do_print_figures:
        print('Figure printed in %s' % pdf_filename)
        g.fig.savefig(pdf_filename)

    if do_statistic_tests is not None and not do_print_figures:
        print(title_str)
        print(p_values)


if __name__ == '__main__':
    # Filter data
    if analysis_type == 'all':
        print('Analyzing all cells')

        rows = np.in1d(data['session_idx'], [-2, -1, 1, 2])
        this_data = data.loc[np.where(rows)[0]].dropna().reset_index(drop=True)
        title_str = 'All cells'
        pdf_filename = os.path.splitext(filename)[0] + '.pdf'
        make_graph(this_data, title_str, pdf_filename)

    elif analysis_type == 'selective':
        stimuli_selectivity = pd.unique(data['selective_to'])

        for stimulus_selective in stimuli_selectivity:
            print('Analyzing cells selective to %s' % stimulus_selective)

            rows = np.in1d(data['session_idx'], [-2, -1, 1, 2]) & \
                   np.in1d(data['selective_to'], stimulus_selective) & \
                   (data['selective_when'] == data['session_idx'])
            this_data = data.loc[np.where(rows)[0]].dropna().reset_index(drop=True)
            title_str = 'Cells selective to %s' % stimulus_selective
            pdf_filename = os.path.splitext(filename)[0] + '_%s.pdf' % stimulus_selective
            make_graph(this_data, title_str, pdf_filename)

    elif analysis_type == 'stable':
        stimuli_selectivity = pd.unique(data['selective_to'])

        for stimulus_selective in stimuli_selectivity:
            print('Analyzing cells stably selective to %s' % stimulus_selective)

            rows = np.in1d(data['selective_to'], stimulus_selective)
            this_data = data.loc[np.where(rows)[0]].dropna().reset_index(drop=True)
            title_str = 'Cells stably selective to %s' % stimulus_selective
            pdf_filename = os.path.splitext(filename)[0] + '_%s.pdf' % stimulus_selective
            make_graph(this_data, title_str, pdf_filename)

    elif analysis_type == 'ensemble':
        print('Analyzing ensemble cells')

        rows = np.in1d(data['session_idx'], [-2, -1, 1, 2])
        this_data = data.loc[np.where(rows)[0]].dropna().reset_index(drop=True)
        title_str = 'Ensemble cells'
        pdf_filename = os.path.splitext(filename)[0] + '.pdf'
        make_graph(this_data, title_str, pdf_filename)

    elif analysis_type == 'specific':
        stimuli_selectivity = pd.unique(data['specific_to'])

        for stimulus_selective in stimuli_selectivity:
            print('Analyzing cells specific to %s' % stimulus_selective)

            rows = np.in1d(data['session_idx'], [-2, -1, 1, 2]) & \
                   np.in1d(data['specific_to'], stimulus_selective) & \
                   (data['specific_when'] == data['session_idx'])
            this_data = data.loc[np.where(rows)[0]].dropna().reset_index(drop=True)
            title_str = 'Cells specific to %s' % stimulus_selective
            pdf_filename = os.path.splitext(filename)[0] + '_%s.pdf' % stimulus_selective
            make_graph(this_data, title_str, pdf_filename)

    # Show figures if they won't get printed
    if not do_print_figures:
        plt.show()

