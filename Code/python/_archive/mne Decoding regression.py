import pandas as pd
import mne
from mne.stats import linear_regression, fdr_correction
from mne.viz import plot_compare_evokeds
from mne.datasets import kiloword

# Load the data
path = kiloword.data_path() + '/kword_metadata-epo.fif'
epochs = mne.read_epochs(path)
print(epochs.metadata.head())

name = "Concreteness"
df = epochs.metadata
df[name] = pd.cut(df[name], 11, labels=False) / 10
colors = {str(val): val for val in df[name].unique()}
epochs.metadata = df.assign(Intercept=1)  # Add an intercept for later
evokeds = {val: epochs[name + " == " + val].average() for val in colors}
# plot_compare_evokeds(evokeds, colors=colors, split_legend=True,
#                      cmap=(name + " Percentile", "viridis"))


names = ["Intercept", name]
res = linear_regression(epochs, epochs.metadata[names], names=names)
# for cond in names:
#     res[cond].beta.plot_joint(title=cond, ts_args=dict(time_unit='s'),
#                               topomap_args=dict(time_unit='s'))

reject_H0, fdr_pvals = fdr_correction(res["Concreteness"].p_val.data)
evoked = res["Concreteness"].beta
evoked.plot_image(mask=reject_H0, time_unit='s')

