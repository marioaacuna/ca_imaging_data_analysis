from decoding.third_party.RandomForest2Embedding import RandomForest2Embedding
from sklearn.cluster import KMeans

rfc = RandomForest2Embedding()
default_clf_parameters = mcc.default_clf_parameters.copy()
default_clf_parameters.pop('oob_score')
default_clf_parameters.pop('class_weight')
rfc.set_params(**default_clf_parameters)
rfc.set_params(**mcc.optimized_parameters)
n_trees = data.shape[1] * 5
rfc.set_params(n_estimators=n_trees)
embedding = rfc.fit_transform(data)
data_shape = PARAMETERS['data_size'][class_names[class_idx]]

from sklearn.metrics import silhouette_score

max_n_clusters = int(data_shape[0])
scores = np.zeros(max_n_clusters + 1, dtype=float)
for n_clusters in range(2, max_n_clusters + 1):
    kmeans_RF = KMeans(init='k-means++', n_clusters=n_clusters,
                       random_state=1234, n_init=20, n_jobs=-1)
    clustered_labels = kmeans_RF.fit_predict(embedding)
    scores[n_clusters] = silhouette_score(embedding, clustered_labels)
plt.clf();
plt.plot(scores, '-ok')

n_clusters = np.argmax(scores)
km = KMeans(init='k-means++', n_clusters=n_clusters, random_state=1234,
            n_init=20, n_jobs=-1)
clustered_labels = km.fit_predict(embedding)
_, l = np.unique(clustered_labels, return_inverse=True)
clusters = pd.DataFrame(columns=['cluster_id', 'frequency'])
clusters['frequency'] = np.bincount(l)
clusters['cluster_id'] = np.arange(clusters.shape[0])
clusters.sort_values(by='frequency', ascending=False, inplace=True)
clusters.reset_index(drop=True, inplace=True)
# clusters2remove = clusters.loc[np.where(clusters['frequency'] < data_shape[0])[0], 'cluster_id'].values
# l = l.astype(np.float64); l[np.in1d(l, clusters2remove)] = np.nan; rel_heatmap = np.zeros(np.prod(data_shape), dtype=float) * np.nan; rel_heatmap[sample_idx_in_class] = l; rel_heatmap = rel_heatmap.reshape(n_trials_per_class[class_idx], -1); plt.classifier(); plt.imshow(rel_heatmap, aspect='auto', interpolation=None)
clusters2show = [6, 8];
_, l = np.unique(clustered_labels, return_inverse=True);
l = l.astype(np.float64);
l[np.logical_not(
    np.in1d(l, clusters.loc[clusters2show, 'cluster_id']))] = np.nan;
rel_heatmap = np.zeros(np.prod(data_shape), dtype=float) * np.nan;
rel_heatmap[sample_idx_in_class] = l;
rel_heatmap = rel_heatmap.reshape(n_trials_per_class[class_idx], -1);
plt.clf();
plt.imshow(rel_heatmap, aspect='auto', interpolation=None)

rel_heatmap = np.zeros(np.prod(data_shape), dtype=float) * np.nan;
rel_heatmap[sample_idx_in_class] = l;
rel_heatmap = rel_heatmap.reshape(n_trials_per_class[class_idx], -1);
plt.clf();
plt.imshow(rel_heatmap, aspect='auto', interpolation=None)














from Utilities.visualization import MidpointNormalize
from scipy.stats import sem
from sklearn.calibration import clone
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from decoding.third_party.Unsupervised_RF import build_synthetic_matrix, proximity_matrix
from decoding.third_party.RandomForest2Embedding import generate_discriminative_dataset
from decoding.third_party.forest_cluster.k_medoids import KMedoids
from decoding.third_party.Pycluster import kcluster

# For each neuron, get the contribution to each stimulus class after
# stimulus onset
reference_class = 'SP'
ref_class_idx = class_names.index(reference_class)
stim_class_idx = np.setdiff1d(range(n_classes), ref_class_idx)
n_stim_classes = len(stim_class_idx)
for i_cell in range(n_cells):
    # Get contributions to a stimulus class
    fig, ax = plt.subplots(nrows=4, ncols=n_classes,
                           figsize=(4 * n_classes, 12))

    for class_idx in range(n_classes):
        # contributions = feature_contributions_testing[LABELS_signif==class_idx, i_cell, class_idx]
        # Re-arrange values in an array of shape [trial x time]
        data_shape = PARAMETERS['data_size'][class_names[class_idx]]
        idx = SAMPLE_IDX.loc[(SAMPLE_IDX['significant_classified']) & (
                    SAMPLE_IDX[
                        'class'] == class_idx), 'sample_idx_in_class'].values
        # contribution_heatmap = np.zeros(np.prod(data_shape), dtype=float) #* np.nan
        # contribution_heatmap[idx] = contributions
        # contribution_heatmap = contribution_heatmap.reshape(n_trials_per_class[class_idx],-1)
        activity = np.zeros(np.prod(data_shape), dtype=float)
        activity[idx] = PARAMETERS['data'][class_names[class_idx]][idx, i_cell]
        activity = activity.reshape(n_trials_per_class[class_idx], -1)

        X = activity.copy().transpose()
        X_merged, y_merged = generate_discriminative_dataset(X,
                                                             method='bootstrap',
                                                             random_state=seed)
        X_merged = StandardScaler(with_mean=True, with_std=True).fit_transform(
            X_merged)

        clf = clone(mcc.base_estimator).named_steps['classifier']
        clf.set_params(n_estimators=1000)
        clf.fit(X_merged, y_merged)
        prox_mat = 1 - proximity_matrix(clf, X, normalize=True)

        # Iteratively assess number of clusters
        max_n_clusters = 8
        sil_score = np.zeros(max_n_clusters + 1, dtype=float) * np.nan
        for n_clusters in range(2, max_n_clusters + 1):
            # cluster_ids, error, n_found = kcluster(prox_mat, nclusters=n_clusters, method="a", npass=10)
            # sil_score[n_clusters] = silhouette_score(prox_mat, cluster_ids, random_state=seed)
            km = KMeans(n_clusters=n_clusters, init='k-means++', n_init=10,
                        random_state=seed, copy_x=True, n_jobs=-1).fit(prox_mat)
            sil_score[n_clusters] = silhouette_score(prox_mat, km.labels_,
                                                     random_state=seed)
        # Get optimal number of clusters from maximal silhouette score
        max_sil_score = np.nanmax(sil_score)
        n_clusters = np.where(sil_score == max_sil_score)[0][-1]
        print('Found %i clusters' % n_clusters)
        cluster_ids = KMeans(n_clusters=n_clusters, init='k-means++', n_init=10,
                             random_state=seed, copy_x=True, n_jobs=-1).fit(
            prox_mat).labels_

        # ax[0, class_idx].imshow(contribution_heatmap, aspect='auto', interpolation=None, cmap='RdBu_r', norm=MidpointNormalize(vmin=np.nanmin(contribution_heatmap), vmax=np.nanmax(contribution_heatmap), midpoint=0))
        # ax[1, class_idx].plot(np.nanmean(contribution_heatmap, axis=0), '-ok')
        ax[0, class_idx].imshow(prox_mat, aspect='auto', interpolation=None)
        ax[1, class_idx].plot(np.nanmean(prox_mat, axis=0),
                              '-ok')  # ; ax[1, class_idx].axhline(0, color='r')
        ax[2, class_idx].fill_between(x=np.arange(activity.shape[1]),
                                      y1=np.mean(activity, axis=0) - sem(
                                          activity, axis=0),
                                      y2=np.mean(activity, axis=0) + sem(
                                          activity, axis=0),
                                      facecolor=(0.04, 0.78, 1), alpha=0.5)
        ax[2, class_idx].plot(np.mean(activity, axis=0), '-ok', markersize=3)
        # ax[3, class_idx].plot(PARAMETERS['data'][class_names[class_idx]][idx, i_cell], contributions, 'ok'); ax[3, class_idx].axhline(0, color='r')
        ax[3, class_idx].imshow(np.atleast_2d(cluster_ids), aspect='auto',
                                interpolation=None, cmap='Dark2')

    plt.tight_layout()









