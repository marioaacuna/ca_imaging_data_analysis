# Numerical packages
import numpy as np
from scipy.spatial.distance import pdist
from sklearn.svm import OneClassSVM
from sklearn.preprocessing import StandardScaler


def train_and_test_model(data, labels, sample_weights, nu, return_metric='decision_function'):
    # Initialize variables
    performance = None

    # Normalize data
    scl = StandardScaler(with_mean=True, with_std=True, copy=True).fit(data[labels == -1, :])
    data = scl.transform(data)
    # Train classifier
    distances = pdist(data, metric='euclidean')
    gamma = (2 ** 4) / (distances.mean() ** 2)
    svm = OneClassSVM(nu=nu, kernel='rbf', gamma=gamma)
    svm.fit(data[labels == -1, :], sample_weight=sample_weights)
    # Compute performance
    if return_metric == 'decision_function':  # Distance from hyperplane
        performance = svm.decision_function(data)

    return performance


def do_significance_test_permute_labels(data, labels, sample_weights, nu):
    # OLD VERSION #
    # Permute labels with replacement across the entire dataset, and keep class size
    # (http://mvpa.blogspot.com/2013/06/mvpa-permutation-schemes-permutation.html)
    # new_data, _ = resample(data, labels, replace=True)
    # Compute accuracy
    # return train_and_test_model(new_data, labels, sample_weights, nu, return_metric='decision_function')

    # NEW VERSION #
    n_samples = data.shape[0]
    # Permute weights
    bootstrapped_weights = np.random.choice(sample_weights, n_samples, replace=False)
    bootstrapped_weights /= bootstrapped_weights.sum()  # Make sure they sum up to 1
    # Compute accuracy
    return train_and_test_model(data, labels, bootstrapped_weights, nu, return_metric='decision_function')


def compute_bootstrap_p_value(null_distribution, empirical_value, sign='lt'):
    # Find where empirical value is more extreme than the null distribution
    if sign == 'lt':
        difference = np.array(null_distribution <= empirical_value, dtype=int)
    else:
        difference = np.array(null_distribution >= empirical_value, dtype=int)
    # Count the number of times the relationship is true
    count_per_resample = np.sum(difference, axis=0)
    # Estimate a one-sided bootstrap p-value
    p = (count_per_resample + 1) / (difference.shape[0] + 1)

    return p

