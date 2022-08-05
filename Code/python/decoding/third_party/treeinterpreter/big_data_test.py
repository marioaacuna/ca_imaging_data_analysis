import os
import pickle
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from decoding.third_party.treeinterpreter import treeinterpreter_old as ti


do_refit = False
clf_file = os.path.join(os.path.expanduser('~'), 'Downloads', 'RF.p')
data_file = os.path.join(os.path.expanduser('~'), 'Downloads', 'data.npy')


if not os.path.exists(data_file) or not os.path.exists(clf_file) or do_refit:
    # Generate large dataset
    n_samples = 50000
    n_features = 16000
    print('Generating dataset with %i features and %i observations' % (n_features, n_samples))
    X, y = make_classification(n_samples=n_samples, n_features=n_features, n_classes=2)

    # Split data in training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5)

    np.save(data_file, X_test)

    # Train RandomForest
    print('Training RandomForest')
    RF = RandomForestClassifier(n_estimators=100, n_jobs=-1)
    RF.fit(X_train, y_train)
    pickle.dump(RF, open(clf_file, 'wb'))

else:
    print('Loading dataset from disk')
    X_test = np.load(data_file)
    print('Reloading RandomForest')
    RF = pickle.load(open(clf_file, 'rb'))


# Compute feature contributions
print('Computing feature contribution')
results = ti.compute_feature_contributions_ensemble(RF, X_test, compute_conditional_contribution=False, n_jobs=None, verbose=True)


# Test equivalence between serial and parallel methods
# results_par = ti.compute_feature_contributions_ensemble(RF, X_test, compute_conditional_contribution=True, n_jobs=-1)
# k = list(results[3].keys())[0]
# np.array_equal(results[3][k], results_par[3][k])
# np.array_equal(results[0], results_par[0])
# np.array_equal(results[1], results_par[1])
# np.array_equal(results[2], results_par[2])

print('done')
