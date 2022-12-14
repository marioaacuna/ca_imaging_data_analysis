import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm


def plot_decision_function(classifier, sample_weight, axis, title):
    # plot the decision function
    xx, yy = np.meshgrid(np.linspace(-4, 5, 500), np.linspace(-4, 5, 500))

    Z = classifier.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)

    # plot the line, the points, and the nearest vectors to the plane
    axis.contourf(xx, yy, Z, alpha=0.75, cmap=plt.cm.bone)
    axis.scatter(X[:, 0], X[:, 1], c='w', s=150 * sample_weight, alpha=0.9,
                 cmap=plt.cm.bone, edgecolors='black')

    axis.axis('off')
    axis.set_aspect('equal')
    axis.set_title(title)


# we create 20 points
np.random.seed(0)
X = np.r_[np.random.randn(10, 2) + [1, 1], np.random.randn(10, 2)]
sample_weight_last_ten = abs(np.random.randn(len(X)))
sample_weight_constant = np.ones(len(X))
# and bigger weights to some outliers
sample_weight_last_ten[:5] *= 3
sample_weight_last_ten[9] *= 5

# for reference, first fit without sample weights
percent_outliers = 5 / X.shape[0]
nu = 1 - percent_outliers

# fit the model
clf_weights = svm.OneClassSVM(nu=nu, kernel='rbf', gamma='scale')
clf_weights.fit(X, sample_weight=sample_weight_last_ten)

clf_no_weights = svm.OneClassSVM(nu=nu, kernel='rbf', gamma='scale')
clf_no_weights.fit(X)

fig, axes = plt.subplots(1, 2, figsize=(17, 9))
plot_decision_function(clf_no_weights, sample_weight_constant, axes[0],
                       "Constant weights")
plot_decision_function(clf_weights, sample_weight_last_ten, axes[1],
                       "Modified weights")

plt.tight_layout()
plt.show()
