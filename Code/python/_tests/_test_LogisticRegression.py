from Utilities import matlab_file
import numpy as np
import warnings
warnings.filterwarnings('ignore')

filename = r'D:\_MATLAB_2PI\FK_5\FK_5_cond1_decoding_data.mat'
PARAMETERS = matlab_file.load(filename)

class_names = list(PARAMETERS['data'].keys())
n_classes = len(class_names)

# Make sure that data is 2D
for cl in class_names:
    if PARAMETERS['data'][cl].ndim == 1:
        PARAMETERS['data'][cl] = np.atleast_2d(
                PARAMETERS['data'][cl]).transpose()

# Prepare data and labels for multi-class classifier
LABELS = []  # array of labels
for ci, cl in enumerate(class_names):
    n_samples = PARAMETERS['data'][cl].shape[0]
    these_labels = np.zeros((n_samples,), dtype=int) + ci
    LABELS.append(these_labels)
# Concatenate the labels array, and make data array
y = np.hstack(LABELS)
X = np.vstack(([PARAMETERS['data'][cl] for cl in class_names]))


from sklearn.linear_model import LogisticRegression, LinearRegression, SGDClassifier, LogisticRegressionCV
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_predict
import pycm
from sklearn.naive_bayes import ComplementNB
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.multiclass import OneVsOneClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
import lightgbm as lgb
import xgboost as xgb

random_state = 2286


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, stratify=y, random_state=random_state)

print()



# nb = ComplementNB().fit(X_train, y_train)
# y_proba = nb.predict_proba(X_test)
# print('%s: ComplementNB' % ['%.2f' % i if i is not 'None' else '0' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])

# ss = StandardScaler(copy=True, with_std=True, with_mean=True).fit(X_train)
ss = MinMaxScaler(copy=True, feature_range=(0, 1)).fit(X_train)
X_train = ss.transform(X_train)
X_test = ss.transform(X_test)


print()
# train_data = lgb.Dataset(X_train, label=y_train)
# param = dict(num_class=4, objective='multiclass', boosting='gbdt', verbosity=0, tree_learner='feature')
# num_round = 10
# bst = lgb.train(param, train_data, num_round, verbose_eval=False)
# y_proba = bst.predict(X_test, num_iteration=-1)
# print('%s: LightGBM' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])


# model = xgb.XGBClassifier(random_state=random_state, learning_rate=0.01, objective='multi:softprob', n_jobs=-1, verbosity=0)
# # model.fit(X_train, y_train)
# # y_proba = model.predict_proba(X_test)
# # print('%s: XGBClassifier' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])
# kfold = StratifiedKFold(n_splits=5, random_state=random_state)
# results = cross_val_predict(model, X, y, cv=kfold)
# print('%s: XGBClassifier' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y, results).AUPR.values())])

# lr = LogisticRegression(solver='lbfgs', penalty='l2', n_jobs=-1, fit_intercept=True, max_iter=10000, multi_class='multinomial', random_state=random_state)
# lr.fit(X_train, y_train)
# y_proba = lr.predict_proba(X_test)
# print('%s: LogisticRegression (lbfgs l2)' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])
#
# lr = LogisticRegression(solver='newton-cg', penalty='l2', n_jobs=-1, fit_intercept=True, max_iter=10000, multi_class='multinomial', random_state=random_state)
# lr.fit(X_train, y_train)
# y_proba = lr.predict_proba(X_test)
# print('%s: LogisticRegression (newton-cg l2)' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])
#
# lr = LogisticRegression(solver='liblinear', penalty='l2', n_jobs=-1, fit_intercept=True, max_iter=10000, random_state=random_state)
# lr.fit(X_train, y_train)
# y_proba = lr.predict_proba(X_test)
# print('%s: LogisticRegression (liblinear)' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])

# lr = LogisticRegression(solver='saga', penalty='l2', n_jobs=-1, fit_intercept=True, max_iter=10000, multi_class='multinomial', random_state=random_state, C=10)
# lr.fit(X_train, y_train)
# y_proba = lr.predict_proba(X_test)
# print('%s: LogisticRegression (saga l2)' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])
# print(lr.get_params())

# lr = LogisticRegression(solver='saga', penalty='l1', n_jobs=-1, fit_intercept=True, max_iter=10000, multi_class='multinomial', random_state=random_state)
# lr.fit(X_train, y_train)
# y_proba = lr.predict_proba(X_test)
# print('%s: LogisticRegression (saga l1)' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])


lr = LogisticRegressionCV(solver='saga', penalty='l2', n_jobs=-1, fit_intercept=True, max_iter=10000, multi_class='multinomial', cv=5, random_state=random_state)
lr.fit(X_train, y_train)
y_proba = lr.predict_proba(X_test)
print('%s: LogisticRegressionCV (saga l2)' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])
# print(lr.get_params())

# lr = LogisticRegressionCV(solver='saga', penalty='l1', n_jobs=-1, fit_intercept=True, max_iter=10000, multi_class='multinomial', cv=5, random_state=random_state)
# lr.fit(X_train, y_train)
# y_proba = lr.predict_proba(X_test)
# print('%s: LogisticRegressionCV (saga l1)' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])

# clf = SGDClassifier(loss='log', penalty='l1', fit_intercept=True, shuffle=True, n_jobs=-1, learning_rate='optimal', average=True, random_state=random_state)
# clf.fit(X_train, y_train)
# y_proba = clf.predict_proba(X_test)
# print('%s: SGDClassifier' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])

clf = SGDClassifier(loss='log', penalty='l1', fit_intercept=True, shuffle=True, n_jobs=-1, learning_rate='optimal', average=True, random_state=random_state)
ovo = OneVsOneClassifier(clf).fit(X_train, y_train)
print('%s: LogisticRegression SGD OVO' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, ovo.predict(X_test)).AUPR.values())])


# rf = RandomForestClassifier(n_estimators=128, criterion='entropy', max_features=int(np.floor(np.log(X_train.shape[1]))), bootstrap=True, n_jobs=-1, random_state=random_state)
# rf.fit(X_train, y_train)
# y_proba = rf.predict_proba(X_test)
# print('%s: RandomForestClassifier' % ['%.2f' % i if i is not 'None' else '    ' for i in list(pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).AUPR.values())])




print()
# pycm.ConfusionMatrix(y_test, np.argmax(y_proba, 1)).print_matrix()
