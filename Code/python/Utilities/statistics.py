import numpy as np
import scipy as sp
import scipy.stats


def confidence_interval(data, confidence=0.95, parametric=True):
    data = data.astype(np.float64)
    a = data[np.logical_not(np.isnan(data))]
    n = a.shape[0]

    if parametric:
        m = np.nanmean(a, axis=0)
        se = scipy.stats.sem(a, axis=0, nan_policy="omit")
        h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
        ci = (m-h, m+h)
    else:
        ci = np.percentile(a, q=[np.round((1-confidence)*100), np.round(confidence*100)], axis=0)

    return ci


def cohen_d(x, y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (x.mean() - y.mean()) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)

