import numpy as np


def raised_cosine_filter(max_time, n_bins, k, nBases=15):
    """Return a cosine bump, kth of nBases, such that the bases tile
    the interval [0, n_bins].

    To plot these bases:
    for i in range(10):
        b =  raised_cosine_filter(250, i, nBases = 10)
        plt.plot(b)
    """
    assert all([isinstance(p, int) for p in [n_bins, k, nBases]])

    t = np.linspace(0, max_time, n_bins)

    nt = np.log(t + 0.1)

    cSt, cEnd = nt[1], nt[-1]
    db = (cEnd - cSt) / (nBases)
    c = np.arange(cSt, cEnd, db)

    bas = np.zeros((nBases, t.shape[0]))

    filt = lambda x: (np.cos(np.maximum(-np.pi, np.minimum(np.pi, (nt - c[x ] )*np.pi/(db))) ) + 1) / 2

    this_filt = filt(k)

    return this_filt/np.nansum(this_filt)
