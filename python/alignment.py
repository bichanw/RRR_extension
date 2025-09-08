import numpy as np
# not tested!!!

def alignment_input(X, W, r=None, C=None):
    """
    Calculate how much the communication weights W align with 
    the principal components of the input X.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_input_neurons)
        Input matrix (stimuli).
    W : ndarray of shape (n_input_neurons, n_output_neurons)
        Communication weights.
    r : int, optional
        Rank of the communication weights W (if not provided,
        estimate from W).
    C : ndarray, optional
        Covariance matrix of X. If None, will compute cov(X).

    Returns
    -------
    aa : float
        Alignment index (0-1), where 1 is maximally aligned.
    p : None
        Placeholder (for compatibility with MATLAB signature).
    aa_rand : None
        Placeholder (for compatibility with MATLAB signature).
    """
    # Covariance of inputs
    if C is None:
        C = np.cov(X, rowvar=False)

    # PCA of covariance matrix
    _, Spcavec, _ = np.linalg.svd(C)

    # SVD of weights
    _, Swvec, _ = np.linalg.svd(W)
    # Pad singular values to length n_input_neurons
    Swvec_padded = np.concatenate([Swvec, np.zeros(W.shape[0] - len(Swvec))])

    # Compute alignment index
    amax = Spcavec @ (Swvec_padded ** 2)             # maximal value
    amin = Spcavec @ (np.flipud(Swvec_padded ** 2))  # minimal value
    araw = np.trace(W.T @ C @ W)                     # test statistic

    aa = (araw - amin) / (amax - amin)

    # Match MATLAB output signature
    return aa, None, None
