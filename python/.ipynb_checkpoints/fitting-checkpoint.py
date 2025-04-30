import numpy as np
from scipy.linalg import sqrtm, inv

def svd_RRR(X, Y, rnk, lambda_=0):
    # Ridge regularization
    XX = X.T @ X + lambda_ * np.eye(X.shape[1])

    # Least squares estimate with ridge
    if np.linalg.cond(XX) < 1e10:
        wridge = np.linalg.solve(XX, X.T @ Y)
    else:
        wridge = np.linalg.pinv(XX) @ (X.T @ Y)

    # SVD of relevant matrix
    _, _, vrrr = np.linalg.svd(Y.T @ X @ wridge)

    # Get the top 'rnk' components
    vrrr = vrrr[:rnk, :].T   # shape: (features, rnk)
    urrr = wridge @ vrrr    # shape: (features, rnk)

    # Construct full RRR estimate
    w0 = urrr @ vrrr.T
    vrrr = vrrr.T  # for compatibility with original code's return

    return w0, urrr, vrrr


def svd_RRR_noniso(X, Y, rnk, C=None):
    # Least squares estimate
    wls = np.linalg.solve(X.T @ X, X.T @ Y)

    # Compute covariance of residuals if C is not provided
    if C is None:
        res_wls = Y - X @ wls
        C = (res_wls.T @ res_wls) / (X.shape[0] - 1)

    # Compute inverse sqrt and sqrt of C
    C_inv_sqrt = inv(sqrtm(C))
    C_sqrt = sqrtm(C)

    # SVD of whitened cross-covariance
    _, _, vrrr = np.linalg.svd(C_inv_sqrt @ Y.T @ X @ wls @ C_inv_sqrt)

    # Adjust for non-isotropic noise
    vrrr = C_sqrt @ vrrr[:, :rnk]
    urrr = np.linalg.solve(X.T @ X, X.T @ Y @ inv(C) @ vrrr)

    # Reconstruct estimate
    w0 = urrr @ vrrr.T

    return w0, urrr, vrrr
    