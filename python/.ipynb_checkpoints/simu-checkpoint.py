import numpy as np

def get_or(d, key, default):
    val = d.get(key, default)
    d[key] = val
    return val, d

def simu_RRR(ops=None):
    if ops is None:
        ops = {}

    # Load or default parameters
    T, ops = get_or(ops, 'T', 100)
    nx, ops = get_or(ops, 'nx', 10)
    ny, ops = get_or(ops, 'ny', 5)
    rnk, ops = get_or(ops, 'rnk', 1)
    signse, ops = get_or(ops, 'signse', 0.1)
    magnitude, ops = get_or(ops, 'magnitude', 1)
    thetas, ops = get_or(ops, 'thetas', [])
    Sigma, ops = get_or(ops, 'Sigma', [])
    U, ops = get_or(ops, 'U', [])
    V, ops = get_or(ops, 'V', [])

    # Generate input X
    X = np.random.randn(T, nx)
    X -= X.mean(axis=0)  # center each column

    # Generate communication matrix
    if U == [] or V == []:
        U = np.random.randn(nx, rnk)
        V = np.random.randn(ny, rnk)
        if len(thetas) > 0:
            U_svd, _, V_svd = np.linalg.svd(U @ V.T, full_matrices=False)
            U = U_svd[:, :rnk]
            V = V_svd[:rnk, :].T / np.sqrt(thetas)
        V = V * magnitude

    B = U @ V.T

    # Generate Y
    if len(Sigma) == 0:
        E = signse * np.random.randn(T, ny)
        # E -= E.mean(axis=0)
        Y = X @ B + E
    else:
        E = np.random.multivariate_normal(np.zeros(ny), Sigma, size=T)
        Y = X @ B + E

    return X, Y, U, V, ops
