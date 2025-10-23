import numpy as np

def rsv_exact(W, weighted=0, direction=1, window_size=None, normalisation=0):
    """
    Compute Node Relative Strength and its variability, matching MATLAB version.

    Parameters
    ----------
    W : ndarray
        Structural connectome, shape (R, R, N) or (R, R) for single network
    weighted : int
        1 for weighted variability, 0 for unweighted
    direction : int
        1 for arithmetic mean (inverted later), 2 for harmonic mean
    window_size : int
        Size of sliding window for variability
    normalisation : int
        If 1, normalize W by max absolute value

    Returns
    -------
    node_rs : ndarray, shape (R, N)
    rsv : ndarray, shape (N,)
    windowed_str_var : ndarray, shape (N, ?)
    h_rsv : float
    """

    # Ensure 3D array
    if W.ndim == 2:
        W = W[:, :, np.newaxis]

    n_nodes, _, n_networks = W.shape
    if window_size is None:
        window_size = n_nodes

    # Normalise if needed
    if normalisation:
        W = W / np.max(np.abs(W))

    # Compute relative strength
    rel_strength_matrix = np.sum(W, axis=1, keepdims=True) - W
    rel_strength_matrix = rel_strength_matrix / np.transpose(rel_strength_matrix, (1, 0, 2))
    
    mask = (W > 0).astype(float)
    average_strength = np.sum(mask * rel_strength_matrix, axis=direction-1) / np.sum(mask, axis=direction-1)

    if direction == 1:
        node_rs = 1 / np.squeeze(average_strength)
    else:
        node_rs = np.squeeze(average_strength)

    if node_rs.ndim == 1:
        node_rs = node_rs[:, np.newaxis]

    # Compute variability
    if weighted == 0:
        rsv = np.std(node_rs, axis=0)
    else:
        total_strength = np.squeeze(np.sum(W, axis=1))
        weightings = total_strength / np.sum(total_strength, axis=0)
        weighted_mean = np.sum(node_rs * weightings, axis=0)
        squared_diff = (node_rs - weighted_mean) ** 2
        rsv = np.sqrt(np.sum(weightings * squared_diff, axis=0))

    # Early return if full window
    if window_size == n_nodes:
        windowed_str_var = None
        h_rsv = np.nan
        return node_rs, rsv, windowed_str_var, h_rsv

    # Sliding window variability
    mean_node_strength = np.mean(np.sum(W, axis=1), axis=1)
    node_strength_order = np.argsort(mean_node_strength)
    node_rs_sorted = node_rs[node_strength_order, :]

    window_size = int(np.floor(n_nodes * window_size))
    windowed_str_var = np.zeros((n_networks, n_nodes - window_size + 1))
    for i in range(n_nodes - window_size + 1):
        windowed_str_var[:, i] = np.std(node_rs_sorted[i:i+window_size, :], axis=0)

    h_rsv = np.nanmean(windowed_str_var)

    return node_rs, rsv, windowed_str_var, h_rsv
