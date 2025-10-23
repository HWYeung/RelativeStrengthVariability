def RSV_exact(W, weighted=0, direction=1, window_size=None, normalisation=0):
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
    NodeRelStrength : ndarray, shape (R, N)
    StrengthVariability : ndarray, shape (N,)
    windowed_StrVar : ndarray, shape (N, ?)
    """

    # Ensure 3D
    if W.ndim == 2:
        W = W[:, :, np.newaxis]

    R, _, N = W.shape
    if window_size is None:
        window_size = R

    # Normalise if needed
    if normalisation:
        W = W / np.max(np.abs(W))

    # Compute relative strength
    W2 = np.sum(W, axis=1, keepdims=True) - W          # sum over rows minus self
    W2 = W2 / np.transpose(W2, (1, 0, 2))             # normalize by column sum
    mask = (W > 0).astype(float)
    Average = np.sum(mask * W2, axis=direction-1) / np.sum(mask, axis=direction-1)

    if direction == 1:
        NodeRelStrength = 1 / np.squeeze(Average)
    else:
        NodeRelStrength = np.squeeze(Average)

    if NodeRelStrength.ndim == 1:
        NodeRelStrength = NodeRelStrength[:, np.newaxis]

    # Variability
    if weighted == 0:
        StrengthVariability = np.std(NodeRelStrength, axis=0)
    else:
        Strength = np.squeeze(np.sum(W, axis=1))
        Weightings = Strength / np.sum(Strength, axis=0)
        WeightedMean = np.sum(NodeRelStrength * Weightings, axis=0)
        SquaredNodeRel = (NodeRelStrength - WeightedMean)**2
        StrengthVariability = np.sqrt(np.sum(Weightings * SquaredNodeRel, axis=0))

    # Early return if full window
    if window_size == R:
        windowed_StrVar = None
        return NodeRelStrength, StrengthVariability, windowed_StrVar

    # Sliding window variability
    MeanNodeStrength = np.mean(np.sum(W, axis=1), axis=1)
    NodeStrengthOrder = np.argsort(MeanNodeStrength)
    NodeRelStrength_sorted = NodeRelStrength[NodeStrengthOrder, :]
    Window_size = int(np.floor(R * window_size))
    windowed_StrVar = np.zeros((N, R - Window_size + 1))
    for i in range(R - Window_size + 1):
        windowed_StrVar[:, i] = np.std(NodeRelStrength_sorted[i:i+Window_size, :], axis=0)
    h_RSV = np.nanmean(windowed_StrVar)
    RSV = StrengthVariability

    return NodeRelStrength, RSV, windowed_StrVar, h_RSV