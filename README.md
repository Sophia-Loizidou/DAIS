# DAIS
Data adaptive change-point detection algorithm

A data-adaptive change-point detection method is proposed, called DAIS (Data Adaptive ISolation). We work under the model

$X_t = f_t + \sigma \epsilon_t, \quad t = 1, 2, \ldots, T$

where $X_t$ are the observed data and $f_t$ is a one-dimensional deterministic signal with structural changes at certain points. The algorithms provided can be applied to piecewise-constant and piecewise-linear signals $f_t$, meaning that the changes are in the mean or the slope. 

The idea behind the data-adaptive behaviour of the algorithm arises from the fact that, for the unobserved true signal $f_t$ in the case that it is piecewise-constant, taking pairwise differences between consecutive time steps will reveal a constant signal at 0, with "spikes"" where the change-points occur. This means that the largest "spike"", in absolute value, is the location of the change-point with the largest jump in the sequence. Similarly, in the case that $f_t$ is piecewise-linear, differencing the signal twice will again reveal a signal with "spikes" at the locations of the change-points. We try to use this fact in the observed signal, $X_t$. In the algorithm, we identify this "spike" test around it for possible change-points. We check in expanding intervals in such a way that we achieve isolation of the change-points and detection in such an interval. Each interval is checked for possible change-points by comparing the value of a contrast function with a threshold.