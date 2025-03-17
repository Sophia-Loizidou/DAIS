# DAIS

A data-adaptive change-point detection method is proposed in Anastasiou & Loizidou (2024) (https://arxiv.org/abs/2404.19344), called DAIS (Data Adaptive ISolation). We work under the model

$X_t = f_t + \sigma \epsilon_t, \quad t = 1, 2, \ldots, T$

where $X_t$ are the observed data and $f_t$ is a one-dimensional deterministic signal with structural changes at certain points. The algorithms provided can be applied to piecewise-constant (changes in the mean) and piecewise-linear (changes in the slope) signals $f_t$. 

### Brief description of the algorithm
The algorithm first identifies the point where the "largest difference" occurs in the interval $[s,e]$, and uses left- and right-expanding intervals around that point. This is done in a sequential way, expanding once either only to the left or only to the right at each step. Each interval is checked by comparing the value of a contrast function with a threshold. This way each change-point is detected in an interval where it is isolated. As soon as a change-point is detected, DAIS restarts on two disjoint intervals, one ending at the start-point of the interval where the detection occurred and one starting from the end-point of that interval.

### Idea behind data-adaptiveness
For the unobserved true signal $f_t$ in the case that it is piecewise-constant, taking pairwise differences between consecutive time steps will reveal a constant signal at 0, with "spikes" where the change-points occur. This means that the largest "spike", which satisfies $\textrm{argmax}\{\lvert X_{t+1} - X_{t} \rvert\}$ is the location of the change-point with the largest jump in the sequence. Similarly, in the case that $f_t$ is piecewise-linear, differencing the signal twice (i.e. using $\textrm{argmax}\{\lvert X_{t+2} - 2X_{t+1} + X_{t} \rvert\}$) will again reveal a signal with "spikes" at the locations of the change-points. We try to use this fact in the observed signal, $X_t$.

