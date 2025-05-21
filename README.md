# DAIS

A data-adaptive change-point detection method is proposed in Anastasiou & Loizidou (2024) (https://arxiv.org/abs/2404.19344), called DAIS (Data Adaptive ISolation). We work under the model

$X_t = f_t + \sigma \epsilon_t, \quad t = 1, 2, \ldots, T$

where $X_t$ are the observed data and $f_t$ is a one-dimensional deterministic signal with structural changes at certain points. The algorithms provided can be applied to piecewise-constant (changes in the mean) and piecewise-linear (changes in the slope) signals $f_t$. 

### Brief description of the algorithm
The algorithm first identifies the point where the "largest difference" occurs in the interval $[s,e]$, and uses left- and right-expanding intervals around that point. This is done in a sequential way, expanding once either only to the left or only to the right at each step. Each interval is checked by comparing the value of a contrast function with a threshold. This way each change-point is detected in an interval where it is isolated. As soon as a change-point is detected, DAIS restarts on two disjoint intervals, one ending at the start-point of the interval where the detection occurred and one starting from the end-point of that interval.

### Files available
The folder 'R code' contains all relevant code. The file 'Algorithms.R' provides the function DAIS, along with the relevant requred functions. This file is needed for the other R files to run. The file 'Simulations.R' provides code to reproduce the simulations as provided in Section 5 of the main paper and Section 1 of the Supplementary Material of the paper. The folder 'Real data' contains the code used to process the real data included in Section 7 of the main paper. 'Real_data_changes_in_mean.R' and 'Real_data_changes_in_slope.R' concern Sections 7.1 and 7.2, respectively.

