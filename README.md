# AnnealingCode
Code to use simulated annealing on random truncated non-commutative geometries. This code was used to generate the two d results in https://arxiv.org/abs/1909.08054

This code relies on three external libaries. 

The highfive, header only implementation of hdf5 https://github.com/BlueBrain/HighFive ( used version 2.0)

The Eigen::MatrixXcd class is implemented in the Eigen Matrix library https://eigen.tuxfamily.org/ and my code uses version 3.24 of it.

The code also requires a random number generator which was used from here class https://www.agner.org/random/
