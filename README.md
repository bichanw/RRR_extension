# RRR_extension
A code base extended from reduced-rank regression (RRR, cite Semedo et al.). We provide 2 simple functions that add regularization and non-isotropic extensions.

## matlab 
- `svd_RRR`: using svd method to perform RRR, with the 4th optional input `lambda` as the regularization strength
- `svd_RRR_noniso`: RRR with non-isotropic noise; either input the covariance function (`C` as the 4th input) or estimate from data
- `main_figures.m`: code to produce the figures in paper

## python
