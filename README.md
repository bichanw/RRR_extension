# RRR_extension
Matlab and Python code providing extensions from reduced-rank regression (RRR, from [Semedo et al. 2019](https://doi.org/10.1016/j.neuron.2019.01.026)). We offer the following functionality:
1. adding ridge regularization
2. adding non-isotropic noise support
3. adding quantification of alignment between activity space and communication subspace, for input and output population respectively

We also add codes to generate random examples and reproduce the figures from the manuscript (in prep).


## Matlab 
### Usage
- Navigate your matlab environment to repo directory, run `startup.m`
- Alternatively, add current directory to your matlab path to run functions freely

### Folder structure
- main functions
    - `svd_RRR.m`: using svd method to perform RRR, with optional ridge regularization 
    - `svd_RRR_noniso.m`: RRR with non-isotropic noise
    - `alignment_input.m`: function that caculates the alignment index between input population activity and communication subspace
    - `alignment_output.m`: function that caculates the alignment index between output population activity and communication subspace
    - `main_figures.m`: contains code that utilizes the above functions and reproduce the figures in paper
- `other_funcs`: folder containing supportive functions for simulation, plotting etc.


## Python
### Usage
- Our code only use very basic functions from `numpy`, `scipy`, and `matplotlib`. So for any Python environment setup newer than 3.7, it should function normally.

### Folder structure
- `figures.ipynb`: the main jupyter notebook that runs the simulation and produces the figures
- `fitting.py`: main script with functions that do RRR inference
    - `svd_RRR`: using svd method to perform RRR, with optional ridge regularization 
    - `svd_RRR_noniso`: RRR with non-isotropic noise
- `alignment.py`: main script with functions that calculate alignment indexes
    - `alignment_input`: function that caculates the alignment index between input population activity and communication subspace
    - `alignment_output`: function that caculates the alignment index between output population activity and communication subspace
- `simu.py`: main script that simulates datasets for figures