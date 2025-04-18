# Aquabreeding2 (test ver)

Aquabreeding2 software is a python tool for simulating aquaculture breeding.

Aquabreeding2 generates a founder population using coalescent simulation implemented in [msprime](https://tskit.dev/msprime/docs/stable/intro.html).  A progeny population is produced by crossing the founder individuals, and the progenies with larger breeding/phenotypic values are selected as parents of the next generation.  Phenotypic values, breeding values, variance components, and inbreeding coefficient are calculated in each generation.

This is a beta version.  Please do not use it for any publications.


## Last update
- 04/18/2025


## Note
- Aquabreeding2 supports python3.13.  


## Requirements
- macOS, Linux
- python3 (>=3.11)  
  (e.g., brew install python3)
- python libraries (automatically installed)
    - numpy
    - scipy
    - [msprime](https://tskit.dev/msprime/docs/stable/intro.html)  
    - pandas (for jupyter notebook)  
    - matplotlib (for jupyter notebook)  
    - seaborn (for jupyter notebook)  
- cmake  
  (e.g., brew install cmake)
- C++ library
    - eigen  
      (e.g., brew install eigen)


## Installation
`git clone --recursive https://github.com/showhey0119/aquabreeding2`  
`cd ./aquabreeding2`  
`git submodule update --init --recursive`  
`pip3 install .`  

If the installation fails, try
`pip3 install -e .`  



## Quick startup
see [aquabreeding\_tutorial.ipynb](https://github.com/showhey0119/aquabreeding/blob/master/aquabreeding_tutorial.ipynb)


## History
- Version 0.0.1 is released (multiple breeding population)
