[![Build Status](https://travis-ci.org/ebranlard/yams.svg?branch=master)](https://travis-ci.org/ebranlard/yams)

# YAMS
YAMS (Yet Another Multibody Solver) - Structural dynamics tools, beam theory, FEM and more.
Both python and matlab tools are available.

This repository started with some example code from the following article:
Branlard,E, *Flexible multibody dynamics using joint coordinates and the Rayleigh‐Ritz approximation: The general framework behind and beyond Flex* , 2019, [link](https://onlinelibrary.wiley.com/doi/abs/10.1002/we.2327). A pre-print of this article is available in the `_doc` folder of this repository.

It has since been extended to include more structural dynamics tools.



## QuickStart
Download, install dependencies and package:
```bash
git clone https://github.com/ebranlard/yams
cd yams
python -m pip install --user -r requirements.txt  
python -m pip install -e .      # install
python -m unittest discover -v  # run test
```

## Packages
The repository contains a set of small packages:

- yams: multibody analyses
- beams: analytical results for beams
- fast: tools to handle OpenFAST models, in particular a linear OpenFAST model
- fem: Finite Element Method tools (beams)
- kalman: kalman filter
- system: tools for dydnamic systems (e.g. LTI, state space) and mechanical systems (M,C,K matrices), eigenvalue analysis, time integration
- ode: tools for time integration of ODE
- ws\_estimator: wind speed estimator for wind energy based on tabulated Cp Ct


Some of the package may have dependency with the [weio](http://github.com/ebranlard/weio/) library to read and write files.





