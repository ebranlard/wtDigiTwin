[![Build Status](https://travis-ci.org/ebranlard/yams.svg?branch=master)](https://travis-ci.org/ebranlard/yams)

# wtDigiTwin
A digital twin model for wind turbine.
Contains a set of structrual dynamics tools, beam theory, FEM and more.



Applications were presented in:

- Branlard,E, Giardina, D., Brown, C. S. , *Augmented Kalman filter with a reduced mechanical model to estimate tower loads on a land-based wind turbine: a step towards digital-twin simulations*, 2020 [link](https://doi.org/10.5194/wes-5-1155-2020)

- Branlard,E, Jonkman, J., Dana, S., Doubrawa, P., *A digital twin based on OpenFAST linearizations for real-time load and fatigue estimation of land-based turbines*, 2020 [link](https://iopscience.iop.org/article/10.1088/1742-6596/1618/2/022030)

The structural model rely wither on OpenFAST linearizations or on the YAMS model described here:

 - Branlard,E, *Flexible multibody dynamics using joint coordinates and the Rayleigh‐Ritz approximation: The general framework behind and beyond Flex* , 2019, [link](https://onlinelibrary.wiley.com/doi/abs/10.1002/we.2327). A pre-print of this article is available in the `_doc` folder of this repository.




## QuickStart
Download, install dependencies and package:
```bash
git clone https://github.com/ebranlard/wtDigiTwin
cd wtDigiTwin
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


Some of the packages have dependency with the [weio](http://github.com/ebranlard/weio/) library to read and write files.





