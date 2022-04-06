# Palmer_CTh_Model
Model of thorium sorption kinetics near Palmer Station, Antarctica. This is the model code used in Stukel, M. R., Schofield, O. M. E., Ducklow, H. W., 2022.  Seasonal variability in carbon:234thorium ratios of suspended and sinking particles in coastal Antarctic waters: Field data and modeling synthesis. Deep-Sea Research I 184: 103764.  doi: 10.1016/j.dsr.2022.103764.  Explanations of the equations used can be found in that manuscript.

The files that start with the name "Run" are the files that are used to run the model for each of the 8 model variations described in Stukel et al. 2022.  Please note that since Bayesian simulations like these should be run for many iterations before the number of iterations (variable "numiter") in the code should be increased substantially to get reliable results.  Alternately, what we do (rather than using a large value for numiter and completing everything in a single run) is to restart the code using the corresponding files that begin with the name "Continue".  These "Continue" files will restart from the saved output of the "Run" files (or from a previous call of the "Continue" file.
