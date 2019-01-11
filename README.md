# Vapor Transport Model 
## Kukla et al., 2019, JGR Atmospheres

### Model use / inquiries

**If you use this model for published analysis** please cite the paper which first describes the model:

*Kukla et al., (2019). The sensitivity of terrestrial d18O gradients to hydroclimate evolution. JGR Atmospheres.*

**If you use this model for educational purposes** and would like scripts for pre-loaded demonstrations or interactive tutorials in R please contact tykukla@stanford.edu. Also, check out this [animation demonstrating model output](https://www.youtube.com/watch?v=_7u0DMuS_DU). 

**For any questions or assistance with running the model** please feel free to contact me (tykukla@stanford.edu).


### Running the model

Make sure all files are located in the same working directory. The only file the user needs to adjust is _**ModelRUN.R**_. Other scripts include constants and functions required to solve the model`s set of equations. Info on these scripts is as follows:

#### build_topo_fxn.R

Function to build an idealized topographic domain with a gaussian or witch of agnesi mountain.

#### constants.R

Set of physical constants required in the model calculations. 

#### hydro_fxns.R

Functions that solve water cycle relationships to initialize the hydroclimate solution.

#### isofrac_fxns.R

Functions to handle equilibrium isotope fractionation (kinetic is handled elsewhere because it is dependent on the hydroclimate solution). 

#### ODE_fxns.R

The ordinary differential equation functions for the hydroclimate and isotope solutions.

#### ModelSolve_fxn.R

The main function that reads in all other scripts and returns the model solution. 

### Required packages to run the model

The following packages are required to run the model and generate figures of output:

*Calculations (absolutely required)*
1. pracma
2. rootSolve

*Figures (only required for plotting code that is already written)*

1. ggplot2
2. ggpubr
3. ggthemes
4. scales
          
**NOTE: the ggplot2 package requires R version $\geq$ 3.1
          
*Accessory (not required but used if loaded)*

1. tictoc
