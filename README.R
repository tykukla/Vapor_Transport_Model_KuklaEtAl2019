#--------------------------------------------------------------------#
#              Vapor transport model -- README                       #
#             ---------------------------------------                #
# T. Kukla (Stanford Univ. 2018)                                     #
#--------------------------------------------------------------------#

# ----- MODEL USE / INQUIRIES ----- # 
# If you use this model for published analysis please cite the 
# paper which first describes the model:
# Kukla et al. "The sensitivity of terrestrial d18O gradients to
# hydroclimate evolution"
# ---- #
# For any questions or assistance with running the model, 
# please contact: 
# Tyler Kukla - tykukla@stanford.edu
# --------------------------------- #

# The folder with this document contains files required to run the vapor 
# transport model described in Kukla et al., "The sensitivity of 
# terrestrial d18O gradients to hydroclimate evolution"

# ------------

# TO RUN THE MODEL:
# --- The only file the user needs to adjust is:
#     ModelRUN.R

# --- other scripts include constants and functions required 
#     to solve the model's set of equations, info on those 
#     scripts is as follows:
#     --build_topo_fxn.R -- function to build an idealized 
#                           topographic domain with a gaussian or 
#                           witch of agnesi mountain 
#     --constants.R -- set of physical constants required
#     --hydro_fxns.R -- functions that solve water cycle relationships
#                       to initialize the hydroclimate solution
#     --isofrac_fxns.R -- functions to handle equilibrium fractionation
#     --ODE_fxns.R -- the ordinary differential equation functions for the
#                     hydroclimate and isotope solutions
#     --ModelSolve_fxn.R -- the main function that reads in all other 
#                           scripts and returns the model solution

# REQUIRED PACKAGES IN R:
# --- The following packages are required to run the model
## CALCULATIONS - absolutely required
# 1) pracma
# 2) rootSolve
## PLOTS - only required for scripted plots to show up in viewer
# 1) ggplot2
# 2) ggpubr
# 3) ggthemes
# 4) scales
# ** NOTE - ggplot2 requires R version >/= 3.1 
## ACCESSORY - not required but used if loaded
# 1) tictoc

