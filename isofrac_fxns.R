#--------------------------------------------------------------------#
#              Vapor transport model -- iso fractionation            #
#             -------------------------------------------            #
# T. Kukla (Stanford Univ. 2018)                                     #
#--------------------------------------------------------------------#

# This script contains necessary functions for isotope fractionation 

## 1: equilibrium of 18/16 O (with condensation fractionation mixing model (see Rowley 2007))
a18_P_lv <- function(temp){  # function for alpha of precipitation from vapor to liquid phase
  1 + (-7.685 + (6.7123*1e3)/(temp) - (1.6664*1e6)/(temp)^2 + (0.35041*1e9)/(temp)^3)/1e3 
}
a18_P_iv <- function(temp){   # function for alpha of precipitation from vapor to ice phase
  1 + ((11839/temp) - 28.224)/1e3
}
a18_P_mix <- function(temp){  # function for mixing line between liq and ice condensate functions (see Rowley 2007)
  -((a18_P_iv(253.15) - a18_P_lv(273.15)) / 20) * temp + 1.104548
}
a18_e_vl <- function(temp){   # function for alpha of evaporation from liquid to vapor phase
  1 - (-7.685 + (6.7123*1e3)/(temp) - (1.6664*1e6)/(temp)^2 + (0.35041*1e9)/(temp)^3)/1e3 
}




