#--------------------------------------------------------------------#
#              Vapor transport model -- constants                    #
#             ---------------------------------------                #
# T. Kukla (Stanford Univ. 2018)                                     #
#--------------------------------------------------------------------#

# script that brings in the physical constants required for the hydroclimate solution

## physical constants
A <- 6.11            # [hPa] clausius clapeyron constant
beta <- 0.067        # [degC-1] clausius clapeyron constant 
B <- 621.9907        # [hPa] mixing ratio constant 
p.surf <- 1013.25    # [hPa] surface pressure 
g <- 9.8076          # [m s-2] gravitational accel constant 
Lh <- 2.5e6          # [J kg-1] latent heat of vaporization for water
Rd <- 287.05         # [J kg-1 degK-1] specific gas constant for dry air
Rm <- 461.495        # [J kg-1 degK-1] specific gas constant for moist air
cpd <- 1003.5        # [J kg-1 degK-1] specific heat dry air
GAM_d <- -(g/cpd)    # [degK m-1] dry adiabatic lapse rate 
