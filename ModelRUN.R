#--------------------------------------------------------------------#
#              Vapor transport model -- model run (DIY)              #
#             -----------------------------------------              #
# T. Kukla (Stanford Univ. 2018)                                     #
#--------------------------------------------------------------------#
library(rootSolve)
library(pracma)

rm(list=ls())  # clear global environment

# This is the main script to run the model with user-prescribed inputs
# Decisions required to initialize the simulation should be defined here
# under "CONSTRUCT DOMAIN" and "CLIMATE CONDITIONS"
# BEFORE RUNNING make sure you change the working directory 
# AND load required packages outlined in the README file
# ---- Running the entire script through will generate example model output
#      and plot the output in 4 panels in the plots viewer in Rstudio

#... identify working directory (must have other model scripts)
setwd('/Users/tylerkukla/Desktop/ForDistribution')

#... bring in the function that solves the model
source('ModelSolve_fxn.R')

#... FIRST a potentially useful function:
# if you need to switch from the more common [mm yr-1] units to required
# [kg m-2 s-1] units for potential evapotranspiration, this is the 
# conversion function 
flux_conv <- function(flux){
  secs_per_min <- 60
  min_per_hr <- 60
  hr_per_day <- 24
  day_per_yr <- 365
  
  secs_per_yr <- secs_per_min * min_per_hr * hr_per_day * day_per_yr
  
  flux * (1/secs_per_yr)
}


## 1: CONSTRUCT DOMAIN ---------------------------------------------- 
n <- 500            # number of nodes in the domain
L <- 1e6            # [m] length of domain
kL <- L/1e3         # [km] length of domain
# source('build_topo_fxn.R')  # provides a function to build your domain--set equal to 'z'
z <- rep(0, n)      # flat domain


## 2: CLIMATE CONDITIONS --------------------------------------------
MAT_Kelv <- 290       # [K] mean annual temperature at 2m above sea level
PET <- flux_conv(1e3) # [kg m-2 s-1] Potential ET (flux_conv converts from mm/yr to kg/m2s)
Peclet <- 100         # peclet number - ratio of advection to diffusion
pw_initial <- 20      # [kg m-2] initial vapor content (column integrated)
rh_0 <- 0.8           # relative humidity of moisture source
u_vel <- 10           # [m s-1] advective velocity
trnsp <- 0.64         # transpired fraction of ET (see Good et al., 2015)
myDI <- 0.5           # starting dryness index
myOMEGA <- 2.6        # prescribed recycling efficiency parameter
    

## 3: RUN MODEL -----------------------------------------------------
df <- VTSolve()



# *************************************************************************************************** #
# --------------------------------------------------------------------------------------------------- #
# *************************************************************************************************** #

## 4: FIGURES --------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(scales)

## PRECIPITATION ISOTOPES
d18P <- ggplot(df) + 
  geom_line(aes(x=x_km, y=d18_P), color='#15959F', size=2) +
  labs(x="Distance [km]", y=expression('δ'^18*'O')) +
  ggtitle("Precipitation isotopes") +
  theme_par() + 
  theme(plot.title=element_text(face='bold', hjust=0),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

## PRECIPITABLE WATER 
pw <- ggplot(df) + 
  geom_line(aes(x=x_km, y=W), color='#43617D', size=2) +
  labs(x="Distance [km]", y=expression("Precipitable water [kg"/'m'^3*']')) +
  ggtitle("Precipitable water") +
  theme_par() + 
  theme(plot.title=element_text(face='bold', hjust=0),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

## TEMPERATURE
tmp <- ggplot(df) + 
  geom_line(aes(x=x_km, y=Temp), color='#C7402D', size=2) +
  labs(x="Distance [km]", y="Temperature [K]") +
  ggtitle("Temperature") +
  theme_par() + 
  theme(plot.title=element_text(face='bold', hjust=0),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

## RELATIVE HUMIDITY
rh <- ggplot(df) + 
  geom_line(aes(x=x_km, y=(rh*100)), color='#C7402D', size=2) +
  labs(x="Distance [km]", y="Relative humidity [%]") +
  ggtitle("Relative humidity") +
  theme_par() + 
  theme(plot.title=element_text(face='bold', hjust=0),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

## BUDYKO 
#... first set up limits
maxDI <- 4
mylimsX <- c(0,1,maxDI)     # lines for budyko limits
mylimsY <- c(0,1,1)     # lines for budyko limits
mylims <- as_tibble(cbind(mylimsX, mylimsY))
colnames(mylims) <- c("X", "Y")
#... find max elevation for diverging color (or set to zero)
maxLoc <- which(z==max(z))
if(length(maxLoc) > 1){maxLoc <- 1.001}
#... now plot it up
budy <- ggplot(df) + 
  # budyko
  geom_line(data=mylims, aes(x=X, y=Y)) + 
  # model output
  geom_point(aes(x=DI, y=(ET/P),color=x_km), alpha=0.5, size=2) +
  # plot adjust
  scale_color_gradientn(name="Distance [km]", colors=c("#03717D", "#F2F8F8", "#D93625"),
                        values=rescale(c(1, maxLoc, n))) + 
  labs(x=expression("Dryness index [E"['o']/'P]'), y=expression("ET"/"P")) +
  ggtitle("Budyko") + 
  theme_par() + 
  theme(plot.title=element_text(face='bold', hjust=0), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position=c(0.8,0.4))

#





