#--------------------------------------------------------------------#
#              Vapor transport model -- hydro fxns                   #
#             ---------------------------------------                #
# T. Kukla (Stanford Univ. 2018)                                     #
#--------------------------------------------------------------------#

# This script contains functions required to solve the hydroclimate fluxes over the domain

## 1: Clausius Clapeyron relationship (input temp in K, output pressure in Pa)
CC <- function(temp_K){(A*exp(beta*(temp_K-273.15)))*100}  

## 2: Dew point temperature function (from Lawrence, 2005)
Tdew_K <- function(temp_K, rh){((temp_K-273.15) - ((1-rh)/.05)) + 273.15} 

## 3: Precipitable water estimate (Smith 1966 Journal of Applied Meteorology)
#  [not used in the model, but available in case the scale height approach 
#   returns unreasonable values, which can happen in the tropics]
#   --see Smith, 1966 for the definition of lambda--
S66 <- function(lambda=0.5, Tdw){
  temp_s <- exp(0.1183 - log(lambda + 1) + 0.0707*(Td-273.15))
  exp(temp_s)
}

## 4: Vapor scale height function (see Smith and Barstad, 2004 appendix)
Hw_fun <- function(Tdew_K, gam_env){-Rm*(Tdew_K**2) / (Lh*gam_env)}

## 5: Moist adiabat function
GAM_m_fun <- function(mix_ratio, temp_K){
  -(g * (1+(Lh*mix_ratio / (Rd*temp_K))) / (cpd+(((Lh**2) * mix_ratio) / (Rm*(temp_K**2))))) }
