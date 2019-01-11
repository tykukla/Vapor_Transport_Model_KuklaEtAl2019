#--------------------------------------------------------------------#
#              Vapor transport model -- main solver                  #
#             ---------------------------------------                #
# T. Kukla (Stanford Univ. 2018)                                     #
#--------------------------------------------------------------------#
library(tictoc)
library(imputeTS)
library(dplyr)
# This script takes care of the necessary calculations to prepare data for the 
# hydroclimate and isotope solutions from the vapor transport ODEs

setwd("/Users/tylerkukla/Desktop")

source('ODE_fxns.R')
source('isofrac_fxns.R')
source('hydro_fxns.R')
source('constants.R')

VTSolve <- function(){
  tic()
  
  ## 1: COLLECT VECTORS -------------------------------------
  if(length(MAT_Kelv) == n){MAT_K <- MAT_Kelv} else{MAT_K <- rep(MAT_Kelv, n)}
  MAT_C <- MAT_K - 273.15
  if(length(rh_0) == n){rh_wb <- rh_0} else{rh_wb <- rep(rh_0, n)}
  if(length(myOMEGA) == n){omega <- myOMEGA} else{omega <- rep(myOMEGA, n)}
  if(length(PET) == n){Eo <- PET} else{Eo <- rep(PET,n)}
  if(length(u_vel) == n){u <- u_vel} else{u <- rep(u_vel,n)}
  if(length(trnsp) == n){t <- trnsp} else{t <- rep(trnsp, n)}
  DI <- myDI[1]
  
  #... general domain info
  dx <- L/n         # meters distance per 1-D grid cell
  x <- seq(dx,L, by=dx)   # distance vector over space
  # calculate slope
  zgrad <- diff(z) / (L/n) ; zgrad[n] <- 0    # calculate gradient to get upslope trajectory of parcel
  for(i in 1:max(n)){
    if(zgrad[i] < 0){zgrad[i] <- 0}
  }
  # modify peclet for Rayleigh-style distillation over orography
  Pe <- zgrad*1e8 + Peclet    # peclet number is high when topography matters, otherwise = Pe climate 
  Win_RBC <- 0      # right boundary condition value of precipitable water (not necessary with Neumann RBC)
  
  
  ## 2: THE ATMOSPHERIC DOMAIN --------------------------------
  ## calculated domain values
  h.atm <- seq(1,13e3,1)     # vector for height of atmosphere (m) (in m above LCL)
  Td <- Tdew_K(MAT_K, rh_wb)   # dewpoint temperature
  SVD_ref <- CC(Td[1]) / (Rm*Td[1])     # saturation vapor density calculated at the lifting condensation level (LCL) (see Smith and Barstad 2004 for eq'n)
  qv <- (B * (CC(MAT_K)/100) / (p.surf - (CC(MAT_K)/100))) / 1e3  # saturated vapor mixing ratio (g/kg) assuming rh=1 (CC divided by 100 because equation takes pressure in hPa, not Pa. whole thing divided by 1e3 to get in kg/kg, not g/kg)
  GAM_m <- GAM_m_fun(qv, Td)         # calculate moist adiabat
  gam_env <- GAM_m / 1               # environmental lapse rate (degK/m)
  Hw <- Hw_fun(Td, gam_env)             # calculate vapor scale height (in m above LCL)
  LCL <- (Td - MAT_K) / GAM_d     # lifting condensation level in (m)
  
  ## set input values
  # calculate integral of Smith model (vert.velocity*exp(-h.atm/Hw))
  f_in <- function(x){(u*zgrad) * exp(-x/Hw)}      # function for the integral
  int_term <- quadv(f_in,min(h.atm),max(h.atm),tol=.Machine$double.eps^(1/2))$Q
  # calculate vapor into domain
  Win_LBC <- SVD_ref * Hw[1]  # this will be overwritten if pw_initial is defined
  if(is.numeric(pw_initial)){Win_LBC <- pw_initial}
  # Wmax - used shortly
  Wmax <- Win_LBC / rh_wb[1] 
  
  # water balance modifiers
  kp <- (Eo[1]/DI) / (Win_LBC**2)  # for P=Mw^2
  
  
  ##### run the model #####
  outwb <- waterbalance.ODE(Pe = Pe, u = u, L = L, int_term = int_term,     
                            n = n, kp = kp, GAM_m = GAM_m, gam_env = gam_env,
                            Win_LBC = Win_LBC, Win_RBC = Win_RBC, omega  = omega, 
                            Eo = Eo, zgrad = zgrad, Hw=Hw)
  
  P.tot <- outwb$P
  ET.tot <- outwb$ET
  w.tot <- outwb$W
  SW.tot <- outwb$soil.w
  
  ## 3: ISOTOPE DOMAIN --------------------------------------
  # set isotope standard and RBC 
  isot.standard <- "VSMOW"  # either "VSMOW" or "VPDB"
  if(isot.standard == "VSMOW"){r.std <- 2005.2 / 1e6}  # vsmow ratio of 18O/16O 
  if(isot.standard == "VPDB"){r.std <- 2067.2 / 1e6}    # isotope standard ratio (VPDB)
  if(isot.standard != "VSMOW" & isot.standard != "VPDB"){print("ERROR: Incorrect isotope standard. Select 'VSMOW' or 'VPDB'")}
  
  d.in_RBC <- -1e3  # right boundary condition. This value must be -1000 for zero gradient out (Neumann)
  
  
  #... isotope fractionation calculations
  # kinetic
  n_MB <- 0.8     # n value from framework of Mathieu and Bariac 1996
  k <- 1-(n_MB*(1-0.9687))    # 18O kinetic fractionation factor (dry-ish soil)
  
  # calculate RH, Tdew and therefore fractionations based on output w
  ## calculate the temperature profile for precipitation 
  ## assume background precip is the temp of the dewpoint at the LCL
  ## assume orographic precip is at LCL or surface topo, whichever is lower
  rh <- w.tot / Wmax
  # rh <- (w.tot/Hw) / (CC(MAT_K)/(Rm*MAT_K))    # relative hum as a function of vapor remaining
  LCL <- (Tdew_K(MAT_K, rh) - MAT_K) / GAM_d
  T_et <- T_prec <- vector()
  for(i in 1:length(z[1:n])){
    if(z[i] <= LCL[i]){T_et[i] <- MAT_K[i] + GAM_d*z[i]
    T_prec[i] <- Tdew_K(MAT_K[i], rh[i])}
    if(z[i] > LCL[i]){T_et[i] <- T_prec[i] <- Tdew_K(MAT_K[i], rh[i]) + GAM_m[i]*(z[i]-LCL[i])}
  }
  # use these temps to calculate fractionation factors 
  # equilibrium
  a_P <- a_E <- vector()
  
  a_P[c(which(T_prec < 253.15))] <- sapply(T_prec[which(T_prec<253.15)], a18_P_iv)
  a_P[c(which(T_prec >= 253.15 & T_prec <= 273.15))] <- sapply(T_prec[which(T_prec >= 253.15 & T_prec <= 273.15)], a18_P_mix)
  a_P[c(which(T_prec > 273.15))] <- sapply(T_prec[which(T_prec > 273.15)], a18_P_lv)
  a_P <- unlist(a_P)
  
  a_E <- sapply(T_et, a18_e_vl) 
  a_E <- unlist(a_E)
  
  ##### calculate isot composition of input vapor #####
  ## source waters
  d.in_LBC <- -12   # [per mille] isotopic composition of vapor at first grid cell
  
  
  ## run the isotope balance 
  outis <- isotopes.ODE(Pe=Pe, u=u, L=L, n=n, kp=kp, rh = rh,
                        d.in_LBC=d.in_LBC, d.in_RBC=d.in_RBC, t = t,
                        r.std=r.std, k=k, a_P=a_P, a_E=a_E, P.tot = P.tot, 
                        ET.tot=ET.tot, SW.tot=SW.tot, w.tot=w.tot)
  
  
  ## 4: RETURN RESULTS --------------------------------------
  # First make a clim dataframe
  outclim <- as_tibble(cbind(T_et, rh, PET/outwb$P, LCL)) 
  colnames(outclim) <- c('Temp', 'rh', 'DI', 'LCL')
  
  # bring it all together
  outdf <- as_tibble(cbind(outwb, outis, outclim))
  outdf$x_km <- x / 1e3
  
  #... track total solve time
  toc()
  
  #... return output 
  return(outdf)

}

