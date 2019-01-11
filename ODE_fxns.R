#--------------------------------------------------------------------#
#              Vapor transport model -- ODE functions                #
#             ----------------------------------------               #
# T. Kukla (Stanford Univ. 2018)                                     #
#--------------------------------------------------------------------#

# precipitation and isotope balance functions

# 1: WATER BALANCE
waterbalance.ODE <- function(Pe = Pe, u = u, L = L, int_term = int_term,     
                             n = n, kp = kp, GAM_m = GAM_m, gam_env = gam_env, 
                             Win_LBC = Win_LBC, Win_RBC = 0, omega  = omega, 
                             Eo = Eo, zgrad = zgrad, Hw=Hw){
  
  #domain parameters
  dx <- L/n
  x <- seq(dx,L,by = dx)
  
  # D <- (u*dx)/Pe # m2/sec
  lscale <- 750*1e3     # length scale of transport (in m)
  D <- (u*lscale)/Pe 
  
  C_i_minus_one <- (D/(dx^2))+(u/(2*dx))
  C_center <- -(2*D)/(dx^2)
  C_i_plus_one <- (D/(dx^2))-(u/(2*dx))
  
  # Construct the matrix
  A <- diag(x=C_center,nrow=n,ncol=n)
  A[row(A)==col(A)-1] <- C_i_plus_one[1:(n-1)]
  A[row(A)==col(A)+1] <- C_i_minus_one[2:n]
  head(A)
  # Substitute appropriate boundary conditions into the corners of A
  #... consistent with Dirichlet boundary conditions on the left side of domain
  A[1,1] = 1 
  A[1,2] = 0
  #... and Neumann conditions on the right side 
  A[n,n] = C_center[n] + C_i_plus_one[n]
  
  # forcing vector
  bb = rep(0,n)           
  # masking vector for the reaction terms f
  # used with dirichlet conditions    
  mask = rep(1,n)
  
  #modify bb
  #... Dirichlet boundary conditions
  bb[1] = -Win_LBC
  mask[1] = 0
  #... Neumann boundary conditions
  bb[n] = -Win_RBC*dx/D[n]*C_i_plus_one[n]
  
  # function for smith model for orographic precip with background rainout
  Pfun = function(w){
    (((w / Hw)*(GAM_m / gam_env))/Hw * int_term) + (kp * (w**2))
  }
  
  # function for the Budyko relationship
  ETfun = function(w) {
    Pfun(w) + Eo - ( (Pfun(w))^omega + Eo^omega)^(1/omega)  
  }
  
  # Define function to solve
  fun <- function(w) {
    A %*% w + mask*(- Pfun(w) + ETfun(w) ) + bb
  }
  
  # Call multiroot to solve it
  W <- multiroot(f=fun,start=rep(0.1,n), 
                 maxiter = 20000, rtol = 1e-16, atol = 1e-10, positive = T)  
  
  w.out <- W$root
  
  
  # Recalculate W, P, ET, and Runoff
  P.tot <- Pfun(w.out) # Precipitation (kg/m2/sec)
  ET.tot <- ETfun(w.out)
  soil.w <- (P.tot-ET.tot) # Runoff (kg/m2/sec) 
  #recycling rate
  m <- P.tot/soil.w
  
  #integrate runoff from right to left (assuming advection to the right) 
  cumulative_discharge = rev(cumsum(rev(soil.w*dx)))
  # can also do sum(runoff)
  # print(cumulative_discharge[1])
  # should be what comes in (Win_LBC).  
  
  # Return Values
  result <- as.data.frame(cbind(w.out,P.tot,ET.tot,soil.w) )
  colnames(result) <- c("W","P","ET","soil.w")
  return(result)
}


# 3: ISOTOPES
isotopes.ODE <- function(Pe=Pe, u=u, L=L, n=n, kp=kp, rh = rh, 
                         d.in_LBC=d.in_LBC, d.in_RBC=-1e3, t = t,
                         r.std=r.std, k=k, a_P=a_P, a_E=a_E, P.tot = P.tot, 
                         ET.tot=ET.tot, SW.tot=SW.tot, w.tot=w.tot){
  
  dx <- L/n
  x <- seq(dx,L,by = dx)
  
  # D <- (u*dx)/Pe # m2/sec
  lscale <- 500*1e3     # length scale of transport (in m)
  D <- (u*lscale)/Pe 
  
  # fill node 1 conditions
  Q <- D/(dx**2)
  R <- u/(2*dx)
  
  #... set up the matrix input vectors
  C_plus <- vector() ; C_plus[1] <- Q[1]*(w.tot[2]/w.tot[1]) - R[1]
  C_plus[n]<- R[n] + Q[n]*(2-(w.tot[(n-1)]/w.tot[n]))
  C_center <- vector() ; C_center[1] <- -2*Q[1]
  C_center[n] <- -2*Q[n]
  C_minus <- vector() ; C_minus[1] <- R[1]+Q[1]*(2 - (w.tot[2]/w.tot[1]))
  C_minus[n] <- Q[n]*(w.tot[(n-1)]/w.tot[(n)]) - R[n]
  for(i in 2:(n-1)){ # fill in coefficients
    C_plus[i] <- Q[i]-R[i] + (Q[i]/(2*w.tot[i]))*(w.tot[i+1]-w.tot[i-1])
    C_center[i] <- -2*Q[i]
    C_minus[i] <- Q[i]+R[i] + (Q[i]/(2*w.tot[i]))*(w.tot[i-1]-w.tot[i+1])
  }
  
  # Construct the matrix
  A <- diag(x=C_center,nrow=n,ncol=n)
  A[row(A)==col(A)-1] <- C_plus[1:(n-1)]
  A[row(A)==col(A)+1] <- C_minus[2:(n)]

  # Substitute appropriate boundary conditions into the corners of A
  #... for Dirichlet boundary conditions on the left of the domain
  A[1,1] = 1 
  A[1,2] = 0
  
  #... and for Neumann on the right
  A[n,(n-1)] <- C_plus[n] + C_minus[n]
  
  # make a forcing vector
  bb = rep(0,n)   
  # make a masking vector for the reaction terms f
  # to use with dirichlet conditions
  mask = rep(1,n)
  
  # quick re-naming
  alpha_P <- a_P
  alpha_E <- a_E
  
  # set boundaries 
  rin_RBC <- r.std*((d.in_RBC/1e3) + 1)
  rin_LBC <- r.std*((d.in_LBC/1e3) + 1)
  
  # modify bb
  #... left side is Dirichlet
  bb[1] = -rin_LBC
  mask[1] = 0
  #... right side is Neumann
  bb[n] <- rin_RBC*C_plus[n]*(2*dx/D[n])
  
  # build intermediate vectors to feed functions
  frac <- vector() ; frac <- (ET.tot/P.tot)    # fraction of water evaporated 
  # build ET, P, and solve functions
  # define ET function
  ETfun = function(r){ # from modified from Craig and Gordon 1965
    r_ET <- (1-t) * ((alpha_E * k * (alpha_P*r)) / (1-rh + (1-t)*k*rh)) + t*(alpha_P*r)*(1 / (1 + (1-t)*k*(rh/(1-rh))))
    dx_P <- (alpha_P*r/r.std - 1) * 1000 # delta notation of precip sample
    dx_E_o <- (r_ET/r.std - 1) * 1000     # ratio of ET 
    eps_E <- (dx_E_o - dx_P) / 1000     # epsilon value for fractionation
    capdelt <- -(1-frac) * ((1000 / (1+eps_E)) * ((1-frac)^eps_E) - 1000)
    dx_E <- dx_P + capdelt
    dx_ET <- (t * dx_P) + (1-t) * dx_E
    rx_ET <- (dx_ET/1000 + 1) * r.std
    (ET.tot/w.tot) * (rx_ET - r)
  }
  
  # define precip function
  Pfun <- function(r){
    (P.tot/w.tot)*(r*alpha_P-r)  
  }
  
  # define function to solve
  fun <- function(r) {
    A %*% r + mask*(- Pfun(r) + ETfun(r) ) + bb
  }
  
  # call multiroot to solve it
  R <- multiroot(f=fun,start=rep(rin_LBC,n), useFortran = FALSE,
                 maxiter = 20000, rtol = 1e-16, atol = 1e-20, positive = FALSE)  
  
  R.out <- R$root
  
  # Recalculate P and ET isotopes
  P_outfun <- function(R.out){alpha_P * R.out}
  ET_outfun <- function(R.out){
    r_ET <- (1-t) * ((alpha_E * k * (alpha_P*R.out)) / (1-rh + (1-t)*k*rh)) + t*(alpha_P*R.out)*(1 / (1 + (1-t)*k*(rh/(1-rh))))
    dx_P <- (alpha_P*R.out/r.std - 1) * 1000 # delta notation of precip sample
    dx_E_o <- (r_ET/r.std - 1) * 1000     # ratio of ET 
    eps_E <- (dx_E_o - dx_P) / 1000     # epsilon value for fractionation
    capdelt <- -(1-frac) * ((1000 / (1+eps_E)) * ((1-frac)^eps_E) - 1000)
    dx_E <- dx_P + capdelt
    dx_ET <- (t * dx_P) + (1-t) * dx_E
    (dx_ET/1000 + 1) * r.std
  }
  
  P.isot <- P_outfun(R.out)
  ET.isot <- ET_outfun(R.out) 
  
  # convert to delta notation  
  d18_vap <- ((R.out/r.std) -1) * 1e3
  d18_P <- ((P.isot/r.std) -1) * 1e3
  d18_ET <- ((ET.isot/r.std) -1) * 1e3
  
  # build data frame and label columns
  isot_out <- as.data.frame(cbind(d18_vap,d18_P,d18_ET))
  colnames(isot_out) <- c('d18_vap', 'd18_P', 'd18_ET')
  
  return(isot_out)
}

