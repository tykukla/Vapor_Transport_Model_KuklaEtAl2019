#--------------------------------------------------------------------#
#              Vapor transport model -- topo builder                 #
#             ---------------------------------------                #
# T. Kukla (Stanford Univ. 2018)                                     #
#--------------------------------------------------------------------#
library(imputeTS)

# This script contains the functions to build one-dimensional idealized topographies
# with a mountain range that is shaped like a gaussian fxn or a 'witch of agnesi' ridge

## 1: Topo builder 
#     Arguments:
#               max_elev -- [m] the peak elevation of the mountain
#               peak_loc -- [node] the node along the domain where the mountain's peak resides
#               hhhw -- [nodes] the half-width of the mountain at half its height (half-height-half-width)
#               nodes -- [] number of nodes along the domain
#               gaussian -- logical: TRUE means construct a gaussian mountain, FALSE means witch of agnesi
#               lee_plateau -- logical: TRUE means construct a plateau behind the mountain (FALSE means don't)
#               plateau_elev -- [m] the elevation of the leeward plateau
topo_builder <- function(max_elev, peak_loc, hhhw, nodes, gaussian=TRUE, lee_plateau=FALSE, plateau_elev=NULL){
  # first... let's just make sure that things make sense
  if(peak_loc > nodes){stop('Peak location cannot be greater than the number of nodes in domain')}
  if(max_elev < 0 ){stop('Max elevation cannot be negative')}
  # now... set a switch based on the shape of the mountain
  ifelse(gaussian==TRUE, shape_switch <- 1, shape_switch <- 2)
  switch(shape_switch, 
         {# if we're building a gaussian range
           my_topo <- max_elev*exp((-(seq(1,nodes,1)-peak_loc)^2) / (2*(hhhw^2)))
           if(lee_plateau==TRUE & is.numeric(plateau_elev)){
             if(plateau_elev < 0 | plateau_elev > max_elev){
               stop('Broken topography! Plateau elevation cannot be < 0 nor > max elevation')
             }
             for(i in 1:nodes){
               if(i < peak_loc){next
               } else if(i >= peak_loc & my_topo[i] > plateau_elev){ next
               } else{my_topo[i] <- plateau_elev}
             }
           } else{stop("if lee_plateau=TRUE, plateau_elev must be defined")}
           },
         {# if we're building a witch of agnesi range
           symlim <- min(nodes-peak_loc, peak_loc)   # symmetry limit (min distance to domain edge)
           woar_z <- vector(length=(symlim*2))
           woar <- max_elev * ((hhhw^2) / ((seq(1,symlim,1)^2) + (hhhw^2)))   # woar as fxn of distance from peak location
           woar_z[(peak_loc+1):(peak_loc+symlim)] <- woar    # topography to right of ridge center
           woar_z[peak_loc:(peak_loc-symlim)] <- woar        # topography to left of ridge center
           z <- rep(NA,nodes) ; z[1:length(woar_z)] <- woar_z    # add topography to the domain, with NAs where the ridge isn't
           if((nodes-peak_loc) > peak_loc){z[nodes] <- 0}   # set boundary of 0m at one end of domain (further from woar peak)
           if((nodes-peak_loc) < peak_loc){z[1] <- 0}
           my_topo <- na.interpolation(z)   # interpolate from woar to boundary to finish the domain
           #... build a plateau if so desired
           if(lee_plateau==TRUE & is.numeric(plateau_elev)){
             if(plateau_elev < 0 | plateau_elev > max_elev){
               stop('Broken topography! Plateau elevation cannot be < 0 nor > max elevation')
             }
             for(i in 1:nodes){
               if(i < peak_loc){next
               } else if(i >= peak_loc & my_topo[i] > plateau_elev){ next
               } else{my_topo[i] <- plateau_elev}
             }
           } else{stop("if lee_plateau=TRUE, plateau_elev must be defined")}
           #... end of the switch
           })
  
  #... return the designed topography
  return(my_topo)
}

