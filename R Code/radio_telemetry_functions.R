antenna_angle_seq <- function(input_angle, num_antennas){
  ant_ang_inc <- 360/num_antennas
  start_angle <- input_angle %% ant_ang_inc
  return(seq(from = start_angle, by = ant_ang_inc, length.out = num_antennas))
}
# encoded_signal_to_actual_watts <- function(Z, p0 = 4.3458*10^-11, b = 0.3012) {
encoded_signal_to_actual_watts <- function(Z, p0 = 4.89*10^-11, b = 0.3013) {
  #convert Z integer value (0-255) to an actual signal strength in watts
  return((exp(atanh(Z/255)/b)-1) * p0)

}

encoded_signal_to_actual_dbm <- function(Z, p0 = 4.89*10^-11, b = 0.3013) {
  #convert Z integer value (0-255) to an actual signal strength in dBm
  signal_watts <- encoded_signal_to_actual_watts(Z, p0, b)
  return(watts_to_dbm(signal_watts))
}

watts_to_dbm <- function(power_watts){
  return(10*log10(power_watts * 1000))
}

dbm_to_watts <- function(power_dbm){
  return((10^(power_dbm/10))/1000)
}

dbm_to_mw <- function(power_dbm){
  return(10^(power_dbm/10))
}


psi_calc <- function(x, y, XT, YT, theta){
  # function to calculate psi for use in antenna strength modeling
  require(pracma)
  degree <- atan2d((y-YT),(x-XT)) # the degrees of the x,y point
  if (degree < 0) {degree <- degree+360}
  psi <- degree - theta  #subtract the antenna angle (theta)
  return(deg2rad(psi)) # radians
}

xi_square_calc_dbm_r_psi <- function(r, psi, z, ant_HT, min_xi_dbm, lambda, D0, p0){
  D <- 10^(D0/10)  #Where D0 is the abosolute antenna gain and D = 7.22*Le/lambda
  le <- D*lambda/7.22
  
  #default antenna params from Janaswamy et al 2018
  k0 <- 2*pi/lambda # wavenumber in freespace
  beta0 <- -(k0+2.94/le) 
  p <- beta0*le/2 
  q <- k0*le/2
  
  # r=radial distance in meters
  # The direct line-of-sight distance between the tower and the bird in m
  R <- sqrt(r^2 + (z - ant_HT)^2) 
  # psi <- angle in radians around antenna
  g_psi <- cos(pi/2*sin(psi))/cos(psi) * (sin(p+q*cos(psi))/(p+q*cos(psi))) #version from Janaswamy - eqn 19 - Antenna field pattern
  
  xi_square <- (g_psi*sin(k0*ant_HT*z/R)/(k0*R))^2  #derived from eqn 40 in  Janawany
  xi_square_dbm <- watts_to_dbm(xi_square)
  
  return(xi_square_dbm)
  
}

#original below for 9 element yagi, default to 5
# xi_square_calc_dbm <- function(x, y, z, XT, YT, HT, theta, lambda=1.8, D0=11.1, p0=4.3458*10^-11, noisy = F){
  
xi_square_calc_dbm <- function(x, y, z, XT, YT, HT, theta, lambda, D0, p0, noisy = F){
  # Function for calculating receiver power at an x,y,z coordinate in dBm
  # 3D coordinates of the animal = x, y, z
  # coordinates of tower and angle of antenna = XT, YT, HT, theta
  # Antenna and receiver specific info = p,q,k0, p0
  # noise value = mu and noisy is whether to include noise or not
  # (r, phi) are the polar coordinates of the point (x,y)
  #lambda = wavelength in free space
  #le the effective length of the overall array
  #from email from R. Janaswamy on 25 Jun 2020
  #   Here is the formula I suggest:
  #     (i) Let D = absolute gain of antenna, D0 = dB gain of the antenna: D = 10^{D0/10). In your case D0 = 9 dB. Then D = 7.94
  #     (ii) D = 7.22*Le/lambda, lambda = wavelength = 1.8m in your case, Le = effective length. Using the value of D from (i) you would get Le = D*lambda/7.22 = approx D/4 = approx 2m. 
  #     (iii) Use this value of Le in equation (19) of my paper.
  #     
  #     IN response to from P Paton.
  #     1) How is the antenna gain included in the equations when estimating the radiation pattern (in the findinters.m script)? I see that Eq. 19 in Rama's paper 
  #        has been used for the radiation pattern (g(psi)) which is the function of p,q, and psi. On the other hand, p, q, and beta0 are also functions of effective 
  #      length (Le) and wave number (K0). As K0 is constant, does it mean that the effect of antenna gain has only been included by the effective length? Then, what about the number of elements on the Yagi?
  #     2) If that is correct, let's say we want to revise the code for a 5-element yagi with 9 dB, should we just replace the effective length for that? (it does not make sense to me)
  # 11 Dec 2020 - default set to 5 element yagi at 166 MHZ,  p0=4.89*10^-11 from 2017 calibration
  
  D <- 10^(D0/10)  #Where D0 is the abosolute antenna gain and D = 7.22*Le/lambda
  le <- D*lambda/7.22
  
  #default antenna params from Janaswamy et al 2018
  k0 <- 2*pi/lambda # wavenumber in freespace
  beta0 <- -(k0+2.94/le) 
  p <- beta0*le/2 
  q <- k0*le/2
  
  # radial distance in meters
  r <- sqrt((x-XT)^2 + (y-YT)^2)
  # The direct line-of-sight distance between the tower and the bird in m
  R <- sqrt(r^2 + (z - HT)^2) 
  # psi <- angle in radians around antenna
  psi <- psi_calc(x,y,XT,YT, theta)
  g_psi <- cos(pi/2*sin(psi))/cos(psi) * (sin(p+q*cos(psi))/(p+q*cos(psi))) #version from Janaswamy - eqn 19
  
  # xi_square <- (r1_g_psi*sin(k0*HT*z/R)/(k0*R))^2 # 2.2 in Hua Bai in his code he had k0*HT*z^2 and
  xi_square <- (g_psi*sin(k0*HT*z/R)/(k0*R))^2  #derived from eqn 40 in  Janawany
  xi_square_dbm <- watts_to_dbm(xi_square)
  
  xi_square_noisy <- (sqrt(xi_square) + sqrt(p0))^2  #+ 2*sqrt(xi_square)*sqrt(p0) + p0*mu^2eqn 27 Ranasmamy
  xi_square_noisy_dbm <- watts_to_dbm(xi_square_noisy)
  
  if (noisy==F) {return(xi_square_dbm)
  } else {return(xi_square_noisy_dbm)}
  
}


r_squared <- function(x, y, z, XT, YT, HT, theta, lambda){
  k0 <- 2*pi/lambda # wavenumber in freespace
  r <- sqrt((x-XT)^2 + (y-YT)^2)
  R <- sqrt(r^2 + (z - HT)^2)
  psi <- psi_calc(x,y,XT,YT, theta)
  g_psi <- cos(pi/2*sin(psi))/cos(psi) * (sin(p+q*cos(psi))/(p+q*cos(psi))) #version from Janaswamy - eqn 19)
  xi <- g_psi*sin(k0*HT*z/R)/(k0*R)
  return((HT*z/abs(xi)*abs(g_psi)))
}

field_radiation_pattern <- function(lambda, D0, xi_min_dbm){
  D <- 10^(D0/10)  #Where D0 is the absolute antenna gain and D = 7.22*Le/lambda
  le <- D*lambda/7.22
  
  #default antenna params from Janaswamy et al 2018
  k0 <- 2*pi/lambda # wavenumber in freespace
  beta0 <- -(k0+2.94/le) 
  p <- beta0*le/2 
  q <- k0*le/2
  
  min_xi_sq_watts <- dbm_to_watts(xi_min_dbm)
  
  gpsi_df <- data.frame()
  for (psi in seq(0,(2*pi), pi/1000)){
    #from eqn 41 Ratnaswamy to approximate the field radiation pattern at a given power
    g_psi <- cos(pi/2*sin(psi))/cos(psi) * (sin(p+q*cos(psi))/(p+q*cos(psi))) #version from Janaswamy - eqn 19
    R <- g_psi/(k0*sqrt(min_xi_sq_watts))
    gpsi_df <- rbind(gpsi_df, cbind(psi=psi, g_psi=g_psi, R=R))
    
  }
  
  return(gpsi_df)
}

detection_limit_distance <- function(z_in, ant_HT, min_xi_dbm, x_max_end = 10000, inc = 10000, seq_inc = 10) { #, x_min_end = -10000
  #determine the distance necessary to estimate along the main antenna lobe 
  #problem is nodes of low antenna strength I need to get across and also issues with rear lobes being at angles extending backwards past 0 axis
  #just calc main lobe at front since this is better behaved. 
  #z is bird height
  #ant_X, ant_Y, ant_HT, ant_theta are antenna params
  #min_xi_dbm is min power for detection
  #x max direction
  x_max_start = 0
  continue_proc = T
  while(continue_proc){
    power_vals <- data.frame()
    for (x_val in seq(x_max_start,x_max_end, seq_inc)){
      power_dbm = xi_square_calc_dbm(x=x_val, y=0, z=z_in, XT=0, YT=0, HT=ant_HT, theta=0)
      power_vals <- rbind(power_vals, cbind("x_val"=x_val,"power_dbm"=power_dbm))
    }
    # print(power_vals)
    x_max_detect_vals <- power_vals[which(power_vals$power_dbm >= min_xi_dbm), ]
    x_max_detect <- max(x_max_detect_vals$x_val) #determine max dist detected
    # print(x_max_detect)
    if (x_max_detect < x_max_end) {
      continue_proc = F
    } else {
      x_max_start <- x_max_end  #reset start
      x_max_end <- x_max_end + inc  #increment end in case we haven't reached max
    }
  }
}

bisect_root <- function (f, a, b, num = 10, eps = 1e-05) {
    h = abs(b - a)/num
    i = 0
    j = 0
    a1 = b1 = 0
    while (i <= num) {
      a1 = a + i * h
      b1 = a1 + h
      if (f(a1) == 0) {
        # print(a1)
        # print(f(a1))
        val <- a1
        # print("A")
      }
      else if (f(b1) == 0) {
        # print(b1)
        # print(f(b1))
        val <- b1
        # print("B")
      }
      else if (f(a1) * f(b1) < 0) {
        repeat {
          if (abs(b1 - a1) < eps) 
            break
          x <- (a1 + b1)/2
          if (f(a1) * f(x) < 0) 
            b1 <- x
          else a1 <- x
        }
        # print(j + 1)
        j = j + 1
        # print("C")
        val <- (a1 + b1)/2
        # print(val)
        # print(f(val))
      }
      i = i + 1
    }
    if (j == 0) 
    {
      # print("finding root is fail")
      return(NA)}
    else {
      # print("finding root is successful")
      return(val)}
  }
  
  
anntenna_detect_pattern_root <- function(z, ant_HT, xi_min_dbm, lambda, D0, p0, maxit=10000){
  require(sf)
  D <- 10^(D0/10)  #Where D0 is the abosolute antenna gain and D = 7.22*Le/lambda
  le <- D*lambda/7.22
  
  #default antenna params from Janaswamy et al 2018
  k0 <- 2*pi/lambda # wavenumber in freespace
  beta0 <- -(k0+2.94/le) 
  p <- beta0*le/2 
  q <- k0*le/2
  
  # xi_square <- (g_psi*sin(k0*HT*z/R)/(k0*R))^2  #derived from eqn 40 in  Janawany
  min_xi_watts <- dbm_to_watts(xi_min_dbm)  #convert dbm to watts
  
  root_R_fx=function(R) {
    (g_psi*sin(k0*ant_HT*z/R)/(k0*R))^2-min_xi_watts
  }
  
  field_radiation_pattern <- data.frame()
  for (psi in seq(0,(2*pi), pi/1000)){
    #from eqn 41 Ratnaswamy to approximate the field radiation pattern at a given power
    #rootSolve deals with finding the roots of n nonlinear (or linear) equations
    g_psi <- cos(pi/2*sin(psi))/cos(psi) * (sin(p+q*cos(psi))/(p+q*cos(psi))) #version from Janaswamy - eqn 19
    #The root of a function f(x) is the value of x for which f(x) = 0
    R1 <- max(rootSolve::uniroot.all(f=root_R_fx, interval=c(1,100000), maxiter=maxit))
    # R <- uniroot(f=root_R_fx, interval=c(10,1000),maxiter=maxit, extendInt="downX", tol = 1e-7)$root
    R2 <- bisect_root(f = root_R_fx, a = 1, b = 100000)
    # R3 <- pracma::brent(fun = root_R_fx, a=-1, b=100000, maxiter = 1000)$root
    R <- max(c(R1, R2), na.rm = T)
    
    # R <- R1
    # R <- R2
    field_radiation_pattern <- rbind(field_radiation_pattern, cbind(psi=psi, R=R))
  }
  
  XT = 0; YT = 0
  
  field_radiation_pattern$x = field_radiation_pattern$R*cos(field_radiation_pattern$psi) + XT
  field_radiation_pattern$y = field_radiation_pattern$R*sin(field_radiation_pattern$psi) + YT
  field_radiation_pattern$stn_ht=ant_HT
  field_radiation_pattern$z=z
  #Some floating point math? issues where  erratic points around x=0, clean
  field_radiation_pattern <- field_radiation_pattern[which(round(field_radiation_pattern$x, 5)!=0),]
  field_radiation_pattern[nrow(field_radiation_pattern)+1,] <- field_radiation_pattern[1,]
  field_radiation_pattern <- field_radiation_pattern[which(!is.infinite(field_radiation_pattern$R)),]
  field_radiation_pattern_sf <-  st_sfc(st_polygon(list(as.matrix(field_radiation_pattern[,c("x", "y")]))))
  
  return(field_radiation_pattern_sf)
}


rotation = function(a){
  r = a * pi / 180 #degrees to radians
  # matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)
  #rotation about axis believe require inverse
  matrix(c(cos(r), -sin(r), sin(r), cos(r)), nrow = 2, ncol = 2)  
  
} 

multi_antenna_pattern <- function(ant_stations, stn_id, ant_angle_st=60, antenna_num, ant_X=0, ant_Y=0, ant_HT, z, min_xi,interval=100, 
                                  single_antenna_detect_poly=NULL, crs_locs = 3857, crs_out=NULL, noisy = F){
  #generate multi-antenna detection pattern for mapping based on single antenna estimate, rotated and shifted from x,y=0,0 and theta=0
  WebMerc <- CRS("+init=epsg:3857")
  # #if single antenna detection poly not passed to function, create the poly pattern
  # if(is.null(single_antenna_detect_poly)){
  #   # min_max_lims <- detection_limit_distance(z_in=z, ant_HT=ant_HT, min_xi_dbm=min_xi, x_max_end = 75000, x_min_end = -15000, inc = 10000, seq_inc = 50)
  #   xmax <- detection_limit_distance(z_in=z, ant_HT=ant_HT, min_xi_dbm=min_xi, x_max_end = 75000, inc = 10000, seq_inc = 50)
  #   #calc single antenna pattern to replicate to all other positions
  #   single_antenna_detect_poly <- antenna_detect_patterns(x_min=-1*xmax, x_max=xmax, y_min=-1*xmax, y_max=1*xmax, z=z, interval=100, 
  #                                                         ant_X=ant_X, ant_Y=ant_Y, ant_HT=ant_HT, ant_theta=0, min_xi = min_xi, crs_locs = 3857, noisy = noisy)}
  #generate angles from start angle
  # browser()
  antenna_angles <- seq(ant_angle_st, 360, 360/antenna_num)
  ant_stations <- spTransform(ant_stations, CRSobj = WebMerc)
  antenna_array_detect_polys <- data.frame()
  for (i in seq(1,length(ant_stations),1)) {
    #get each antenna station
    antenna <- ant_stations[i,]
    #cycle through each location and antenna placement
    # antenna_x <- antenna@coords[1]
    # antenna_y <- antenna@coords[2]
    
    for (j in (1:antenna_num)) {
      theta <- antenna_angles[j]
      #rotate about theta
      single_antenna_rotated <- single_antenna_detect_poly * rotation(theta)
      
      #affine shift to location of antenna
      single_antenna_rotated <- single_antenna_rotated + antenna@coords[,1:2]
      antenna_array_detect_polys <- rbind(antenna_array_detect_polys, cbind("station" = stn_id, "theta"=theta, "geometry"=single_antenna_rotated))
    }
  }
  
  multiple_antennas_detect_polys <- st_sf(antenna_array_detect_polys, crs=crs_locs)
  if(!is.null(crs_out)){multiple_antennas_detect_polys <- st_transform(multiple_antennas_detect_polys, crs=4326)}
  return(multiple_antennas_detect_polys) 
  
}

antenna_angle_optim_effecient <- function(proposed_stn_points, n_antennas, ant_ang_inc, single_antenna_optim, detect_locs, debug.out){
  require(sf)
  require(data.table)
  require(nngeo)
  require(dplyr)
  # 1) determine nearest station pairs
  # 2) determine angles of least overlap for adjacent stations
  # 3) Combine all min. overlap station to find optimal solution
  # 4) look at angles +/- an angle increment to find best
  
  n_stations <- length(proposed_stn_points)
  proposed_stn_points_sf <- st_transform(st_as_sf(proposed_stn_points), WebMerc)
  #get closest 3 points (k=3)
  near_pts <- nngeo::st_nn(proposed_stn_points_sf, proposed_stn_points_sf, k=3, returnDist = T, progress=F)  

  near_pts_idx <- rbindlist(lapply(near_pts$nn, function(x) as.data.table(cbind(x[2:3])))) #get index of the second/third closest (first is itself)
  near_pts_dist <- rbindlist(lapply(near_pts$dist, function(x) as.data.table(cbind(x[2:3])))) #get distance of the second/third closest (first is itself)
  
  point_pairs <- data.frame(stn1=sort(rep(1:n_stations, 2)))
  point_pairs <- cbind(point_pairs, near_pts_idx)
  point_pairs <- cbind(point_pairs, near_pts_dist)
  colnames(point_pairs) <- c("stn1","stn2","stn1_2_dist")
  point_pairs$pairs <- apply(point_pairs, 1, function(x) paste(sort(c(x[["stn1"]], x[["stn2"]])), collapse = ","))
  
  #Since pairs of points - can remove half because closest pair already found
  point_pairs_no_dups <- point_pairs[!duplicated(point_pairs[,c("pairs", "stn1_2_dist")]),] #remove dups
  point_pairs_no_dups <- point_pairs_no_dups[order(point_pairs_no_dups$stn1_2_dist),]
  #select top n stations -1 for comparison
  # point_pairs_optim <-  point_pairs_no_dups[1:(n_stations-1),]
  point_pairs_optim <-  point_pairs_no_dups[1:n_stations,]
  
  # st_area(st_intersection())
  ant_angle_seq <- seq(0,floor(360/n_antennas), ant_ang_inc)
  n_ant_ang_seq <- length(ant_angle_seq)
  
  #get set of all possible antenna angle permutations to test between pairs of stations
  ant_perms <- as.data.frame(gtools::permutations(n_ant_ang_seq, 2, ant_angle_seq, ant_ang_inc, repeats.allowed = T))

  ###JUST CHANGED BELOW on 11.2.20 due to different CRS error only in shinyapps.io
  points_to_detect <- st_transform(detect_locs, crs=3857)
  # points_to_detect <- detect_locs
  
  #cycle through list of proposed station pairs to get best coverage for pairs
  station_pair_list=list()
  for (i in 1:nrow(point_pairs_optim)) {
    max_detected = 0
    stn1_val <- point_pairs_optim[i, ]$stn1
    stn2_val <- point_pairs_optim[i, ]$stn2
    pt1 <- proposed_stn_points[stn1_val,]
    pt2 <- proposed_stn_points[stn2_val,]
    close_stns <- maptools::spRbind(pt1, pt2)

    for (k in 1:nrow(ant_perms)){ 
      angle_set <- ant_perms[k,]  
      #build antenna stations based on angle set and determine how many point covered
      #step through each station to generate different pattern and combine
      station_set <- list()
      for (j in 1:ncol(angle_set)){
        ant_angle_j <- angle_set[,j] 
        stn <- close_stns[j,]
        if (j==1) {station_set[[paste0("stn_", j)]]$station <- stn1_val}
        else {station_set[[paste0("stn_", j)]]$station <- stn2_val}
        station_set[[paste0("stn_", j)]] <- multi_antenna_pattern(ant_stations = stn, stn_id = station_set[[paste0("stn_", j)]]$station,ant_angle_st = ant_angle_j, 
                                                                  antenna_num = n_antennas, single_antenna_detect_poly = single_antenna_optim)

      }

      station_set <- data.table::rbindlist(station_set)

      #combine polys to for intersection of one surface, otherwise given more points than possible because intersection by each poly
      detected_points <- sum(st_intersects(st_union(station_set$geometry), points_to_detect, sparse = F)) #total number of detected points for this array

      if (max_detected < detected_points) {
        max_detected <- detected_points
        max_angles <- station_set
      }
      
      # detected_points_sum <- sum(detected_points)  #total number of detected points for this array
      debug.out$debug <- renderPrint(paste("TO HERE D7", str(max_detected)))
      
    }
    # return(build_stations)
    station_pair_list[[i]] <- list("max_detected"=max_detected, "max_angles"=max_angles)
  }
  # station_pair_list_combine <- data.frame()
  station_pair_list_combine <- lapply(station_pair_list, function(x) {
    tbl <- as.data.frame(x$max_angles)
    # station_pair_list_combine <- rbind(station_pair_list_combine, tbl)
  })
  
  station_pair_list_combine <- data.table::rbindlist(station_pair_list_combine)
  station_pair_list_combine$station <- unlist(station_pair_list_combine$station)
  station_pair_list_combine$theta <- unlist(station_pair_list_combine$theta)
  
  #remove full dups first
    station_pair_list_combine <- unique(station_pair_list_combine, by=c("station", "theta"))
  station_pair_list_combine$id <- unlist(lapply(1:(nrow(station_pair_list_combine)/n_antennas), function(x) rep(x, n_antennas)))
  #determine which stations are duplicate and use if same or use mean (rounded to nearest 15) for angles

  #generates optimal list of start angles from which we can get pattern
  stn_dups <- station_pair_list_combine %>% 
    group_by(station, id) %>% 
    summarise(theta=min(theta)) %>% 
    group_by(station) %>%
    summarise(theta=DescTools::RoundTo(mean(theta), 15)) %>% #round to 15 multiple 
    mutate(theta_min=theta-15, theta_plus=theta+15) 
  
  #generate combos of angles
  stn_dups_list <- as.list(as.data.frame(t(stn_dups[,2:ncol(stn_dups)])))
  stn_dups_exp <- expand.grid(stn_dups_list)
  
  max_detected = 0
  for (i in 1:nrow(stn_dups_exp)) {
    angle_opts <- stn_dups_exp[i,]
    #build antenna stations based on angle set and determine how many point covered
    #step through each station to generate different pattern and combine
    station_set <- list()
    for (j in 1:ncol(angle_opts)){
      ant_angle_j <- angle_opts[,j] 
      stn <- proposed_stn_points[j, ]
      station_set[[paste0("stn_", j)]] <- multi_antenna_pattern(ant_stations = stn, stn_id = station_set[[paste0("stn_", j)]]$station, ant_angle_st = ant_angle_j, 
                                                                antenna_num = n_antennas, single_antenna_detect_poly = single_antenna_optim)
      station_set[[paste0("stn_", j)]]$station <- j
    }
    
    # points_to_detect
    
    station_set <- data.table::rbindlist(station_set)

    #combine polys to for intersection of one surface, otherwise given more points than possible because intersection by each poly
    detected_points <- sum(st_intersects(st_union(station_set$geometry), points_to_detect, sparse = F)) #total number of detected points for this array

    #look for max detected of points to find optimal arrangement record best set if more detected
    if (max_detected < detected_points) {
      max_detected <- detected_points
      ant_optimized_stn_polys_sf <- station_set
    }
  }
  #added below to create regular df from lists
  ant_optimized_stn_polys_sf$station <- unlist(ant_optimized_stn_polys_sf$station)
  ant_optimized_stn_polys_sf$theta <- unlist(ant_optimized_stn_polys_sf$theta)
  # ant_optimized_stn_polys_sf <- sf::st_sf(ant_optimized_stn_polys_sf)
  
  return(ant_optimized_stn_polys_sf)
}

