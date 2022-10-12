#' -----------------------------------------------------------------------------------
#' DTVGM_PML.R
#' 
#' Functions for hydrological simulation based on the DTVGM-PML 
#' the Distributed Time Variant Gain Model coupled with Penman–Monteith–Leuning (DTVGM-PML)
#' Models include HBV_snowmelt, Gash_interception, PML_evapotranspiration, DTVGM_runoff
#'  
#' Created by Songzh, 2022-10-10
#' 
#' @references 
#' [1] Song, Z., Xia, J., Wang, G., She, D., Hu, C., Hong, S., 2022. Regionalization 
#'     of hydrological model parameters using gradient boosting machine. Hydrology 
#'     and Earth System Sciences 26, 505–524. https://doi.org/10.5194/hess-26-505-2022
#' [2] Zhang, Y., Kong, D., Gan, R., Chiew, F.H.S., McVicar, T.R., Zhang, Q., Yang, Y., 2019. 
#'     Coupled estimation of 500 m and 8-day resolution global evapotranspiration and gross 
#'     primary production in 2002–2017. Remote Sensing of Environment 222, 165–182. 
#'     https://doi.org/10.1016/j.rse.2018.12.031
#' [3] https://github.com/gee-hydro/gee_PML
#' -----------------------------------------------------------------------------------

#' DTVGM_PML
#' Hydrological simulation by DTVGM-PML
#' @param data   data.frame, forcing data including prec, temp, pres, shum, wind, srad, lrad, lai, alb
#' @param para   list, model parameters fer0	sls	hc	gsx	cmelt	g1	g2	ks	kr	kg	wm
#' @param W_init double, Degree-day factor, a melt factor of daily snowmelt when warm enough to melt.
#' @param G_init double, Initial snowpack 
#' @return list, water balance variables including 'P', 'E', 'R', 'P2', 'Rs', 'Rss', 'Rg', 'Ei', 'Ec', 'Es', 'W', 'G', 'PET'  
#' DTVGM_PML
DTVGM_PML <- function(data, para, W_init = 0, G_init = 0) {
  # initialization 
  prec <- data$prec
  temp <- data$temp
  pres <- data$pres
  shum <- data$shum
  wind <- data$wind
  srad <- data$srad
  lrad <- data$lrad
  lai <- data$lai
  alb <- data$alb
  
  fer0 <- para$fer0
  sls <- para$sls
  hc <- para$hc
  gsx <- para$gsx
  cmelt <- para$cmelt
  g1 <- para$g1
  g2 <- para$g2
  ks <- para$ks
  kr <- para$kr
  kg <- para$kg
  wm <- para$wm
  
  n <- nrow(data)  # data length
  W <- G <- vector(length = n + 1)  # soil moisture
  W[1] <- W_init
  G[1] <- G_init
  Rs <- Rss <- Rg <- vector(length = n)  # runoff 
  Es <- vector(length = n)  # soil evaporation 
  
  # snowmelt routine
  # separating precipitation
  tt <- 0  # temperature threshold 
  pr <- ifelse(temp < tt, 0, prec)    # rainfall
  ps <- ifelse(temp < tt, prec, 0)    # snowfall
  sm <- getSnowmelt(ps, temp, cmelt, tt = tt)  # snowmelt
  
  # interception routine
  Ei <- getInterceptEvap(pr, lai, fer0, sls)
  
  # actual evapotranspiration routine
  e <- getEvapotranspiration(temp, pres, shum, wind, srad, lrad, lai, alb, gsx, hc)
  Ec <- e$Ec
  Es_eq <- e$Es_eq
  PET <- e$PET
  
  # runoff routine
  pe <- pr - Ei + sm  # rainfall and snowmelt down into soil
  for (i in 1:n) {
    # soil evaporation
    fval_soil <- W[i] / wm
    Es[i] <- fval_soil * Es_eq[i]
    et <- Es[i] + Ec[i]  # total evapotranspiration
    
    # surface runoff
    g <- g1 * fval_soil ^ g2
    g <- min(g, 1)  # ifelse(g > 1, 1, g)
    # Rs[i] <- ifelse(pe[i] < et, 0, pe[i] * g)
    if (pe[i] < et)
      Rs[i] <- 0
    else
      Rs[i] <- pe[i] * g
    
    inf <- pe[i] - Rs[i]                  # infiltration
    Rss[i] <- ks * inf * fval_soil        # interflow
    Rg[i] <- kg * G[i]                    # base flow
    Gr <- kr * fval_soil * (inf - Rss[i]) # groundwater recharge
    Sr <- inf - Rss[i] - Gr               # soil moisture storage recharge
    W[i + 1] <- W[i] + Sr - et            # soil moisture storage change
    G[i + 1] <-  G[i] + Gr - Rg[i]        # groundwater storage change
    
    if (W[i + 1] < 0) {
      G[i + 1] <- G[i + 1] + W[i + 1]
      W[i + 1] <- 0
      if (G[i + 1] < 0) {
        Es[i] <- Es[i] * (1 +  G[i + 1] / et)
        Ec[i] <- Ec[i] * (1 +  G[i + 1] / et)
        G[i + 1] <- 0
      }
    } else if (W[i + 1] > wm) {
      Rss[i] <- Rss[i] + (W[i + 1] - wm)
      W[i + 1] <- wm
    }
  }
  
  R <- Rs + Rss + Rg  # total runoff
  E <- Ei + Ec + Es   # total evapotranspiration 
  
  out <- data.frame(P=prec, E, R, P2=pr+sm, Rs, Rss, Rg, Ei, Ec, Es, W=W[1:n], G=G[1:n], PET)
  return(out)
}


#' snowMelt
#' HBV snowmelt model
#' @param ps       vector, Snowfall (mm)
#' @param temp     vector, Temperatures (deg C)
#' @param cmelt    double, Degree-day factor, a melt factor of daily snowmelt when warm enough to melt.
#' @param spk_init double, Initial snowpack 
#' @param rsm_init double, Initial retained snowmelt in snowpack
#' @param tt       double, Threshold temperature (deg C)
#' @param cwh      double, Retain limits of snowmelt in snowpack
#' @param cfr      double, Refreezing coefficient
#' @return vector, Snowmelt down into soil 
#' snowMelt()
getSnowmelt <- function(ps, temp, cmelt, spk_init=0, rsm_init=0, tt=0, cwh=0.1, cfr=0.05) {
  n <- length(ps)
  spk <- rsm <- vector(length = n + 1)
  smtd <- vector(length = n)  # snowmelt down to soil moisture routine
  
  spk[1] <- spk_init
  rsm[1] <- rsm_init

  for (i in 1:n) {
    if (temp[i] > tt) {  # snowmelt
      smt <- min(spk[i] + ps[i], cmelt * (temp[i] - tt))
      spk[i + 1] <- spk[i] + ps[i] - smt
      rsm[i + 1] <- rsm[i] + smt
      # smtd[i] <- ifelse(rsm[i + 1] > cwh * spk[i + 1], rsm[i + 1] - cwh * spk[i + 1],  0)
      tmp <- rsm[i + 1] - cwh * spk[i + 1]
      if (tmp > 0)
        smtd[i] <- tmp
      else
        smtd[i] <- 0
      rsm[i+1] <- rsm[i+1] - smtd[i]
    } else {  # refreezing meltwater
      rmw <- min(rsm[i], cfr * cmelt * (tt - temp[i]))
      rsm[i + 1] <- rsm[i] - rmw
      spk[i + 1] <- spk[i] + ps[i] + rmw
      smtd[i] <- 0
    }
  }
  return(smtd)
}


#' interceptEvap
#' Gash interception model
#' @param pr     vector, Rainfall (mm)
#' @param lai    vector, Leaf area index
#' @param fER0   double, Average ratio of wet canopy evaporation rate and rainfall rate for full canopy cover
#' @param s_sls  double, Canopy storage capacity per unit leaf area (mm)
#' @param lairef double, Reference LAI determining canopy cover, setting as 5
#' @return vector, Intercept evaporation  
#' interceptEvap()
getInterceptEvap <- function(pr, lai, fER0, s_sls, lairef = 5) {
  fveg <- 1 - exp(-lai / lairef)
  sveg <- s_sls * lai
  fER <- fveg * fER0
  prec_wet <- -log(1 - fER0) / fER * sveg
  Ei <- ifelse(pr < prec_wet, fveg * pr, fveg * prec_wet + fER * (pr - prec_wet))
  Ei <- ifelse(is.na(Ei), 0, Ei)
  return(Ei)
}


#' getEvapotranspiration
#' PML (Penman-Monteith-Leuning) model 
#' Reference to https://github.com/gee-hydro/gee_PML
#' @param temp vector, Temperature (deg C)
#' @param pres vector, Pressure (KPa)
#' @param shum vector, Specific humidity (kg/kg)
#' @param wind vector, Wind speed at 10m height (m/s)
#' @param srad vector, Downward short-wave radiation (MJ/m2/day)
#' @param lrad vector, Downward long-wave radiation (MJ/m2/day)
#' @param lai  vector, Leaf area index
#' @param alb  vector, Albedo
#' @param gsx  double, Maximum stomatal conductance of leaves at the top of the canopy [0.003-0.01] (m/s)
#' @param hc   double, Canopy height
#' @return list, Transpiration, soil evaporation at equilibrium, and potential evapotranspitation
#' The soil evaporation is calculated as fval_soil * Es_eq with fval_soil = W/Wm in DTVGM water balance procedure
#' getEvapotranspiration()
getEvapotranspiration <- function(temp, pres, shum, wind, srad, lrad, lai, alb, gsx, hc) {
  # Constants 
  sigma <- 4.9e-9 # Stefan-Boltzmann constant [MJ K-4 m-2 day-1].
  kmar <- 0.4     # von Karman's constant
  Zob <- 10       # m, making sure higher than hc, reference height where μ is measured (m).
  Cp <- 0.001013  # specific heat at constant pressure  [MJ kg-1 K-1] equal to [MJ kg-1 ℃-1]
  kQ <- 0.4488    # extinction coefficient
  kA <- 0.7       # the attenuation of net all-wave irradicance, typically about 0.6-0.8 (Denmend, 1976, Kelliher FM et al., (1995))
  Q50 <- 2.6      # the value of absorbed PAR when gs<-gsx/2, [MJ m-2 day-1]
  D50 <- 0.8      # the humidity deficit at which stomatal conductance is half its maximum value [kpa] 
  lai_ref <- 5
  
  # Atmosphere parameters
  lambda <- 2.501 - 2.361e-3 * temp  # Latent heat of vaporization
  va <- shum * pres / (0.622 + 0.378 * shum)  # actual vapor pressure
  vs <- 0.6108 * exp(17.27 * temp / (temp + 237.3))  # saturated vapour pressure
  vpd <- vs - va
  rou_a <- 1.293 * 273 / (273 + temp) * pres / 101.3  # density of air [kg/m3]
  gamma <- 0.00163 * pres / lambda  # psychrometric constant (S2.9)
  delta <- 4098 * vs / ((temp + 237.3)^2)  # slope of vapour pressure curve (S2.4)
  
  # Radiation calculation
  fv <- 1 - exp(-lai / lai_ref)       # vegetation fractional cover
  se <- 0.97 * fv + 0.92 * (1 - fv)   # surface emissivity
  Rns  <- srad * (1 - alb)
  Rnl <- lrad - se * sigma * (temp + 273.2)^4
  Rn <- Rns + Rnl
  Rn[Rn<0] <- 0
  
  # estimate daily FAO-56 reference crop evapotranspiration
  u2 <- wind * 4.87 / (log(67.8 * Zob - 5.42))
  pet <- (0.408 * delta * Rn + gamma * 900 / (temp + 273) * u2 * vpd) / (delta+gamma*(1+0.34*u2))
  
  # PML procedure 
  Qh <- srad * 0.45 # flux density of visible radiation at the top of the canopy (approximately half of incoming solar radiation)
  
  Gc <- gsx / kQ * log((Qh + Q50) / (Qh * exp(-kQ * lai) + Q50)) / (1 + delta / D50) # canopy conductance (m/s)
  Gc[Gc < 1e-6] <- 1e-6
  # 
  # hc <- 8 # canopy height
  d <- hc * 2 / 3 # zero plane displacement height (m)
  zom <- 0.123 * hc # roughness lengths governing transfer
  zov <- 0.1 * zom # of momentum and water vapor (m)
  # uz <- log(67.8 * Zob - 5.42) / 4.85 * u2 # wind speed at height zr (m/s)
  uz <- wind
  Ga <- kmar * kmar * uz / (log((Zob - d) / zom) * log((Zob - d) / zov)) # aerodynamic conductance (m/s)
  Ga[Ga > 0.033] <- 0.033
  
  # the fraction of the total available energy (Rn - G) that is absorbed by the soil surface
  tou <- exp(-kA * lai)
  # Transpiration from plant cause by radiation water transfer
  LEcr <- delta / gamma * Rn * (1 - tou) / (delta / gamma + 1 + Ga / Gc)
  # Transpiration from plant cause by aerodynamic water transfer
  LEca <- rou_a * Cp * vpd * Ga / gamma / (delta / gamma + 1 + Ga / Gc)
  LEc <- LEcr + LEca
  # soil evaporation at equilibrium 
  LEs_eq <- delta * Rn * tou / (delta + gamma)
  # evaporation components [mm] [MJ m-2 day-1] / [MJ kg-1] <- [kg m-2 day-1]
  Es_eq <- LEs_eq / lambda
  Ecr <- LEcr / lambda
  Eca <- LEca / lambda
  Ec <- LEc / lambda
  # Es <- Es_eq * fval_soil
  return(list(Ec = Ec, Es_eq = Es_eq, PET = pet))
}
