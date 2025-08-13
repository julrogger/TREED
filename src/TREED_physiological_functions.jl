###########################################################
########## TREED Julia implementation V0.0 ################
########## Function file ##################################
###########################################################

# Status 05.12.2024


### Plant allometry 
plant_allometry = function(; tr, pars)

    D = (tr.H / pars.k_allom2)^(1 / pars.k_allom3)

    CA = pars.k_allom1 * (D^pars.k_rp)
    CA = min(max(CA, pars.CA_min),pars.CA_max)

    SLA =  (2e-4) * (1/0.4763) * 10^(2.25 - 0.5*log10(tr.a_ll * 12)) # Intermediate of new version Schaphoff.

    SA = (tr.C_leaf * SLA)/pars.k_la_sa
    
    LAI = (tr.C_leaf * SLA)/CA
    
    # Account for differences in light interception depending on SLA/a_ll, ranges as in Schaphoff et al. and Zhang et al.
    # if tr.a_ll > 4
    #   light_extinction_coefficient = 0.4
    # elseif tr.seasonality == 0 && tr.a_ll < 1
    #   light_extinction_coefficient = 0.6
    # else
    #   light_extinction_coefficient = 0.5
    # end
    light_extinction_coefficient = 0.6 - (0.2/5) * tr.a_ll
    
    FPC = (1-exp(-light_extinction_coefficient*LAI)) # P * CA # Foliage projective cover: area of ground covered by foliage (aka. how much light is intercepted) (Lambert-Beer law); *P*CA not applied as I don't consider population densities as [individuals ha-1]
    
    C_sapwood = SA * tr.H * pars.WD_sapwood # I calculate C_sapwood even if SA > totalD (which is not realistic), but in that case it accounts for total supporting structure C

    C_stem_total = pars.k_density_intercept * ((CA * tr.H)^pars.k_power) * 0.47 # g of C per individual
    
    if C_sapwood >= C_stem_total
        # C_sapwood = C_stem_total
        C_sapwood = C_sapwood
        C_heartwood = 0
    else  # For big enough tree assuming sapwood remains cylindrical and rest forms a cone (or something similar) around it
        C_sapwood = C_sapwood
        C_heartwood = C_stem_total - C_sapwood
    end
    
    C_fineroot = tr.C_leaf * tr.r_s_r

    C_coarseroot = (C_heartwood + C_sapwood) * 0.25 * tr.r_s_r # C coarseroot scales with woody C, only as soon as tree gets big enuogh to have D > SA

    tr = (
    C_leaf = tr.C_leaf,
    r_s_r = tr.r_s_r,
    a_ll = tr.a_ll,
    H = tr.H,
    seasonality = tr.seasonality,
    Tave_optim = tr.Tave_optim,
    Tmax_optim = tr.Tmax_optim,
    Tmin_optim = tr.Tmin_optim,
    Pave_optim = tr.Pave_optim,
    C_fineroot = C_fineroot,
    C_sapwood = C_sapwood,
    C_coarseroot = C_coarseroot,
    C_heartwood = C_heartwood,
    D = D,
    CA = CA,
    SLA = SLA,
    SA = SA,
    LAI = LAI,
    FPC = FPC
    )

    return(tr)
end

### GPP and AET function
GPP_function = function(; env, tr, pars)
  
  daytime_integrated_rsds = env.rsds_monthly # * (24/env.daylength) # [MJ m-2 daytime-1/day-1]
  PAR = daytime_integrated_rsds * 4.6 * 0.5 # Photosynthetically active radiation, Schlaphoff version 
  
  # Climatic niche change stress
  Tmean_diff = abs(env.tair_annual - tr.Tave_optim)
  if maximum(env.tair_monthly) > tr.Tmax_optim
    Tmax_diff = abs(maximum(env.tair_monthly) - tr.Tmax_optim)
  else
    Tmax_diff = 0
  end

  if minimum(env.tair_monthly) < tr.Tmin_optim
    Tmin_diff = abs(tr.Tmin_optim - minimum(env.tair_monthly))
  else
    Tmin_diff = 0
  end

  tair_mean_niche_change_stress =  exp( -pars.temp_niche_breadth_parameter * (Tmean_diff)^2 )
  tair_max_niche_change_stress =  exp( -pars.temp_niche_breadth_parameter * (Tmax_diff)^2 )
  tair_min_niche_change_stress =  exp( -pars.temp_niche_breadth_parameter * (Tmin_diff)^2 )
  tair_niche_change_stress = min(tair_mean_niche_change_stress, tair_max_niche_change_stress, tair_min_niche_change_stress)
  
  # phenology_status = 1 # In current model, always considered fully grown plants 
  
  FPAR = tr.FPC .* tair_niche_change_stress # phenology_status # Fraction PAR absorved by vegetation  
  APAR = FPAR .* PAR .* pars.alpha_leaf_to_stand # [mol E m-2 day-1] Effectively absorved photosynthetically active radiation APAR
  
  ##### Continue translation here. 
  # If plant is seasonal, derive favorable growing season adjusted for leaf longevity 
  if tr.seasonality == 0 # 0 represents a seasonal plant
    favorable_months = repeat([1.], 12)
    while sum(favorable_months) > tr.a_ll * 12
      if sum(favorable_months) - (tr.a_ll * 12) >= 1
        index = findall(env.tair_monthly .== minimum(env.tair_monthly[favorable_months .== 1]))
        favorable_months[index] .= 0.0
        # print("first case")
      elseif (sum(favorable_months)-(tr.a_ll * 12)) < 1
        index = findall(env.tair_monthly .== minimum(env.tair_monthly[favorable_months .== 1]))
        favorable_months[index] .= 1 - (sum(favorable_months) - (tr.a_ll * 12))
        # print("second case")
      end
    end
  end
  
  
  # Atmospheric CO2 partial pressure 
  p_a = (env.CO2_ppm * pars.p_atm)/1e+6 # Partial pressure of atmospheric CO2 [Pa]
  
  
  ############# Maximum potential photosynthesis only light and rubisco limited, no water limitation
  # Calculate Vm with maximum p_i/lambdamc3, not actual lambda (see Schaphoff)
  # Currently only low tairerature limitation, following Haxeltine & Prentice 1996, could also be trait dependent
  # f_tair_limit_low = 1/(1 + exp(0.2*(15-env.tair_monthly)))
  # f_tair_limit_low = 1/(1 + exp(0.35*(13.5-env.tair_monthly)))
  f_tair_limit_low = 1 ./ (1 .+ exp.( 0.25 .* (12 .- env.tair_monthly)))
  f_tair_limit_crit = 1 .- (1 ./ (1 .+ exp.(.-(env.tair_monthly .- 40.85))))
  f_tair_limit = f_tair_limit_low .* f_tair_limit_crit
  f_tair_limit[env.tair_monthly .<= pars.temp_crit_photosynthesis] .= 0
  
  p_i_max = pars.lambdamc3 * p_a
  
  Tau = 2600   .*  0.57.^( (env.tair_monthly  .- 25 )  ./  10 )
  
  Tau_star = (pars.pO2)  ./   (2   .*  Tau) # [Pa]
  
  c1 = pars.alpha_c3   .*  f_tair_limit   .*  ( ( p_i_max  .-  Tau_star )  ./   ( p_i_max  .+  ( 2  .*  Tau_star) ) )
  
  k_o = 30   .*  1000   .*  (1.2 .^ ((env.tair_monthly  .-  25)  ./   10)) # [Pa]
  
  k_c = 30   .*  (2.1 .^ ((env.tair_monthly .- 25) ./  10)) # [Pa]
  
  c2 = (p_i_max  .-  Tau_star) ./  (p_i_max  .+  (k_c   .*  (1  .+  (pars.pO2 ./  k_o))))
  
  # Calculate s
  s = (24 ./  env.daylength)   .*  pars.b_c3
  
  # Calculate sigma_c
  sigma = sqrt.( 1 .- ((c2 .- s) ./  (c2 .- pars.theta  .* s)) )
  
  # Calculate Vm (maximum rubisco capacity) [g C m .- 2 day .- 1]
  V_m = (1 ./  pars.b_c3)  .* (c1 ./  c2)  .* ((2  .* pars.theta  .-  1)  .* s  .-  (2  .* pars.theta  .* s  .-  c2)  .* sigma)  .* (APAR)  .* pars.c_mass
  
  # Calculate light .- limited photosynthesis on hourly basis [g C m .- 2 h .- 1]
  Je_hourly = c1   .*  (APAR ./  env.daylength)   .*  pars.c_mass
  
  # Calculate Rubisco .- limited photosynthesis rate [g C m .- 2 h .- 1]
  Jc_hourly = c2   .*  V_m   .*  (1 ./  24)
  
  # Calculate daytime gross photosynthesis Agd [g C m .- 2 daytime .- 1]
  Agd = ( Je_hourly  .+  Jc_hourly  .-  sqrt.((Je_hourly  .+  Jc_hourly) .^ 2  .-  4  .* pars.theta  .* Je_hourly  .* Jc_hourly) )  ./   (2   .*  pars.theta)   .*  env.daylength
  
  # Leaf respiration 
  R_leaf = pars.b_c3   .*  V_m
  
  # Net daily photosynthesis (minus leaf respiration)
  And = Agd  .-  R_leaf
  
  # Daily net daytime photosynthesis (gross photosynthesis  .-  daytime dark respiration)
  Adt = And  .+  ((1  .-  env.daylength ./  24)  .* R_leaf)
  #############

  ############# Water limitation
  # Consider water limitation  .- > actual lambda
  lambda = repeat([pars.lambdamc3], 12) # initialize lambda for later adjustment with bisection
  
  # Check water balance and calculate maximum photosynthesis possible under water stress
  Emax = 5 # [mm day .- 1] maximum evapotranspiration rate
  # Esupply = Emax   .*  (env.aridity_monthly) .^ ( 1  .-  (0.08   .*  tr.r_s_r))
  # Esupply[env.aridity_monthly <= 0.2] = Emin
  Esupply = min.(5, env.precip_monthly  .* 1000) # in mm ./  day
  Esupply = max.(0.1, Esupply) # 0.1  .* (tr.r_s_r ./  5)

  # # Calculate potential canopy conductance with unlimited photosynthesis rate
  gmin = 0.3 #0.3 # minimum water loss in [mm s .- 1]
  Adt_mol = Adt  ./   pars.c_mass # mol C m .- 2 day .- 1
  Adt_mol_s = Adt_mol  ./   (env.daylength   .*  60   .*  60) # mol C m .- 2 s .- 1
  c_a = (p_a ./  pars.p_atm) # ambient mole fraction of co2 (check Sitch 2003): p_a = p   .*  ca (p = atmospheric pressure)
  gc_pot_mol = ((1.6   .*  Adt_mol_s) ./  (c_a   .*  (1  .-  pars.lambdamc3))) # mol H2O m .- 2 s .- 1
  gc_pot_mm = gmin  .+  (gc_pot_mol   .*  8.314   .*  (env.tair_monthly  .+  273.15)  ./   pars.p_atm)   .*  1000 # m3 H2O m .- 2 s .- 1   .*  1000  .-  .- > mm H2O s .- 1 (using ideal gas constant)

  # # Slope of saturation vapour pressure curve [Pa K .- 1]
  sa = (2.502   .*  10 .^ 6)   .*  (exp.(17.269   .*  env.tair_monthly  ./   (237.3  .+  env.tair_monthly)) ./  (237.3  .+  env.tair_monthly) .^ 2)
  # latent heat vaporization [J kg .- 1]
  lat_heat_vap = (2.495   .*  10 .^ 6)  .-  2380   .*  env.tair_monthly
  # Daily surface net radiation [J m .- 2 day .- 1]
  Rn_day = (env.rss_monthly  .+  env.rls_monthly   .*  (env.daylength ./  24))   .*  (1e+6) # Consider only daytime longwave radiation # from MJ day .- 1 to J day .- 1 (daytime evapotranspiration)
  # Psychrometric constant [Pa K .- 1]
  psi_constant = 65.05  .+  (0.064   .*  env.tair_monthly) # with tair in °C?
  # Equilibrium evapotranspiration [mm day .- 1] or [kg H2O m .- 2 day .- 1]
  Eeq = (sa ./  (sa  .+  psi_constant))   .*  (Rn_day ./  lat_heat_vap)
  Eeq[Eeq .<= 0] .= 0
  Edemand = Eeq   .*  1.391 ./  (1  .+  (3.26 ./  gc_pot_mm))
  
  # ######### with chelsa data  .-  look up potential evapotranspiration
  # Eeq = env.pet_monthly    .*  1000 # in mm ./  day
  # Eeq[Eeq <= 0] = 0
  # Edemand = Eeq #    .*  1.391 ./  (1  .+  (3.26 ./  gc_pot_mm))

  # Are there months with water limitation?
  water_stressed = findall(Esupply .< Edemand)
  
  if length(water_stressed) > 0
    # print("water limitation!")
    # browser()
    for i in water_stressed

      find_lambda = function(lambda)
          
        # Calculate Adt with photosynthesis module
        p_i = lambda   .*  p_a # [Pa]
        c1 = pars.alpha_c3   .*  f_tair_limit[i]   .*  ( ( p_i  .-  Tau_star[i] )  ./   ( p_i  .+  ( 2  .* Tau_star[i]) ) )
        c2 = (p_i  .-  Tau_star[i]) ./  (p_i  .+  (k_c[i]   .*  (1  .+  (pars.pO2 ./  k_o[i]))))
          
        # Calculate light .- limited photosynthesis on hourly basis [g C m .- 2 h .- 1]
        Je_hourly = c1   .*  (APAR[i] ./  env.daylength[i])   .*  pars.c_mass
          
        # Calculate Rubisco .- limited photosynthesis rate [g C m .- 2 h .- 1]
        Jc_hourly = c2   .*  V_m[i]   .*  (1 ./  24)
          
        # Calculate daytime gross photosynthesis Agd [g C m .- 2 daytime .- 1]
        Agd = ( Je_hourly  .+  Jc_hourly  .-  sqrt((Je_hourly  .+  Jc_hourly) .^ 2  .-  4  .* pars.theta  .* Je_hourly  .* Jc_hourly) )  ./   (2   .*  pars.theta)   .*  env.daylength[i]
        # Leaf respiration 
        R_leaf = pars.b_c3   .*  V_m[i]
        # Net daily photosynthesis (minus leaf respiration)
        And = Agd  .-  R_leaf
        # Daily net daytime photosynthesis (gross photosynthesis  .-  daytime dark respiration)
        Adt_v1 = And  .+  ((1  .-  env.daylength[i] ./  24)  .* R_leaf)
          
        # Calculate canopy conductance that meets water supply limitation
        gc = 3.26 ./  (((Eeq[i]   .*  1.391) ./  Esupply[i]) .- 1)
        if gc < 0 
          gc = 0
        end

        Adt_v2_mm_s = c_a   .*  (gc  .-  gmin)   .*  ((1  .-  lambda) ./  1.6) # in mm s .- 1 
        Adt_v2_m_s = Adt_v2_mm_s   .*  (1 ./  1000) # m s .- 1  .-  .- > m3 m .- 2 s .- 1
        Adt_v2_mol_day = Adt_v2_m_s   .*  (pars.p_atm  ./   (8.314   .*  (env.tair_monthly[i]  .+  273.15)))   .*  (env.daylength[i]   .*  60   .*  60) # convert volume co2 gas to mol m .- 2 day .- 1
        Adt_v2 = Adt_v2_mol_day   .*  pars.c_mass # convert to g C m .- 2 day .- 1
        Adt_v2 = max.(- 0.5, Adt_v2) # Ensure root in a negative number to avoid failure of bisection algorithm
          
        return(Adt_v1  .-  Adt_v2)
          
      end
        
      # Check for opposite signs: 
      sign_1 = sign(find_lambda(0))
      sign_2 = sign(find_lambda(0.8))
      if sign_1 == sign_2 # In case there is no intersect  .-  no prodcutivity, set lambda to zero
        lambda[i] = 0
      else
        lambda[i] = find_zero(find_lambda, (0, 0.8), Bisection(), maxiters = 100)
      end

    end

  end
    
  # Recalculate Photosynthesis with updated ./  actual lambda considerin water limitation
  p_i = lambda   .*  p_a # [Pa]
  c1 = pars.alpha_c3   .*  f_tair_limit   .*  ( ( p_i  .-  Tau_star )  ./   ( p_i  .+  ( 2  .* Tau_star) ) )
  c2 = (p_i  .-  Tau_star) ./  (p_i  .+  (k_c   .*  (1  .+  (pars.pO2 ./  k_o))))
  
  # Calculate light .- limited photosynthesis on hourly basis [g C m .- 2 h .- 1]
  Je_hourly = c1   .*  (APAR ./  env.daylength)   .*  pars.c_mass
  
  # Calculate Rubisco .- limited photosynthesis rate [g C m .- 2 h .- 1]
  Jc_hourly = c2   .*  V_m   .*  (1 ./  24)
  
  # Calculate daytime gross photosynthesis Agd [g C m .- 2 daytime .- 1]
  Agd = ( Je_hourly  .+  Jc_hourly  .-  sqrt.((Je_hourly  .+  Jc_hourly) .^ 2  .-  4  .* pars.theta  .* Je_hourly  .* Jc_hourly) )  ./   (2   .*  pars.theta)   .*  env.daylength
  
  
  ################ Calculate AET  .-  at the moment just zero during non .- growing season
  # Leaf respiration
  R_leaf = pars.b_c3   .*  V_m
  # Net daily photosynthesis (minus leaf respiration)
  And = Agd  .-  R_leaf
  # Daily net daytime photosynthesis (gross photosynthesis  .-  daytime dark respiration)
  Adt = And  .+  ((1  .-  env.daylength ./  24)  .* R_leaf)

  # Calculate canopy conductance
  gmin = 0.3 #0.3 # minimum water loss in [mm s .- 1]
  Adt_mol = Adt  ./   pars.c_mass # mol C m .- 2 day .- 1
  Adt_mol_s = Adt_mol  ./   (env.daylength   .*  60   .*  60) # mol C m .- 2 s .- 1
  c_a = (p_a ./  pars.p_atm) # ambient mole fraction of co2 (check Sitch 2003): p_a = p   .*  ca (p = atmospheric pressure)
  gc_act_mol = ((1.6   .*  Adt_mol_s) ./  (c_a   .*  (1  .-  lambda))) # mol H2O m .- 2 s .- 1
  gc_act_mm = gmin  .+  (gc_act_mol   .*  8.314   .*  (env.tair_monthly  .+  273.15)  ./   pars.p_atm)   .*  1000 # m3 H2O m .- 2 s .- 1   .*  1000  .-  .- > mm H2O s .- 1 (using ideal gas constant)
  
  # # Slope of saturation vapour pressure curve [Pa K .- 1]
  sa = (2.502   .*  10 .^ 6)   .*  (exp.(17.269   .*  env.tair_monthly  ./   (237.3  .+  env.tair_monthly)) ./  (237.3  .+  env.tair_monthly) .^ 2)
  # latent heat vaporization [J kg .- 1]
  lat_heat_vap = (2.495   .*  10 .^ 6)  .-  2380   .*  env.tair_monthly
  # Daily surface net radiation [J m .- 2 day .- 1]
  Rn_day = (env.rss_monthly  .+  env.rls_monthly   .*  (env.daylength ./  24))   .*  (1e+6) # Consider only daytime longwave radiation # from MJ day .- 1 to J day .- 1 (daytime evapotranspiration)
  # Psychrometric constant [Pa K .- 1]
  psi_constant = 65.05  .+  (0.064   .*  env.tair_monthly) # with tair in °C?
  # Equilibrium evapotranspiration [mm day .- 1] or [kg H2O m .- 2 day .- 1]
  Eeq = (sa ./  (sa  .+  psi_constant))   .*  (Rn_day ./  lat_heat_vap)
  Eeq[Eeq .<= 0] .= 0
  Edemand_act = Eeq   .*  1.391 ./  (1  .+  (3.26 ./  gc_act_mm))
  
  AET_monthly = min.(Esupply, Edemand_act)

  # Integration
  Agd = max.(0, Agd) # Take care of rare cases with negative daily photosynthesis (e.g., lambda == 0)
  Agd = Agd    .*  pars.month_days 
  GPP = sum(Agd)
  AET = AET_monthly   .*  pars.month_days
  AET_sum = sum(AET)
  V_m = V_m 
  
  if tr.seasonality == 0
    Agd = Agd   .*  favorable_months
    GPP = sum(Agd)
    AET = AET   .*  favorable_months
    AET_sum = sum(AET)
    V_m = V_m   .*  favorable_months
  end
  
  #######
  # ### Test Results with different approach with first deriving And then adding nighttime respiration
  # And = APAR[favorable_months]   .*  (c1 ./  c2)   .*  (c2  .-  (2  .* pars.theta  .- 1)  .* s  .-  2  .* (c2 .- pars.theta  .* s)  .* sigma)   .*  pars.c_mass
  # Agd = And  .+  (pars.b_c3   .*  V_m)
  ########
  
  
  GPP_out = (
  GPP = GPP,
  V_m = V_m,
  AET = AET_sum
  )
  
  return(GPP_out)

end

### GPP function for optimizaiton - avoid bisection for lambda & AET calculation
GPP_function_for_optimization = function(; env, tr, pars)
  
  daytime_integrated_rsds = env.rsds_monthly # * (24/env.daylength) # [MJ m-2 daytime-1/day-1]
  PAR = daytime_integrated_rsds * 4.6 * 0.5 # Photosynthetically active radiation, Schlaphoff version 
  
  # Climatic niche change stress
  Tmean_diff = abs(env.tair_annual - tr.Tave_optim)
  if maximum(env.tair_monthly) > tr.Tmax_optim
    Tmax_diff = abs(maximum(env.tair_monthly) - tr.Tmax_optim)
  else
    Tmax_diff = 0
  end

  if minimum(env.tair_monthly) < tr.Tmin_optim
    Tmin_diff = abs(tr.Tmin_optim - minimum(env.tair_monthly))
  else
    Tmin_diff = 0
  end

  tair_mean_niche_change_stress =  exp( -pars.temp_niche_breadth_parameter * (Tmean_diff)^2 )
  tair_max_niche_change_stress =  exp( -pars.temp_niche_breadth_parameter * (Tmax_diff)^2 )
  tair_min_niche_change_stress =  exp( -pars.temp_niche_breadth_parameter * (Tmin_diff)^2 )
  tair_niche_change_stress = min(tair_mean_niche_change_stress, tair_max_niche_change_stress, tair_min_niche_change_stress)
  
  # phenology_status = 1 # In current model, always considered fully grown plants 
  
  FPAR = tr.FPC .* tair_niche_change_stress # phenology_status # Fraction PAR absorved by vegetation  
  APAR = FPAR .* PAR .* pars.alpha_leaf_to_stand # [mol E m-2 day-1] Effectively absorved photosynthetically active radiation APAR
  
  ##### Continue translation here. 
  # If plant is seasonal, derive favorable growing season adjusted for leaf longevity 
  if tr.seasonality == 0 # 0 represents a seasonal plant
    favorable_months = repeat([1.], 12)
    while sum(favorable_months) > tr.a_ll * 12
      if sum(favorable_months) - (tr.a_ll * 12) >= 1
        index = findall(env.tair_monthly .== minimum(env.tair_monthly[favorable_months .== 1]))
        favorable_months[index] .= 0.0
        # print("first case")
      elseif (sum(favorable_months)-(tr.a_ll * 12)) < 1
        index = findall(env.tair_monthly .== minimum(env.tair_monthly[favorable_months .== 1]))
        favorable_months[index] .= 1 - (sum(favorable_months) - (tr.a_ll * 12))
        # print("second case")
      end
    end
  end
  
  
  # Atmospheric CO2 partial pressure 
  p_a = (env.CO2_ppm * pars.p_atm)/1e+6 # Partial pressure of atmospheric CO2 [Pa]
  
  
  ############# Maximum potential photosynthesis only light and rubisco limited, no water limitation
  # Calculate Vm with maximum p_i/lambdamc3, not actual lambda (see Schaphoff)
  # Currently only low tairerature limitation, following Haxeltine & Prentice 1996, could also be trait dependent
  # f_tair_limit_low = 1/(1 + exp(0.2*(15-env.tair_monthly)))
  # f_tair_limit_low = 1/(1 + exp(0.35*(13.5-env.tair_monthly)))
  f_tair_limit_low = 1 ./ (1 .+ exp.( 0.25 .* (12 .- env.tair_monthly)))
  f_tair_limit_crit = 1 .- (1 ./ (1 .+ exp.(.-(env.tair_monthly .- 40.85))))
  f_tair_limit = f_tair_limit_low .* f_tair_limit_crit
  f_tair_limit[env.tair_monthly .<= pars.temp_crit_photosynthesis] .= 0
  
  p_i_max = pars.lambdamc3 * p_a
  
  Tau = 2600   .*  0.57.^( (env.tair_monthly  .- 25 )  ./  10 )
  
  Tau_star = (pars.pO2)  ./   (2   .*  Tau) # [Pa]
  
  c1 = pars.alpha_c3   .*  f_tair_limit   .*  ( ( p_i_max  .-  Tau_star )  ./   ( p_i_max  .+  ( 2  .*  Tau_star) ) )
  
  k_o = 30   .*  1000   .*  (1.2 .^ ((env.tair_monthly  .-  25)  ./   10)) # [Pa]
  
  k_c = 30   .*  (2.1 .^ ((env.tair_monthly .- 25) ./  10)) # [Pa]
  
  c2 = (p_i_max  .-  Tau_star) ./  (p_i_max  .+  (k_c   .*  (1  .+  (pars.pO2 ./  k_o))))
  
  # Calculate s
  s = (24 ./  env.daylength)   .*  pars.b_c3
  
  # Calculate sigma_c
  sigma_argument = 1 .- ((c2 .- s) ./  (c2 .- pars.theta  .* s))
  sigma_argument[sigma_argument .< 0] .= 0 
  sigma = sqrt.( sigma_argument )
  
  # Calculate Vm (maximum rubisco capacity) [g C m .- 2 day .- 1]
  V_m = (1 ./  pars.b_c3)  .* (c1 ./  c2)  .* ((2  .* pars.theta  .-  1)  .* s  .-  (2  .* pars.theta  .* s  .-  c2)  .* sigma)  .* (APAR)  .* pars.c_mass
  
  # Calculate light .- limited photosynthesis on hourly basis [g C m .- 2 h .- 1]
  Je_hourly = c1   .*  (APAR ./  env.daylength)   .*  pars.c_mass
  
  # Calculate Rubisco .- limited photosynthesis rate [g C m .- 2 h .- 1]
  Jc_hourly = c2   .*  V_m   .*  (1 ./  24)
  
  # Calculate daytime gross photosynthesis Agd [g C m .- 2 daytime .- 1]
  Agd = ( Je_hourly  .+  Jc_hourly  .-  sqrt.((Je_hourly  .+  Jc_hourly) .^ 2  .-  4  .* pars.theta  .* Je_hourly  .* Jc_hourly) )  ./   (2   .*  pars.theta)   .*  env.daylength
  
  # Leaf respiration 
  R_leaf = pars.b_c3   .*  V_m
  
  # Net daily photosynthesis (minus leaf respiration)
  And = Agd  .-  R_leaf
  
  # Daily net daytime photosynthesis (gross photosynthesis  .-  daytime dark respiration)
  Adt = And  .+  ((1  .-  env.daylength ./  24)  .* R_leaf)
  #############

  ############# Water limitation
  # Consider water limitation  .- > actual lambda
  lambda = repeat([pars.lambdamc3], 12) # initialize lambda for later adjustment with bisection
  
  # Check water balance and calculate maximum photosynthesis possible under water stress
  Emax = 5 # [mm day .- 1] maximum evapotranspiration rate
  # Esupply = Emax   .*  (env.aridity_monthly) .^ ( 1  .-  (0.08   .*  tr.r_s_r))
  # Esupply[env.aridity_monthly <= 0.2] = Emin
  Esupply = min.(Emax, env.precip_monthly  .* 1000) # in mm ./  day
  #Esupply = max.(0.1, Esupply) # 0.1  .* (tr.r_s_r ./  5)

  # # Calculate potential canopy conductance with unlimited photosynthesis rate
  gmin = 0.5 #0.3 # minimum water loss in [mm s .- 1]
  Adt_mol = Adt  ./   pars.c_mass # mol C m .- 2 day .- 1
  Adt_mol_s = Adt_mol  ./   (env.daylength   .*  60   .*  60) # mol C m .- 2 s .- 1
  c_a = (p_a ./  pars.p_atm) # ambient mole fraction of co2 (check Sitch 2003): p_a = p   .*  ca (p = atmospheric pressure)
  gc_pot_mol = ((1.6   .*  Adt_mol_s) ./  (c_a   .*  (1  .-  pars.lambdamc3))) # mol H2O m .- 2 s .- 1
  gc_pot_mm = gmin  .+  (gc_pot_mol   .*  8.314   .*  (env.tair_monthly  .+  273.15)  ./   pars.p_atm)   .*  1000 # m3 H2O m .- 2 s .- 1   .*  1000  .-  .- > mm H2O s .- 1 (using ideal gas constant)

  # # Slope of saturation vapour pressure curve [Pa K .- 1]
  sa = (2.502   .*  10 .^ 6)   .*  (exp.(17.269   .*  env.tair_monthly  ./   (237.3  .+  env.tair_monthly)) ./  (237.3  .+  env.tair_monthly) .^ 2)
  # latent heat vaporization [J kg .- 1]
  lat_heat_vap = (2.495   .*  10 .^ 6)  .-  2380   .*  env.tair_monthly
  # Daily surface net radiation [J m .- 2 day .- 1]
  Rn_day = (env.rss_monthly  .+  env.rls_monthly   .*  (env.daylength ./  24))   .*  (1e+6) # Consider only daytime longwave radiation # from MJ day .- 1 to J day .- 1 (daytime evapotranspiration)
  # Psychrometric constant [Pa K .- 1]
  psi_constant = 65.05  .+  (0.064   .*  env.tair_monthly) # with tair in °C?
  # Equilibrium evapotranspiration [mm day .- 1] or [kg H2O m .- 2 day .- 1]
  Eeq = (sa ./  (sa  .+  psi_constant))   .*  (Rn_day ./  lat_heat_vap)
  Eeq[Eeq .<= 0] .= 0
  Edemand = Eeq   .*  1.391 ./  (1  .+  (3.26 ./  gc_pot_mm))
  
  # ######### with chelsa data  .-  look up potential evapotranspiration
  # Eeq = env.pet_monthly    .*  1000 # in mm ./  day
  # Eeq[Eeq <= 0] = 0
  # Edemand = Eeq #    .*  1.391 ./  (1  .+  (3.26 ./  gc_pot_mm))

  # Are there months with water limitation?
  water_stressed = findall(Esupply .< Edemand)
  
  if length(water_stressed) > 0
    # print("water limitation!")
    # browser()
    for i in water_stressed
          performance = Vector{Float64}() 
          test_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
          for w in test_values
          
              # Calculate Adt with photosynthesis module
              p_i = w   .*  p_a # [Pa]
              c1 = pars.alpha_c3   .*  f_tair_limit[i]   .*  ( ( p_i  .-  Tau_star[i] )  ./   ( p_i  .+  ( 2  .* Tau_star[i]) ) )
              c2 = (p_i  .-  Tau_star[i]) ./  (p_i  .+  (k_c[i]   .*  (1  .+  (pars.pO2 ./  k_o[i]))))
              
              # Calculate light .- limited photosynthesis on hourly basis [g C m .- 2 h .- 1]
              Je_hourly = c1   .*  (APAR[i] ./  env.daylength[i])   .*  pars.c_mass
              
              # Calculate Rubisco .- limited photosynthesis rate [g C m .- 2 h .- 1]
              Jc_hourly = c2   .*  V_m[i]   .*  (1 ./  24)
              
              # Calculate daytime gross photosynthesis Agd [g C m .- 2 daytime .- 1]
              Agd = ( Je_hourly  .+  Jc_hourly  .-  sqrt((Je_hourly  .+  Jc_hourly) .^ 2  .-  4  .* pars.theta  .* Je_hourly  .* Jc_hourly) )  ./   (2   .*  pars.theta)   .*  env.daylength[i]
              # Leaf respiration 
              R_leaf = pars.b_c3   .*  V_m[i]
              # Net daily photosynthesis (minus leaf respiration)
              And = Agd  .-  R_leaf
              # Daily net daytime photosynthesis (gross photosynthesis  .-  daytime dark respiration)
              Adt_v1 = And  .+  ((1  .-  env.daylength[i] ./  24)  .* R_leaf)
              
              # Calculate canopy conductance that meets water supply limitation
              gc = 3.26 ./  (((Eeq[i]   .*  1.391) ./  Esupply[i]) .- 1)
              if gc < 0 
                  gc = 0
              end
      
              Adt_v2_mm_s = c_a   .*  (gc  .-  gmin)   .*  ((1  .-  w) ./  1.6) # in mm s .- 1 
              Adt_v2_m_s = Adt_v2_mm_s   .*  (1 ./  1000) # m s .- 1  .-  .- > m3 m .- 2 s .- 1
              Adt_v2_mol_day = Adt_v2_m_s   .*  (pars.p_atm  ./   (8.314   .*  (env.tair_monthly[i]  .+  273.15)))   .*  (env.daylength[i]   .*  60   .*  60) # convert volume co2 gas to mol m .- 2 day .- 1
              Adt_v2 = Adt_v2_mol_day   .*  pars.c_mass # convert to g C m .- 2 day .- 1
              Adt_v2 = max.(- 0.5, Adt_v2) # Ensure root in a negative number to avoid failure of bisection algorithm
              
              push!(performance, Adt_v1  .-  Adt_v2)
          end
          
          # Check for opposite signs: 
          sign_1 = sign(performance[1])
          sign_2 = sign(performance[length(performance)])
          if sign_1 == sign_2 # In case there is no intersect  .-  no prodcutivity, set lambda to zero
              lambda[i] = 0
          else
              lambda[i] = test_values[argmin(abs.(performance))]
          end
    end

  end
    
  # Recalculate Photosynthesis with updated ./  actual lambda considerin water limitation
  p_i = lambda   .*  p_a # [Pa]
  c1 = pars.alpha_c3   .*  f_tair_limit   .*  ( ( p_i  .-  Tau_star )  ./   ( p_i  .+  ( 2  .* Tau_star) ) )
  c2 = (p_i  .-  Tau_star) ./  (p_i  .+  (k_c   .*  (1  .+  (pars.pO2 ./  k_o))))
  
  # Calculate light .- limited photosynthesis on hourly basis [g C m .- 2 h .- 1]
  Je_hourly = c1   .*  (APAR ./  env.daylength)   .*  pars.c_mass
  
  # Calculate Rubisco .- limited photosynthesis rate [g C m .- 2 h .- 1]
  Jc_hourly = c2   .*  V_m   .*  (1 ./  24)
  
  # Calculate daytime gross photosynthesis Agd [g C m .- 2 daytime .- 1]
  Agd = ( Je_hourly  .+  Jc_hourly  .-  sqrt.((Je_hourly  .+  Jc_hourly) .^ 2  .-  4  .* pars.theta  .* Je_hourly  .* Jc_hourly) )  ./   (2   .*  pars.theta)   .*  env.daylength
  

  ################ Calculate AET  .-  at the moment just zero during non .- growing season
  # Leaf respiration
  R_leaf = pars.b_c3   .*  V_m
  # Net daily photosynthesis (minus leaf respiration)
  And = Agd  .-  R_leaf
  # Daily net daytime photosynthesis (gross photosynthesis  .-  daytime dark respiration)
  Adt = And  .+  ((1  .-  env.daylength ./  24)  .* R_leaf)

  # Calculate canopy conductance
  gmin = 0.3 #0.3 # minimum water loss in [mm s .- 1]
  Adt_mol = Adt  ./   pars.c_mass # mol C m .- 2 day .- 1
  Adt_mol_s = Adt_mol  ./   (env.daylength   .*  60   .*  60) # mol C m .- 2 s .- 1
  c_a = (p_a ./  pars.p_atm) # ambient mole fraction of co2 (check Sitch 2003): p_a = p   .*  ca (p = atmospheric pressure)
  gc_act_mol = ((1.6   .*  Adt_mol_s) ./  (c_a   .*  (1  .-  lambda))) # mol H2O m .- 2 s .- 1
  gc_act_mm = gmin  .+  (gc_act_mol   .*  8.314   .*  (env.tair_monthly  .+  273.15)  ./   pars.p_atm)   .*  1000 # m3 H2O m .- 2 s .- 1   .*  1000  .-  .- > mm H2O s .- 1 (using ideal gas constant)
  
  # # Slope of saturation vapour pressure curve [Pa K .- 1]
  sa = (2.502   .*  10 .^ 6)   .*  (exp.(17.269   .*  env.tair_monthly  ./   (237.3  .+  env.tair_monthly)) ./  (237.3  .+  env.tair_monthly) .^ 2)
  # latent heat vaporization [J kg .- 1]
  lat_heat_vap = (2.495   .*  10 .^ 6)  .-  2380   .*  env.tair_monthly
  # Daily surface net radiation [J m .- 2 day .- 1]
  Rn_day = (env.rss_monthly  .+  env.rls_monthly   .*  (env.daylength ./  24))   .*  (1e+6) # Consider only daytime longwave radiation # from MJ day .- 1 to J day .- 1 (daytime evapotranspiration)
  # Psychrometric constant [Pa K .- 1]
  psi_constant = 65.05  .+  (0.064   .*  env.tair_monthly) # with tair in °C?
  # Equilibrium evapotranspiration [mm day .- 1] or [kg H2O m .- 2 day .- 1]
  Eeq = (sa ./  (sa  .+  psi_constant))   .*  (Rn_day ./  lat_heat_vap)
  Eeq[Eeq .<= 0] .= 0
  Edemand_act = Eeq   .*  1.391 ./  (1  .+  (3.26 ./  gc_act_mm))
  
  AET_monthly = min.(Esupply, Edemand_act)

  # Integration
  Agd = max.(0, Agd) # Take care of rare cases with negative daily photosynthesis (e.g., lambda == 0)
  Agd = Agd    .*  pars.month_days 
  GPP = sum(Agd)
  AET = AET_monthly   .*  pars.month_days
  AET_sum = sum(AET)
  V_m = V_m 
  
  if tr.seasonality == 0
    Agd = Agd   .*  favorable_months
    GPP = sum(Agd)
    AET = AET   .*  favorable_months
    AET_sum = sum(AET)
    V_m = V_m   .*  favorable_months
  end
  
  #######
  # ### Test Results with different approach with first deriving And then adding nighttime respiration
  # And = APAR[favorable_months]   .*  (c1 ./  c2)   .*  (c2  .-  (2  .* pars.theta  .- 1)  .* s  .-  2  .* (c2 .- pars.theta  .* s)  .* sigma)   .*  pars.c_mass
  # Agd = And  .+  (pars.b_c3   .*  V_m)
  ########
  
  
  GPP_out = (
  GPP = GPP,
  V_m = V_m,
  AET = AET_sum
  )
  
  return(GPP_out)

end

### Maintenance respiration 
R_maintenance_function =  function(; env, tr, pars, GPP_out)
  
  # If plant is seasonal, derive favorable growing season adjusted for leaf longevity 
  if tr.seasonality == 0 # 0 represents a seasonal plant
    favorable_months = repeat([1.], 12)
    while sum(favorable_months) > tr.a_ll * 12
      if sum(favorable_months) - (tr.a_ll * 12) >= 1
        index = findall(env.tair_monthly .== minimum(env.tair_monthly[favorable_months .== 1]))
        favorable_months[index] .= 0.0
        # print("first case")
      elseif (sum(favorable_months)-(tr.a_ll * 12)) < 1
        index = findall(env.tair_monthly .== minimum(env.tair_monthly[favorable_months .== 1]))
        favorable_months[index] .= 1 - (sum(favorable_months) - (tr.a_ll * 12))
        # print("second case")
      end
    end
  end
  
  # Derive maximum carboxylation rates to calculate respiration (as in LPJ)
  V_m =  GPP_out.V_m
  
  g_T =  exp.(308.56  .*  ((1 ./ 56.02)  .-  (1 ./ (env.tair_monthly  .+  46.02))))
  
  # Simplified soil temperature  .-  just damped trajectory
  temp_soil =  mean(env.tair_monthly)  .+  (env.tair_monthly  .-  mean(env.tair_monthly)) ./ 2
  g_T_soil =  exp.(308.56  .*  ((1 ./ 56.02)  .-  (1 ./ (temp_soil  .+  46.02))))
  
  g_T[env.tair_monthly .<=  .- 20] .=  0
  g_T_soil[env.tair_monthly .<=  .- 20] .=  0
  
  # env dependent temperature acclimation of respiration  .-  linear acclimation factor (max. 0.5 reduction)
  acclimation_factor =  1.35  .-  0.035  .*  mean(env.tair_monthly)
  if acclimation_factor < 0.3
    acclimation_factor =  0.3
  elseif acclimation_factor > 1
    acclimation_factor =  1
  end
  
  # Active carbon pools
  r_leaf =  sum(pars.b_c3  .*  V_m  .*  pars.month_days) # Is already adjusted for growing season, if plant is seasonal (see V_m)
  
  r_sapwood =  sum((pars.r   .*  pars.k  .*  acclimation_factor  .*  (tr.C_sapwood ./ pars.cn_wood)  .*  g_T)  .*  pars.month_days)  ./  tr.CA
  r_fineroot =  sum((pars.r  .*  pars.k  .*  acclimation_factor  .*  (tr.C_fineroot ./ pars.cn_root)  .*  g_T_soil)  .*  pars.month_days)  ./  tr.CA
  
  # Structural carbon pools  .-  no ./ reduced respiration
  r_heartwood =  sum((pars.r   .*  pars.k  .*  (0 ./ 15) .*  acclimation_factor  .*  (tr.C_heartwood ./ pars.cn_wood)  .*  g_T)  .*  pars.month_days)  ./  tr.CA
  r_coarseroot =  sum((pars.r  .*  pars.k .*  (0 ./ 15) .*  acclimation_factor  .*  (tr.C_coarseroot ./ pars.cn_root)  .*  g_T_soil)  .*  pars.month_days)  ./  tr.CA
  

  if tr.seasonality == 0
    r_sapwood =  sum((pars.r   .*  pars.k  .*  acclimation_factor  .*  (tr.C_sapwood ./ pars.cn_wood)  .*  g_T)  .*  pars.month_days  .*  favorable_months)  ./  tr.CA
    r_fineroot =  sum((pars.r  .*  pars.k  .*  acclimation_factor  .*  (tr.C_fineroot ./ pars.cn_root)  .*  g_T_soil)  .*  pars.month_days  .*  favorable_months)  ./  tr.CA

    # Structural carbon pools  .-  no ./ reduced respiration
    r_heartwood =  sum((pars.r   .*  pars.k  .*  (0 ./ 15) .*  acclimation_factor  .*  (tr.C_heartwood ./ pars.cn_wood)  .*  g_T)  .*  pars.month_days  .*  favorable_months)  ./  tr.CA
    r_coarseroot =  sum((pars.r  .*  pars.k .*  (0 ./ 15) .*  acclimation_factor  .*  (tr.C_coarseroot ./ pars.cn_root)  .*  g_T_soil)  .*  pars.month_days  .*  favorable_months)  ./  tr.CA
  end
  
  R_maintenance_total =  r_leaf  .+  r_sapwood  .+  r_fineroot  .+  r_heartwood  .+  r_coarseroot
  
  return(R_maintenance_total)

end

### Tissue turnover 
C_turnover_function =  function(; env, tr, pars)

  # Heat and frost mortaility - frost mortatility reduced for evergreens compared to deciduous
  # Increase turnover in mortality-prone environments
  # if tr.seasonality == 1
  #   heat_damage = 2 * ( min(sum(env.tair_monthly .> 50), 5) / 5 ) 
  #   frost_damage = 2 * (  min(sum(env.tair_monthly .< -20), 5) / 5 ) 
  # elseif  tr.seasonality == 0
  #   heat_damage = 2 * ( min(sum(env.tair_monthly .> 50), 5) / 5 )  
  #   frost_damage = 2 * (  min(sum(env.tair_monthly .< -10), 5) / 5 )  # stronger effect of low temperatures on deciduous plants
  # end
  heat_damage = 0
  frost_damage = 0
  
  # Calculate C turnover 
  # Depends on annual or perennial leaf cycle
  if tr.seasonality == 1 # perennial plants

    if tr.a_ll >= 1
      C_turnover_fineroot =  tr.C_fineroot   .*  (1 ./ tr.a_ll)
      C_turnover_coarseroot =  tr.C_coarseroot  .*  pars.f_coarseroot  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_sapwood =  tr.C_sapwood  .*  pars.f_sapwood  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_heartwood =  tr.C_heartwood  .*  pars.f_heartwood  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_leaves =  tr.C_leaf  .*  (1 ./ tr.a_ll)
      
      C_turnover_total =  C_turnover_fineroot  .+  C_turnover_coarseroot  .+  C_turnover_sapwood  .+  C_turnover_heartwood  .+  C_turnover_leaves

    elseif tr.a_ll < 1
      
      C_turnover_fineroot =  tr.C_fineroot  .*  1  #(2  .-  tr.a_ll) # Increased root costs for more acquisitive strategies is unknown
      C_turnover_coarseroot =  tr.C_coarseroot  .*  pars.f_coarseroot  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_sapwood =  tr.C_sapwood  .*  pars.f_sapwood  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_heartwood =  tr.C_heartwood  .*  pars.f_heartwood  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_leaves =  tr.C_leaf  .*  (1 ./ tr.a_ll)
      
      C_turnover_total =  C_turnover_fineroot  .+  C_turnover_coarseroot  .+  C_turnover_sapwood  .+  C_turnover_heartwood  .+  C_turnover_leaves
    end

  elseif tr.seasonality == 0 # seasonal plants
    
    if tr.a_ll >= 1 
      C_turnover_fineroot =  tr.C_fineroot   .*  (1 ./ tr.a_ll)
      C_turnover_coarseroot =  tr.C_coarseroot  .*  pars.f_coarseroot  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_sapwood =  tr.C_sapwood  .*  pars.f_sapwood  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_heartwood =  tr.C_heartwood  .*  pars.f_heartwood  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_leaves =  tr.C_leaf  .*  (1 ./ tr.a_ll)
      
      C_turnover_total =  C_turnover_fineroot  .+  C_turnover_coarseroot  .+  C_turnover_sapwood  .+  C_turnover_heartwood  .+  C_turnover_leaves
      
    elseif tr.a_ll < 1
      
      C_turnover_fineroot =  tr.C_fineroot  .*  1 
      C_turnover_coarseroot =  tr.C_coarseroot  .*  pars.f_coarseroot  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_sapwood =  tr.C_sapwood  .*  pars.f_sapwood  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_heartwood =  tr.C_heartwood  .*  pars.f_heartwood  .*  (1  .+  heat_damage  .+  frost_damage)
      C_turnover_leaves =  tr.C_leaf  .*  1
      
      C_turnover_total =  C_turnover_fineroot  .+  C_turnover_coarseroot  .+  C_turnover_sapwood  .+  C_turnover_heartwood  .+  C_turnover_leaves
    end
    
  end
  
  # Area normalization
  C_turnover_total =  C_turnover_total ./ tr.CA

  return(C_turnover_total)

end

### Calculate NPP 
calc_NPP =  function(; GPP_out, R_maintenance, pars)
  NPP =  (1  .-  pars.r_gr)  .*  (GPP_out.GPP  .-  R_maintenance)
  return(NPP)
end

### Calculate carbon balance - net C gain after respiration maintenance/growth, reproduction carbon investments and tissue turnover 
calc_net_C_gain =  function(; NPP, C_turnover_total)
  net_C_gain =  ((1 - pars.r_repr) .* NPP)  .-  C_turnover_total # area normalized (CA) - > [g C m-2 yr-1]
  return(net_C_gain)
end

### Objective function to be minimized
trait_optimization_function =  function(x, pass_params)

  # Unpack passed arguments 
  env = pass_params.env
  # tr = pass_params.tr # New optimization, complete decoupling of optimization from current trait values
  pars = pass_params.pars

  # If seasonality is an option (based on env)  --> enter seasonal optimization 
  # This is done by a preliminary assignment of seasonal habit. Will be changed in trait evolution depending on actual a_ll
  # old: if (sum(env.tair_monthly .<= pars.temp_threshold_growing_season) >= 3. && tr.a_ll >= 1.) || (tr.seasonality == 0. && sum(env.tair_monthly .<= pars.temp_threshold_growing_season) >= 3.)
  # new: The decoupling between optimization and current traits might give some non-optimal behavior for the seasonal locations under slow evolution: evergreen might be more optimal but plant will tend to evolve towards seasonal, resulting in a transient non-optimal state
  if sum(env.tair_monthly .<= pars.temp_threshold_growing_season) >= 3.
    seasonality_optim =  0
  else
    seasonality_optim =  1
  end

  # Define optimization trait variables 
  tr_optim = (

  # Optimization variables
  a_ll =  x[1], 
  C_leaf =  x[2],
  H = x[3],

  #r_s_r_optim =  par_vec[3]
  r_s_r =  1,  # Assume constant r_s_r
  
  # Derive optimum fineroot C based on water balance optimization 
  C_fineroot =  x[2]  *  1, # Assume constant r_s_r
  seasonality =   seasonality_optim,

  # Transfer current climatic niche to optimization function
  Tave_optim =  mean(env.tair_monthly),
  Tmax_optim =  maximum(env.tair_monthly),
  Tmin_optim =  minimum(env.tair_monthly),
  Pave_optim =  mean(env.precip_monthly)
  )
  
  tr_optim =  plant_allometry(tr = tr_optim,  pars = pars) # Hypothetical plant
  GPP_out_optim =  GPP_function_for_optimization(env = env, tr = tr_optim, pars = pars)
  #GPP_out_optim =  GPP_function_for_optimization(env = env, tr = tr_optim, pars = pars)
  R_maintenance_optim =  R_maintenance_function(env = env, tr = tr_optim, pars = pars, GPP_out = GPP_out_optim)
  NPP_optim =  calc_NPP(GPP_out = GPP_out_optim, R_maintenance = R_maintenance_optim, pars = pars)
  C_turnover_total_optim =  C_turnover_function(env = env, tr = tr_optim, pars = pars)
  Net_C_gain_optim =  calc_net_C_gain(NPP = NPP_optim, C_turnover_total = C_turnover_total_optim)

  # Rule out all situations in which there is no positive NPP
  if NPP_optim < 0 
    Net_C_gain_optim = 100*tr_optim.H
  end

  # For all situations with productivity, find the ones where Net_C_gain is zero at the maximum height
  return abs(Net_C_gain_optim) - 10 * tr_optim.H

end

### Run Optimizer 
trait_optimizer = function(; env, pars,trait_optimization_function, iters)

    # Leaf C and longevity optimization 
    pass_params = (
      env = env, 
      pars = pars,
    )
  
    f = OptimizationFunction(trait_optimization_function)
    x0 = [1.0, 500, 1]
    prob = Optimization.OptimizationProblem(f, x0, pass_params, lb = [0.2, 50.0, 1], ub = [5.0, 5000.0, 50])
    a_ll_optimized, C_leaf_optimized, H_optimized = Optimization.solve(prob, ECA(), maxiters = iters, abstol = 50)
      
    # Seasonality optimization 
    if (sum(env.tair_monthly .<= pars.temp_threshold_growing_season) >= 3. && a_ll_optimized .< 1.)
      seasonality_optimized =  0.
    else
      seasonality_optimized =  1.
    end
  
    # r_s_r 
    r_s_r_optimized = 1.
  
  
    optimized_traits = (
      a_ll_optimized = a_ll_optimized,
      C_leaf_optimized = C_leaf_optimized, 
      r_s_r_optimized = r_s_r_optimized,
      H_optimized = H_optimized, 
      seasonality_optimized = seasonality_optimized
    )
  
    return(optimized_traits)
end

### Trait evolution and height function 
trait_evolution = function(; optimized_traits, env, tr, pars, evorate)

  new_C_leaf = tr.C_leaf + evorate * (optimized_traits.C_leaf_optimized - tr.C_leaf)

  new_r_s_r = tr.r_s_r + evorate * (optimized_traits.r_s_r_optimized - tr.r_s_r)

  # Is seasonality an option? 
  if sum(env.tair_monthly .<= pars.temp_threshold_growing_season) >= 3.
    potential_seasonality =  0
  else
    potential_seasonality =  1
  end

  new_a_ll = tr.a_ll + evorate * (optimized_traits.a_ll_optimized - tr.a_ll)

  if new_a_ll < 1 && potential_seasonality == 0
    new_seasonality = 0
  else 
    new_seasonality = 1
  end


  # Update temperature niche 
  new_Tave_optim = tr.Tave_optim + evorate * (env.tair_annual - tr.Tave_optim)
  new_Tmin_optim = tr.Tmin_optim + evorate * (minimum(env.tair_monthly) - tr.Tmin_optim)
  new_Tmax_optim = tr.Tmax_optim + evorate * (maximum(env.tair_monthly) - tr.Tmax_optim)
  new_Pave_optim = tr.Pave_optim + evorate * (env.precip_annual - tr.Pave_optim)

  tr_new = (
    H = NaN, 
    a_ll = new_a_ll, 
    C_leaf = new_C_leaf, 
    seasonality = new_seasonality, 
    r_s_r = new_r_s_r, 
    Tave_optim = new_Tave_optim, 
    Tmax_optim = new_Tmax_optim, 
    Tmin_optim = new_Tmin_optim, 
    Pave_optim = new_Pave_optim
  )

  return tr_new

end

### Daylength calculation 
# Check this function
daylength_calculation = function(latitude, DOI)

  delta = 23.44 * sin( 2*pi / 365 * (DOI - 81))

  P = -tan(latitude * pi / 180) * tan(delta * pi / 180)

  if P >= 1
    daylength = 0
  elseif P <= -1
    daylength = 24
  else 
    daylength = 24/pi * acos(P)
  end

  return daylength 

end

### Climate field rotation 
rotate_raster = function(input_raster)

  x_resolution = length(lookup(input_raster, X))
  y_resolution = length(lookup(input_raster, Y))
  interval_x = 360 / x_resolution
  interval_y = 180 / y_resolution

  ti = Ti(1:length(input_raster[1,1,:]))

  lon = X(-((x_resolution/2)-0.5 * interval_x):interval_x:((x_resolution/2)-0.5 * interval_x))
  lat = Y(-((y_resolution/2)-0.5 * interval_y):interval_y:((y_resolution/2)-0.5 * interval_y))


  left = parent(input_raster[1:Int(x_resolution/2), :, :])
  right = parent(input_raster[Int(x_resolution/2+interval_x):Int(x_resolution), :, :])
  switched = cat(right, left, dims = 1)

  new_raster = Raster(switched, dims = (lon, lat, ti))

  # lon, lats are now correct - still need to reverse by Y dimension for correct indexing

  new_raster = reverse(new_raster, dims = 2)

  
  return(new_raster)

end

# ### Approximate radiation components 
# tair_mar = tair[Ti = 3]

# latitude_grid = deepcopy(rls[Ti = 3])
# longitude_grid = deepcopy(rls[Ti = 3])
# for i = 1:length(lookup(latitude_grid, X))
#   for j = 1:length(lookup(latitude_grid, Y))

#     lon = lookup(latitude_grid, X)[i]
#     lat = lookup(latitude_grid, Y)[j]

#     latitude_grid[i, j] = lat
#     longitude_grid[i, j] = lon


#   end
# end

# # empirical constant
# b = 0.2
# A = 107
# cloudiness = 50 .* precip[Ti = 3]
# ni = 1 .- cloudiness # ni = 1 - cloudiness
# ni = ni[Y(-40 .. 40)]

# # Net longwave radiation 
# # positive outgoing 
# rl = (b .+ (1 - b) .* ni) .* (A .- tair_mar)
# plot(rl)

# # Instantaneous solar irradiance at surface 
# c = 0.25 
# d = 0.5
# Q00 = 1360 # Solar constant 1360 W m-2
# i = 73 # Day of the year
# Q0 = Q00  * (1 + 2 * 0.01675 * cos(2 * pi * i / 365))

# del = -23.4 * (pi / 180) * cos(2 * pi * (i + 10) / 365)

# t = 0 # 12 am
# cos_z = sin.(latitude_grid .* pi/180) .* sin.(del) .+ cos.(latitude_grid .* pi/180) .* cos(del) .* cos(2 * pi * t/24) # not sure about last cos(0)

# # Instantaneous solar irradiance, depends on time of the day t
# rs = (c .+ d .* ni) .* Q0 .* cos_z
# plot(rs)

# # To get daytime solar irradiance, rs needs to be integrated from sunrise to sunset 

# h05 = acos.(-(sin.(latitude_grid .* pi/180) .* sin.(del)) ./ (cos.(latitude_grid .* pi/180) .* cos.(del)))

# k = 12 / pi * 3600
# rs_day = 2 * (c .+ d .* ni) .* Q0 .* (sin.(latitude_grid .* pi/180) .* sin(del) .* h05 + cos.(latitude_grid .* pi/180) .* cos(del) .* sin.(h05)) * k
# plot(rs_day)

# # Par
# plot(rs_day .* 0.5 .* (4.6e-6))
# plot(rsds[Ti = 3] .* 0.5 .* 4.6)

# plot(rsds[Ti = 3] ./ (60 .* 60 .* 24 .* 1e-6))

raster_area = function(raster)

  area_raster = deepcopy(raster)
  for x in lookup(raster, X)
    for y in lookup(raster, Y)

      # Get resolution 
      x_resolution, y_resolution = abs.(step.(dims(raster, (X,Y))))

      # Latitude distance is constant 
      latitude_distance = ( (pi * 6371.008) / 180 ) * y_resolution

      # Longitude distance depends on latitude 
      longitude_distance = ( (pi * 6371.008) / 180 ) * cos(y * (pi/180)) * x_resolution

      area_raster[At(x), At(y)] = latitude_distance * longitude_distance
    end
  end

  return(area_raster)
end