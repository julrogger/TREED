###########################################################
########## TREED Julia implementation V0.0 ################
########## Parameter file #################################
###########################################################

pars = (

    # Fixed model parameters: 
    # Favorable growing season 
    temp_threshold_growing_season = 3, # min [°C]

    # Plant allometry 
    WD_sapwood = 2.5e+5,      # Constant average density of sapwood
    WD_heartwood = 2.5e+5,     # Constant average density of heartwood 
    k_la_sa = 4000, #8000 # 4000 Schaphoff 2018 #8000 (Sitch 2003) # Ratio of sapwood area needed to maintain leaf area 
    k_allom1 = 75.00093,                            
    k_allom2 = 49.99979,                          
    k_allom3 = 0.6000011,          # Sitch #0.67 # Schaphoff
    k_rp =  1.599998,   
    k_density_intercept = 194.0783,        
    k_power = 1.2255,               
    CA_max = 15, 
    CA_min = 1,

    # Photosynthesis model 
    p_atm = 101300, # Atmospheric pressure [Pa]
    pO2 = 20900, # O2 partial pressure [Pa]
    T_base = 5,
    alpha_c3 = 0.08, # Intrinsic quantum efficiency of CO2 uptake in C3 plants 
    alpha_leaf_to_stand = 0.55, # Efficiency of light absorption at ecosystem level compared to leaf level 
    c_mass = 12.01, # Molar mass of C, g/mol 
    lambdamc3 = 0.8,
    theta = 0.7, # shape parameter of photosynthesis model 
    month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
    temp_crit_photosynthesis = 0,

    # Respiration 
    cn_leaf = 29,
    cn_wood = 330,
    cn_root = 29,
    b_c3 = 0.015, # Leaf respiration as a fraction of Rubisco capacity in C3 plants 
    k = 0.0548, # Sprugel et al. 1995
    r = 1.2, # Base respiration rate, will be lower in warmer temperature (i.e. trop. trees have value ~ 0.2) --> implemented as a function of temp
    r_gr = 0.25, # % of NPP that is growth respiration

    # Carbon turnover 
    f_sapwood = (1/15),
    f_heartwood = (1/15), # turnover times in agreement with vegetation models: Pugh et al. 2020 Biogeosciences
    f_coarseroot = (1/15),
    r_repr = 0.10, # Fixed fraction of NPP that is invested into reproduction, as well as other unaccounted carbon investments (i.e., symbionts)


    # Niche change stress parameter
    temp_niche_breadth_parameter = 0.02, # Approx. 15°C deviation from habitat before loss of productivity
)

