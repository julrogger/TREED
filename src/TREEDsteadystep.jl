###########################################################
########## TREEDsteadystep function #######################
###########################################################

# Author: Julian Rogger
# Startdate: 04.2025
# Log: 

"""
    TREEDsteadystep(; tair::Raster, precip::Raster, clt::Raster, rsds::Raster, topo::Raster, CO2::Float64, res::Float64, FDsampling::Bool, RIsampling::Bool, outputdir::String)

    Calculates a steady state vegetation distribution for the given topography and climatic inputs.

    Arguments: 
    - tair: Raster of monthly average temperatures in K
    - precip: Raster of monthly average precipitation in m/s 
    - clt: Raster of monthly average cloud cover in fraction 0-1
    - rsds: Raster of monthly average downwelling shortwave radiation in W m-2
    - topo: Raster of topography 
    - CO2: Current atmospheric CO2 concentration in ppm 
    - res: Target resolution of the model in arcdegrees, will be applied in longitude and latitude
    - FDsampling: If true, assess volume of functional space at each location. Computationally heavy.
    - RIsampling: If true, assess richness potential based on functional diversity and landscape complexity. Computationally heavy. 
    - RI_landscape_window: Size of the landscape window for calculating diversity indices. In km. Default to 300 km. 
    - outputdir: path to store outputs
    - kwargs...: Optional named tuple for changes in the TREED.pars parameter list.

    All Raster inputs need to match in resolution and orientation. Orientation of rasters [X = longitude, Y = latitude].
"""
function TREEDsteadystep(;tair, precip, clt, rsds, topo, CO2, res, FDsampling, RIsampling, RI_landscape_window=300.0, outputdir,  kwargs...)

    ############################################################
    ### 0) Include changes in paramter file  ###################
    ############################################################
    pars2 = merge(pars, kwargs)

    ############################################################
    ### 1) Get climate/topo inputs  ############################
    ############################################################

    climate = create_TREED_climate_input(tair, precip, clt, rsds, topo, CO2, res)

    ############################################################
    ### 2) Initialize traits  ##################################
    ############################################################

    traits_start = initialize_TREED_traits(climate)

    ############################################################
    ### 3) Run trait optimization  #############################
    ############################################################

    println(string("Starting optimization function"))

    traits_optimized = run_TREED_optimization(traits_start, climate, pars2)

    println("Done with optimization")

    ############################################################
    ### ASSUME TRAITS OPTIMIZED = TRAITS EVOLVED  ##############
    ############################################################

    ############################################################
    ### 4) calculate fluxes with final traits distribution  ####
    ############################################################

    # Skipping ecological and evolutionary functionalities 
    traits_end = traits_optimized
    vegetation_output = run_TREED_final_distribution(traits_end, climate, pars2)

    
    if FDsampling == true || RIsampling == true 

        ############################################################
        ### 5) sample functional diversity   #######################
        ############################################################

        println("Starting functional trait sampling")

        height_space = range(1, 40, length = 10)
        C_leaf_space = range(50, 5000, length = 8)
        a_ll_space = range(0.2, 5, length = 5)

        functional_diversity_record = run_TREED_functional_diversity_sampling(traits_end, climate, height_space, C_leaf_space, a_ll_space, pars2)

        println("Done with functional trait sampling")

        if RIsampling == true

            ############################################################
            ### 6) Assess species richness potential   #################
            ############################################################

            println("Starting richness estimation")

            richness_estimation = run_TREED_richness_assessment(functional_diversity_record, vegetation_output, climate, RI_landscape_window)

            println("Done with richness estimation")
        end

    end


    ############################################################
    ### Save TREED output ######################################
    ############################################################

    if FDsampling == true 
        functional_diversity = (functional_diversity = functional_diversity_record, )
        vegetation_output = merge(vegetation_output, functional_diversity)
    end

    if RIsampling == true 
        vegetation_output = merge(vegetation_output, richness_estimation)
    end

    topography = (topography = climate.topo, )
    vegetation_output = merge(vegetation_output, topography)
    TREED_output = RasterStack(vegetation_output)

    if !isdir(outputdir)
        mkpath(outputdir)
    end

    write(string(outputdir,"/TREED_output.nc"), TREED_output, missingval=-Inf32, force=true)
    println("Done, everything ok!")

    return(TREED_output)

end