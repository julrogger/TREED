###########################################################
########## TREEDnonsteadystep function ####################
###########################################################

# Author: Julian Rogger
# Startdate: 04.2025
# Log: 

"""
    TREEDnonsteadystep(; tair::Raster, precip::Raster, clt::Raster, rsds::Raster, topo::Raster, CO2::Float64, evorate::Float64, dispersal::Float64, res::Float64, outputdir::String, traitdir::String)

    Calculate vegetation distribution for a given climatic and topographic input, considering eco-evolutionary adaptation constraints and a starting trait distribution. 
    Made for model coupling.

    Arguments: 
    - tair: Raster of monthly average temperatures in K for each timestep 
    - precip: Raster of monthly average precipitation in m/s for each timestep
    - clt: Raster of monthly average cloud cover in fraction 0-1 for each timestep 
    - rsds: Raster of monthly average downwelling shortwave radiation in W m-2 for each timestep
    - topo: Raster of topography for each timestep
    - CO2: Atmospheric CO2 concentration in ppm 
    - evorate: Rate of trait evolution, fraction. Trait evolved = alpha * (trait optimum - trait current).
    - dispersal: Rate of dispersal. Defines moving window radius to be searched for competing vegetation units. In km. 
    - res: Target resolution of the model in arcdegrees, will be applied in longitude and latitude
    - outputdir: path to store outputs
    - traitdir: path to access initial trait distribution. If no trait object (traitdir/traits_restart.jld2) exists, a default will be created.

    All Raster inputs need to match in resolution and orientation. Orientation of rasters [X = longitude, Y = latitude].
"""
function TREEDnonsteadystep(; tair, precip, clt, rsds, topo, CO2, evorate, dispersal, res, outputdir, traitdir)

    ############################################################
    ### 1) Get climate/topo inputs  ############################
    ############################################################

    climate = create_TREED_climate_input(tair, precip, clt, rsds, topo, CO2, res)

    ############################################################
    ### 2) Initialize traits  ##################################
    ############################################################

    # Create a trait directory if it does not exist yet
    mkpath(traitdir)

    # If the trait directory is empty, create a default
    if !isfile(string(traitdir,"/traits_restart.jld2"))
        println("No starting trait distribution exists. Creating default.")
        traits_start = initialize_TREED_traits(climate)
    else
        # Otherwise read from trait directory
        println("Reading initial trait distribution from restart file")
        data = load(string(traitdir,"/traits_restart.jld2"))
        traits_start = data["traits_end"]
    end


    ############################################################
    ### 3) Run trait optimization  #############################
    ############################################################

    println(string("Starting optimization function for timestep ", timestep))

    traits_optimized = run_TREED_optimization(traits_start, climate)

    println("Done with optimization")

    ############################################################
    ### 4) Run trait evolution  ################################
    ############################################################

    traits_evolved = run_TREED_evolution(traits_start, traits_optimized, climate, evorate)

    ############################################################
    ### 5) Run dispersal and competition  ######################
    ############################################################

    traits_end = run_TREED_ecology(traits_start, traits_evolved, traits_optimized, climate, dispersal)

    # Save end trait distribution for restart 
    @save string(traitdir,"/traits_restart.jld2") traits_end

    ############################################################
    ### 6) calculate fluxes with final traits distribution  ####
    ############################################################

    vegetation_output = run_TREED_final_distribution(traits_end, climate)


    ############################################################
    ### Save TREED output ######################################
    ############################################################

    topography = (topography=climate.topo,)
    vegetation_output = merge(vegetation_output, topography)
    TREED_output = RasterStack(vegetation_output)

    if !isdir(outputdir)
        mkpath(outputdir)
    end

    write(string(outputdir, "/TREED_output_timestep_", timestep, ".nc"), TREED_output, missingval=-Inf32, force=true)
    println("Done with timestep ", timestep, ", everything ok!")

    return(TREED_output)

end
