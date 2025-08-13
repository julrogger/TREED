###########################################################
########## TREEDnonsteadycontinuous function ##############
###########################################################

# Author: Julian Rogger
# Startdate: 04.2025
# Log: 

"""
    TREEDnonsteadycontinuous(; tairvec::Vector{Raster}, precipvec::Vector{Raster}, cltvec::Vector{Raster}, rsdsvec::Vector{Raster}, topovec::Vector{Raster}, CO2vec::Vector{Float64}, evorate::Float64, startinsteady::Bool, dispersal::Float64, res::Float64, FDsampling::Bool, outputdir::String)

    Calculate vegetation distribution across a series of topographic and climatic inputs, considering eco-evolutionary adaptation constraints. 

    Arguments: 
    - tairvec: Vector of rasters of monthly average temperatures in K for each timestep 
    - precipvec: Vector of rasters of monthly average precipitation in m/s for each timestep
    - cltvec: Vector of rasters of monthly average cloud cover in fraction 0-1 for each timestep 
    - rsdsvec: Vector of rasters of monthly average downwelling shortwave radiation in W m-2 for each timestep
    - topovec: Vector of rasters of topography for each timestep
    - CO2vec: Vector of atmospheric CO2 concentrations in ppm 
    - evorate: Rate of trait evolution, fraction. Trait evolved = alpha * (trait optimum - trait current).
    - startinsteady: Run first model iteration with evorate = 1, resulting in start trait distribution = optimum distribution.
    - dispersal: Rate of dispersal. Defines moving window radius to be searched for competing vegetation units. In km. 
    - res: Target resolution of the model in arcdegrees, will be applied in longitude and latitude
    - outputdir: path to store outputs

    All Raster inputs need to match in resolution and orientation. Orientation of rasters [X = longitude, Y = latitude].
"""
function TREEDnonsteadycontinuous(;tairvec, precipvec, cltvec, rsdsvec, topovec, CO2vec, evorate, startinsteady, dispersal, res, outputdir)

    traits_end = ()
    for timestep in eachindex(tairvec)
        
        ############################################################
        ### 1) Get climate/topo inputs  ############################
        ############################################################
        tair = tairvec[timestep]
        precip = precipvec[timestep]
        clt = cltvec[timestep]
        rsds = rsdsvec[timestep]
        topo = topovec[timestep]
        CO2 = CO2vec[timestep]

        climate = create_TREED_climate_input(tair, precip, clt, rsds, topo, CO2, res)

        ############################################################
        ### 2) Initialize traits  ##################################
        ############################################################
        
        if timestep == 1
            # If new simulation - initialize traits
            traits_start = initialize_TREED_traits(climate)
        else 
            # If continued simulation - initial traits of timestep from final trait distribution at t-1
            traits_start = traits_end
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

        if timestep == 1 & startinsteady == true
            # First timestep set evorate to 1 -> update start trait distribution to initial optimum steady state
            traits_evolved = run_TREED_evolution(traits_start, traits_optimized, climate, 1)
        else 
            traits_evolved = run_TREED_evolution(traits_start, traits_optimized, climate, evorate)
        end

        ############################################################
        ### 5) Run dispersal and competition  ######################
        ############################################################

        traits_end = run_TREED_ecology(traits_start, traits_evolved, traits_optimized, climate, dispersal)

        ############################################################
        ### 6) calculate fluxes with final traits distribution  ####
        ############################################################

        vegetation_output = run_TREED_final_distribution(traits_end, climate)


        ############################################################
        ### Save TREED output ######################################
        ############################################################

        topography = (topography = climate.topo, )
        vegetation_output = merge(vegetation_output, topography)
        TREED_output = RasterStack(vegetation_output)

        if !isdir(outputdir)
            mkpath(outputdir)
        end

        write(string(outputdir,"/TREED_output_timestep_",timestep,".nc"), TREED_output, missingval=-Inf32, force=true)
        println("Done with timestep ",timestep,", everything ok!")

    end
end
