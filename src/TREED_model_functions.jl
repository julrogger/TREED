    
###########################################################
########## TREED model functions  #########################
###########################################################

# Climate input struct
struct TREED_climate_input

    precip::Raster
    tair::Raster
    rsds::Raster
    rss::Raster
    rls::Raster
    daylength::Raster
    CO2_ppm::Float64
    habitable::Raster
    topo::Raster

end

# Plant traits struct 
struct traits 

    H::Raster
    a_ll::Raster
    C_leaf::Raster
    seasonality::Raster
    r_s_r::Raster
    Tave_optim::Raster
    Tmax_optim::Raster
    Tmin_optim::Raster
    Pave_optim::Raster

end

# Get and format climate and topographic inputs
create_TREED_climate_input = function(tair, precip, clt, rsds, topo, CO2, res)

    # Create template raster
    lon = X(Sampled(-180+(res/2):res:180-(res/2); sampling=Intervals(Center())))
    lat = Y(Sampled(90-(res/2):-res:-90+(res/2); sampling=Intervals(Center())))
    template_raster = Raster(rand(lon, lat))

    # Handle tair
    tair = replace_missing(tair, -Inf32)
    tair = convert.(Float32, tair)
    tair = resample(tair; to = template_raster)
    tair = tair .- 273.15 # °K --> °C

    # Handle precip 
    precip = replace_missing(precip, -Inf32)
    precip = convert.(Float32, precip)
    precip = resample(precip; to = template_raster)
    precip = precip .* (60*60*24) # m/s --> m/day

    # Handle clt 
    clt = replace_missing(clt, -Inf32)
    clt = convert.(Float32, clt)
    clt = resample(clt; to = template_raster)
    
    # Handle rsds 
    rsds = replace_missing(rsds, -Inf32)
    rsds = convert.(Float32, rsds)
    rsds = resample(rsds; to = template_raster)
    
    # Handle topo 
    topo = replace_missing(topo, -Inf32)
    topo = convert.(Float32, topo)
    topo = resample(topo; to=template_raster)
    
    #### Calculate radiation components 
    # Land-sea albedo map (I do not consider ice here)
    albedo = deepcopy(topo)
    albedo[topo .< 0] .= 0.069
    albedo[topo .> 0] .= 0.15

    rss = deepcopy(rsds)
    for m in 1:12
        rss[Ti = m] .= rsds[Ti = m] .* (1 .- albedo)
    end

    rls = deepcopy(rsds)
    # Using formula from LPJ Schaphoff et al. 2018
    b = 0.2; A = 107;
    for m in 1:12
        ni = 1 .- clt[Ti = m]
        rls[Ti = m] = (b .+ (1-b) .* ni) .* (A .- tair[Ti = m])
    end

    rss = rss .* (60*60*24*1e-6) # J s-1 m-2 --> MJ day-1 m-2
    rsds = rsds .* (60*60*24*1e-6) # J s-1 m-2 --> MJ day-1 m-2
    rls = rls .* (60*60*24*1e-6) .* (-1) # J s-1 m-2 --> MJ day-1 m-2, .* (-1) --> downward positive

    # Create daylength input
    daylength = deepcopy(precip)
    DOIs = [15, 46, 73, 104, 134, 166, 196, 227, 259, 288, 319, 349]
    for m in 1:12
        for j = 1:length(lookup(daylength, Y))
            latitude = lookup(daylength, Y)[j]
            hours = daylength_calculation(latitude, DOIs[m])
            if hours < 1
                hours = 1
            end
            daylength[:, j, m] .= hours
        end
    end

    # Create land sea mask 
    lsm = deepcopy(topo)
    lsm[topo .< 0] .= 0
    lsm[topo .>= 0] .= 1
    lsm[isinf.(tair[:,:,1])] .= 0


    # Create habitability mask for speed up 
    habitable = deepcopy(topo)
    habitable[habitable .< 0] .= NaN
    habitable[habitable .>= 0] .= 1
    habitable[isinf.(tair[:,:,1])] .= NaN

    climate = TREED_climate_input(precip, tair, rsds, rss, rls, daylength, CO2, habitable, topo)

    return(climate)

end

# Initialize trait map 
initialize_TREED_traits = function(climate)

    H_initial = deepcopy(climate.habitable)
    H_initial[climate.habitable .== 1] .= 0.5
    
    a_ll_initial = deepcopy(climate.habitable)
    a_ll_initial[climate.habitable .== 1] .= 1.0

    C_leaf_initial = deepcopy(climate.habitable)
    C_leaf_initial[climate.habitable .== 1] .= 500

    seasonality_initial = deepcopy(climate.habitable)
    seasonality_initial[climate.habitable .== 1] .= 1.0

    r_s_r_initial = deepcopy(climate.habitable)
    r_s_r_initial[climate.habitable .== 1] .= 1.0

    Tave_optim_initial = mean(climate.tair, dims=Ti)

    Tmax_optim_initial = maximum(climate.tair, dims = Ti)

    Tmin_optim_initial = minimum(climate.tair, dims = Ti)

    Pave_optim_initial = mean(climate.precip, dims = Ti)


    traits_start = traits(H_initial, a_ll_initial, C_leaf_initial, seasonality_initial, r_s_r_initial,
                          Tave_optim_initial, Tmax_optim_initial, Tmin_optim_initial, Pave_optim_initial)


end

# Run the optimization loop 
run_TREED_optimization = function(traits_start, climate, pars2)

    # Initialize object to write optimization results 
    traits_optimized = deepcopy(traits_start)

    # Run loop on multiple threads 
    #     
 Threads.@threads for i = 1:length(lookup(traits_start.H, X))
        for j = 1:length(lookup(traits_start.H, Y))
            if climate.habitable[i, j] == 1

                # Extract climate at location 
                env = (precip_monthly = parent(climate.precip[i, j, :]), 
                tair_monthly = parent(climate.tair[i, j, :]), 
                tair_annual = mean(parent(climate.tair[i, j, :])), 
                rsds_monthly = parent(climate.rsds[i, j, :]), 
                rss_monthly = parent(climate.rss[i, j, :]), 
                rls_monthly = parent(climate.rls[i, j, :]), 
                daylength = parent(climate.daylength[i, j, :]), 
                CO2_ppm = climate.CO2_ppm, 
                precip_annual = mean(parent(climate.precip[i, j, :]))
                )

                # Run optimization to get maximimising H, C_leaf and a_ll
                optimized_traits = trait_optimizer(env = env, par = pars2, trait_optimization_function = trait_optimization_function, iters = 20)

                # Save optimization results 
                traits_optimized.H[i, j] = optimized_traits.H_optimized
                traits_optimized.a_ll[i, j] = optimized_traits.a_ll_optimized
                traits_optimized.C_leaf[i, j] = optimized_traits.C_leaf_optimized
                traits_optimized.r_s_r[i, j] = optimized_traits.r_s_r_optimized
                traits_optimized.seasonality[i, j] = optimized_traits.seasonality_optimized
                traits_optimized.Tave_optim[i, j] = env.tair_annual
                traits_optimized.Tmax_optim[i, j] = maximum(env.tair_monthly)
                traits_optimized.Tmin_optim[i, j] = minimum(env.tair_monthly)
                traits_optimized.Pave_optim[i, j] = env.precip_annual

            end
        end
    end

    # Return optimized trait landscape 
    return(traits_optimized)
end


# Run the evolution loop
run_TREED_evolution = function(traits_start, traits_optimized, climate, evorate, pars2)

    # Initialize object to write optimization results 
    traits_evolved = deepcopy(traits_start)

    # Run loop on multiple threads 
    Threads.@threads for i = 1:length(lookup(traits_start.H, X))
        for j = 1:length(lookup(traits_start.H, Y))
            if climate.habitable[i, j] == 1

                # Extract current traits at location 
                tr = (H=traits_start.H[i, j],
                    a_ll=traits_start.a_ll[i, j],
                    C_leaf=traits_start.C_leaf[i, j],
                    seasonality=traits_start.seasonality[i, j],
                    r_s_r=traits_start.r_s_r[i, j],
                    Tave_optim=traits_start.Tave_optim[i, j],
                    Tmax_optim=traits_start.Tmax_optim[i, j],
                    Tmin_optim=traits_start.Tmin_optim[i, j],
                    Pave_optim=traits_start.Pave_optim[i, j]
                )

                # Extract optimized traits at location 
                optimized_traits = (
                    a_ll_optimized=traits_optimized.a_ll[i, j],
                    C_leaf_optimized=traits_optimized.C_leaf[i, j],
                    r_s_r_optimized=traits_optimized.r_s_r[i, j],
                    H_optimized=traits_optimized.H[i, j],
                    seasonality_optimized=traits_optimized.seasonality[i, j]
                )

                # Extract climate at location 
                env = (precip_monthly=parent(climate.precip[i, j, :]),
                    tair_monthly=parent(climate.tair[i, j, :]),
                    tair_annual=mean(parent(climate.tair[i, j, :])),
                    rsds_monthly=parent(climate.rsds[i, j, :]),
                    rss_monthly=parent(climate.rss[i, j, :]),
                    rls_monthly=parent(climate.rls[i, j, :]),
                    daylength=parent(climate.daylength[i, j, :]),
                    CO2_ppm=climate.CO2_ppm,
                    precip_annual=mean(parent(climate.precip[i, j, :]))
                )

                # Run evolution
                new_traits = trait_evolution(optimized_traits=optimized_traits, env=env, tr=tr, par=pars2, evorate=evorate)

                # Write new traits to output array 
                traits_evolved.H[i, j] = new_traits.H
                traits_evolved.a_ll[i, j] = new_traits.a_ll
                traits_evolved.C_leaf[i, j] = new_traits.C_leaf
                traits_evolved.seasonality[i, j] = new_traits.seasonality
                traits_evolved.r_s_r[i, j] = new_traits.r_s_r
                traits_evolved.Tave_optim[i, j] = new_traits.Tave_optim
                traits_evolved.Tmax_optim[i, j] = new_traits.Tmax_optim
                traits_evolved.Tmin_optim[i, j] = new_traits.Tmin_optim
                traits_evolved.Pave_optim[i, j] = new_traits.Pave_optim

            end
        end
    end

    # Return evolved trait landscape 
    return(traits_evolved)
end



# Run the dispersal and competition loop
run_TREED_ecology = function(traits_start, traits_evolved, traits_optimized, climate, dispersal, pars2)

    # Initialize object to write optimization results 
    traits_end = deepcopy(traits_evolved)

    # Circular arrays for competition 
    traits_evolved_circular = (H=CircularArray(traits_evolved.H),
        a_ll=CircularArray(traits_evolved.a_ll),
        C_leaf=CircularArray(traits_evolved.C_leaf),
        seasonality=CircularArray(traits_evolved.seasonality),
        r_s_r=CircularArray(traits_evolved.r_s_r),
        Tave_optim=CircularArray(traits_evolved.Tave_optim[Ti=1]), ## To do: why do I need to put Ti = 1 to remove time dimension, should already be without time dimensions from initialization
        Tmax_optim=CircularArray(traits_evolved.Tmax_optim[Ti=1]),
        Tmin_optim=CircularArray(traits_evolved.Tmin_optim[Ti=1]),
        Pave_optim=CircularArray(traits_evolved.Pave_optim[Ti=1])
    )

    traits_optimized_circular = (H=CircularArray(traits_optimized.H),
        a_ll=CircularArray(traits_optimized.a_ll),
        C_leaf=CircularArray(traits_optimized.C_leaf),
        seasonality=CircularArray(traits_optimized.seasonality),
        r_s_r=CircularArray(traits_optimized.r_s_r),
        Tave_optim=CircularArray(traits_optimized.Tave_optim[Ti=1]),
        Tmax_optim=CircularArray(traits_optimized.Tmax_optim[Ti=1]),
        Tmin_optim=CircularArray(traits_optimized.Tmin_optim[Ti=1]),
        Pave_optim=CircularArray(traits_optimized.Pave_optim[Ti=1])
    )

    traits_start_circular = (H=CircularArray(traits_start.H),
        a_ll=CircularArray(traits_start.a_ll),
        C_leaf=CircularArray(traits_start.C_leaf),
        seasonality=CircularArray(traits_start.seasonality),
        r_s_r=CircularArray(traits_start.r_s_r),
        Tave_optim=CircularArray(traits_start.Tave_optim[Ti=1]),
        Tmax_optim=CircularArray(traits_start.Tmax_optim[Ti=1]),
        Tmin_optim=CircularArray(traits_start.Tmin_optim[Ti=1]),
        Pave_optim=CircularArray(traits_start.Pave_optim[Ti=1])
    )

    Threads.@threads for i = 1:length(lookup(traits_optimized.H, X))
        for j = 1:length(lookup(traits_optimized.H, Y))
            if climate.habitable[i, j] == 1
    
                # lon/lat resolution in °
                longitude_resolution = 360 / length(lookup(traits_optimized.H, X))
                latitude_resolution = 180 / length(lookup(traits_optimized.H, Y))
    
                # Get latitude of current location 
                current_latitude = lookup(traits_optimized.H, Y)[j]
    
                # cell to cell longitude distance at current latitude
                cell_to_cell_longitude = (6378 * cos(current_latitude * pi / 180) * pi / 180) * longitude_resolution
                cell_to_cell_latitude = 111.11 * latitude_resolution

                # dispersal window 
                dispersal_cell_range_west = round(Int, dispersal / cell_to_cell_longitude)
                dispersal_cell_range_north = round(Int, dispersal / cell_to_cell_latitude)

                # Check for maximum one surrounding
                if dispersal_cell_range_west >= round(Int, length(lookup(climate.habitable, X)) / 2)
                    dispersal_cell_range_west = round(Int, length(lookup(climate.habitable, X)) / 2)
                end

                # Current height within dispersal buffer 
                current_H_buffer = traits_start_circular.H[i-dispersal_cell_range_west:i+dispersal_cell_range_west, j-dispersal_cell_range_north:j+dispersal_cell_range_north]
                buffer_center_index_i = dispersal_cell_range_west + 1
                buffer_center_index_j = dispersal_cell_range_north + 1

                # Search for most optimal trait combination from geographically limited trait space 
                # Leaf longevity
                a_ll_buffer = traits_evolved_circular.a_ll[i-dispersal_cell_range_west:i+dispersal_cell_range_west, j-dispersal_cell_range_north:j+dispersal_cell_range_north]
                a_ll_buffer[current_H_buffer .<= 0] .= NaN # Only consider vegetation units that are productive at the moment
                a_ll_buffer[buffer_center_index_i,buffer_center_index_j] = traits_evolved_circular.a_ll[i,j] # Make sure there is a value in case in completely unproductive landscape
                a_ll_minimizer = argmin(replace!(abs.(a_ll_buffer .- traits_optimized_circular.a_ll[i, j]), NaN => Inf))
                traits_end.a_ll[i, j] = a_ll_buffer[a_ll_minimizer]
    
                # Seasonality (bound to leaf longevity)
                seasonality_buffer = traits_evolved_circular.seasonality[i-dispersal_cell_range_west:i+dispersal_cell_range_west, j-dispersal_cell_range_north:j+dispersal_cell_range_north]
                if a_ll_buffer[a_ll_minimizer] <= 1.0 && seasonality_buffer[a_ll_minimizer] < 1.0
                    traits_end.seasonality[i,j] = 0.0
                else
                    traits_end.seasonality[i,j] = 1.0
                end


               # Leaf carbon 
                C_leaf_buffer = traits_evolved_circular.C_leaf[i-dispersal_cell_range_west:i+dispersal_cell_range_west, j-dispersal_cell_range_north:j+dispersal_cell_range_north]
                C_leaf_buffer[current_H_buffer .<= 0] .= NaN # Only consider vegetation units that are productive at the moment
                C_leaf_buffer[buffer_center_index_i,buffer_center_index_j] = traits_evolved_circular.C_leaf[i,j] # Make sure there is a value in case in completely unproductive landscape
                C_leaf_minimizer = argmin(replace!(abs.(C_leaf_buffer .- traits_optimized_circular.C_leaf[i, j]), NaN=>Inf))
                traits_end.C_leaf[i, j] = C_leaf_buffer[C_leaf_minimizer]
    
                # Tave_optim
                Tave_optim_buffer = traits_evolved_circular.Tave_optim[i-dispersal_cell_range_west:i+dispersal_cell_range_west, j-dispersal_cell_range_north:j+dispersal_cell_range_north]
                Tave_optim_buffer[current_H_buffer .<= 0] .= NaN # Only consider vegetation units that are productive at the moment
                Tave_optim_buffer[buffer_center_index_i,buffer_center_index_j] = traits_evolved_circular.Tave_optim[i,j] # Make sure there is a value in case in completely unproductive landscape
                Tave_optim_minimizer = argmin(replace!(abs.(Tave_optim_buffer .- traits_optimized_circular.Tave_optim[i, j]), NaN=>Inf))
                traits_end.Tave_optim[i, j] = Tave_optim_buffer[Tave_optim_minimizer]
    
                # Tmin_optim
                Tmin_optim_buffer = traits_evolved_circular.Tmin_optim[i-dispersal_cell_range_west:i+dispersal_cell_range_west, j-dispersal_cell_range_north:j+dispersal_cell_range_north]
                Tmin_optim_buffer[current_H_buffer .<= 0] .= NaN # Only consider vegetation units that are productive at the moment
                Tmin_optim_buffer[buffer_center_index_i,buffer_center_index_j] = traits_evolved_circular.Tmin_optim[i,j] # Make sure there is a value in case in completely unproductive landscape
                Tmin_optim_minimizer = argmin(replace!(abs.(Tmin_optim_buffer .- traits_optimized_circular.Tmin_optim[i, j]), NaN=>Inf))
                traits_end.Tmin_optim[i, j] = Tmin_optim_buffer[Tmin_optim_minimizer]
    
                # Tmax_optim
                Tmax_optim_buffer = traits_evolved_circular.Tmax_optim[i-dispersal_cell_range_west:i+dispersal_cell_range_west, j-dispersal_cell_range_north:j+dispersal_cell_range_north]
                Tmax_optim_buffer[current_H_buffer .<= 0] .= NaN # Only consider vegetation units that are productive at the moment
                Tmax_optim_buffer[buffer_center_index_i,buffer_center_index_j] = traits_evolved_circular.Tmax_optim[i,j] # Make sure there is a value in case in completely unproductive landscape
                Tmax_optim_minimizer = argmin(replace!(abs.(Tmax_optim_buffer .- traits_optimized_circular.Tmax_optim[i, j]), NaN=>Inf))
                traits_end.Tmax_optim[i, j] = Tmax_optim_buffer[Tmax_optim_minimizer]
    
                # Pave_optim
                Pave_optim_buffer = traits_evolved_circular.Pave_optim[i-dispersal_cell_range_west:i+dispersal_cell_range_west, j-dispersal_cell_range_north:j+dispersal_cell_range_north]
                Pave_optim_buffer[current_H_buffer .<= 0] .= NaN # Only consider vegetation units that are productive at the moment
                Pave_optim_buffer[buffer_center_index_i,buffer_center_index_j] = traits_evolved_circular.Pave_optim[i,j] # Make sure there is a value in case in completely unproductive landscape
                Pave_optim_minimizer = argmin(replace!(abs.(Pave_optim_buffer .- traits_optimized_circular.Pave_optim[i, j]), NaN=>Inf))
                traits_end.Pave_optim[i, j] = Pave_optim_buffer[Pave_optim_minimizer]
    
    
    
                # Evaluate maximum vegetation height with resulting trait combination
                #####################################################################
    
                # Extract climate at location 
                env = (precip_monthly = parent(climate.precip[i, j, :]), 
                tair_monthly = parent(climate.tair[i, j, :]), 
                tair_annual = mean(parent(climate.tair[i, j, :])), 
                rsds_monthly = parent(climate.rsds[i, j, :]), 
                rss_monthly = parent(climate.rss[i, j, :]), 
                rls_monthly = parent(climate.rls[i, j, :]), 
                daylength = parent(climate.daylength[i, j, :]), 
                CO2_ppm = climate.CO2_ppm, 
                precip_annual = mean(parent(climate.precip[i, j, :]))
                )
    
                find_H_new = function(h)
    
                    tr_int = (
                      H = h, 
                      a_ll = traits_end.a_ll[i,j], 
                      C_leaf = traits_end.C_leaf[i,j], 
                      seasonality = traits_end.seasonality[i,j], 
                      r_s_r = traits_end.r_s_r[i,j],
                      Tave_optim = traits_end.Tave_optim[i,j],
                      Tmax_optim = traits_end.Tmax_optim[i,j],
                      Tmin_optim = traits_end.Tmin_optim[i,j], 
                      Pave_optim = traits_end.Pave_optim[i,j]
                      )
              
                      tr_int = plant_allometry(tr = tr_int, par = pars2)
                      GPP_out_int = GPP_function_for_optimization(env = env, tr = tr_int, par = pars2)
                      R_maintenance_int = R_maintenance_function( env = env, tr  = tr_int, par = pars2, GPP_out = GPP_out_int)
                      C_turnover_total_int = C_turnover_function( env = env, tr = tr_int, par = pars2)
                      NPP_int = calc_NPP(GPP_out = GPP_out_int, R_maintenance = R_maintenance_int, par = pars2)
                      Net_C_gain_int = calc_net_C_gain(NPP = NPP_int, C_turnover_total = C_turnover_total_int, par=pars2)
    
                      if Net_C_gain_int <= -80
                        Net_C_gain_int = -1
                      else
                        Net_C_gain_int = 1
                      end
              
                      return(Net_C_gain_int)
              
                end
                
                test_heights = collect(0:0.25:50)
                evaluations = zeros(length(test_heights))
                evaluations = find_H_new.(test_heights)
                evaluations[ evaluations .> 0 ] .= 1
                evaluations[ evaluations .<= 0 ] .= 0
                
                # If there is no productive trait combination available: 
                # Set height to zero, but keep other traits as from traits end
                # H = 0 will result in non productive location (NPP < 0), and it will not be involved in outwards dispersal
                # Keeping trait values in these unproductive lands is assuming there are still species that are evolving and slowly adapting
                # If environmental conditions return - they will bloom again, in the style of persisting in a refuge
                if sum(evaluations) == 0
                    traits_end.H[i,j] = 0.
                elseif findlast(Int8.([1, 0]), Int8.(evaluations)) === nothing
                    new_H = 50
                    traits_end.H[i,j] = new_H
                else
                    index = findlast(Int8.([1, 0]), Int8.(evaluations))[1]
                    new_H = test_heights[index]
                    traits_end.H[i,j] = new_H
                end
                #####################################################################

            end
        end
    end


    # Return final trait landscape 
    return(traits_end)
end


# Calculate final allometry and fluxes from resulting trait distribution 
run_TREED_final_distribution = function(traits_end, climate, pars2)

    # Initialize output collectors
    output_raster = deepcopy(climate.habitable)

    allometry = (H = deepcopy(output_raster), 
    a_ll = deepcopy(output_raster), 
    C_leaf = deepcopy(output_raster),
    Tave_optim = deepcopy(output_raster), 
    Tmin_optim = deepcopy(output_raster), 
    Tmax_optim = deepcopy(output_raster), 
    C_fineroot = deepcopy(output_raster), 
    C_coarseroot = deepcopy(output_raster),
    C_heartwood = deepcopy(output_raster), 
    C_sapwood = deepcopy(output_raster), 
    CA = deepcopy(output_raster),
    SLA = deepcopy(output_raster), 
    LAI = deepcopy(output_raster), 
    FPC = deepcopy(output_raster),
    seasonality = deepcopy(output_raster),
    r_s_r = deepcopy(output_raster))

    fluxes = (GPP = deepcopy(output_raster), 
    R_maintenance = deepcopy(output_raster),
    NPP = deepcopy(output_raster), 
    Net_C_gain = deepcopy(output_raster), 
    AET = deepcopy(output_raster))

    for i = 1:length(lookup(traits_end.H, X))
        for j = 1:length(lookup(traits_end.H, Y))
            if climate.habitable[i, j] == 1
            

                #print(i / length(lookup(traits_start.H, X)))

                # Extract traits at location 
                tr = (H = traits_end.H[i, j], 
                a_ll = traits_end.a_ll[i, j], 
                C_leaf = traits_end.C_leaf[i, j], 
                seasonality = traits_end.seasonality[i, j], 
                r_s_r = traits_end.r_s_r[i, j], 
                Tave_optim = traits_end.Tave_optim[i, j], 
                Tmax_optim = traits_end.Tmax_optim[i, j], 
                Tmin_optim = traits_end.Tmin_optim[i, j], 
                Pave_optim = traits_end.Pave_optim[i, j])

                # Extract climate at location 
                env = (precip_monthly = parent(climate.precip[i, j, :]), 
                tair_monthly = parent(climate.tair[i, j, :]), 
                tair_annual = mean(parent(climate.tair[i, j, :])), 
                rsds_monthly = parent(climate.rsds[i, j, :]), 
                rss_monthly = parent(climate.rss[i, j, :]), 
                rls_monthly = parent(climate.rls[i, j, :]), 
                daylength = parent(climate.daylength[i, j, :]), 
                CO2_ppm = climate.CO2_ppm, 
                precip_annual = mean(parent(climate.precip[i, j, :]))
                )

                tr_complete = plant_allometry(tr = tr, par = pars2)
                #GPP_output = GPP_function(env = env, tr = tr_complete, par = pars2)
                GPP_output = GPP_function_for_optimization(env = env, tr = tr_complete, par = pars2) # For now I use the same procedure for lambda as in optimization - is more consistent! 
                R_maintenance_output = R_maintenance_function(env = env, tr = tr_complete, par = pars2, GPP_out = GPP_output)
                NPP_output = calc_NPP(GPP_out = GPP_output, R_maintenance = R_maintenance_output, par = pars2)
                C_turnover_total_output = C_turnover_function(env = env, tr = tr_complete, par = pars2)
                Net_C_gain_output = calc_net_C_gain(NPP = NPP_output, C_turnover_total = C_turnover_total_output, par=pars2)



                # Save output
                # Allometry 
                allometry.H[i, j] = tr_complete.H
                allometry.a_ll[i, j] = tr_complete.a_ll
                allometry.C_leaf[i, j] = tr_complete.C_leaf
                allometry.Tave_optim[i, j] = tr_complete.Tave_optim
                allometry.Tmax_optim[i, j] = tr_complete.Tmax_optim
                allometry.Tmin_optim[i, j] = tr_complete.Tmin_optim
                allometry.C_fineroot[i, j] = tr_complete.C_fineroot
                allometry.C_coarseroot[i,j] = tr_complete.C_coarseroot
                allometry.C_heartwood[i, j] = tr_complete.C_heartwood
                allometry.C_sapwood[i, j] = tr_complete.C_sapwood
                allometry.CA[i, j] = tr_complete.CA
                allometry.SLA[i, j] = tr_complete.SLA
                allometry.LAI[i, j] = tr_complete.LAI
                allometry.FPC[i, j] = tr_complete.FPC
                allometry.seasonality[i,j] = tr_complete.seasonality
                allometry.r_s_r[i,j] = tr_complete.r_s_r

                # Fluxes
                fluxes.GPP[i, j] = GPP_output.GPP
                fluxes.R_maintenance[i, j] = R_maintenance_output
                fluxes.NPP[i, j] = NPP_output
                fluxes.Net_C_gain[i, j] = Net_C_gain_output
                fluxes.AET[i, j] = GPP_output.AET

            end
        end
    end

    allometry.H[isnan.(allometry.H).&&climate.habitable.==1] .= 0
    allometry.a_ll[isnan.(allometry.a_ll).&&climate.habitable.==1] .= 0
    allometry.C_leaf[isnan.(allometry.C_leaf).&&climate.habitable.==1] .= 0
    allometry.C_fineroot[isnan.(allometry.C_fineroot).&&climate.habitable.==1] .= 0
    allometry.C_heartwood[isnan.(allometry.C_heartwood).&&climate.habitable.==1] .= 0
    allometry.C_sapwood[isnan.(allometry.C_sapwood).&&climate.habitable.==1] .= 0
    allometry.CA[isnan.(allometry.CA).&&climate.habitable.==1] .= 0
    allometry.LAI[isnan.(allometry.LAI).&&climate.habitable.==1] .= 0
    allometry.FPC[isnan.(allometry.FPC).&&climate.habitable.==1] .= 0

    fluxes.GPP[isnan.(fluxes.GPP).&&climate.habitable.==1] .= 0
    fluxes.R_maintenance[isnan.(fluxes.R_maintenance).&&climate.habitable.==1] .= 0
    fluxes.NPP[isnan.(fluxes.NPP).&&climate.habitable.==1] .= 0
    fluxes.Net_C_gain[isnan.(fluxes.Net_C_gain).&&climate.habitable.==1] .= 0
    fluxes.AET[isnan.(fluxes.AET).&&climate.habitable.==1] .= 0

    # Filter for unproductive land
    allometry.H[fluxes.NPP.<=0 .|| allometry.H.<= 0] .= 0
    allometry.a_ll[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    allometry.C_leaf[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    allometry.C_fineroot[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    allometry.C_heartwood[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    allometry.C_sapwood[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    allometry.CA[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    allometry.LAI[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    allometry.FPC[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0

    fluxes.GPP[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    fluxes.R_maintenance[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    fluxes.Net_C_gain[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    fluxes.AET[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0
    fluxes.NPP[fluxes.NPP.<=0 .|| allometry.H.<=0] .= 0

    vegetation_output = merge(allometry, fluxes)

    return(vegetation_output)
end

# Funcional diversity sampling
run_TREED_functional_diversity_sampling = function(traits_end, climate, H_range, C_leaf_range, a_ll_range, pars2)

    # Initialize output collector 
    functional_diversity_record = deepcopy(climate.habitable)

    # Get combinations from defined range
    combinations = [(h, c, a) for h in H_range,  c in C_leaf_range,  a in a_ll_range ]
    combination_matrix = hcat(map(collect, combinations)...)
    combination_matrix = combination_matrix'

    Threads.@threads for i = 1:length(lookup(traits_end.H, X))
        for j = 1:length(lookup(traits_end.H, Y))
            if climate.habitable[i, j] == 1
            
                #print(i / length(lookup(traits_start.H, X)))

                # Extract climate at location 
                env = (precip_monthly = parent(climate.precip[i, j, :]), 
                tair_monthly = parent(climate.tair[i, j, :]), 
                tair_annual = mean(parent(climate.tair[i, j, :])), 
                rsds_monthly = parent(climate.rsds[i, j, :]), 
                rss_monthly = parent(climate.rss[i, j, :]), 
                rls_monthly = parent(climate.rls[i, j, :]), 
                daylength = parent(climate.daylength[i, j, :]), 
                CO2_ppm = climate.CO2_ppm, 
                precip_annual = mean(parent(climate.precip[i, j, :]))
                )

                output = Vector()
                for a in 1:length(combination_matrix[:,1])

                    # Initialize test traits 
                    tr = (H = combination_matrix[a, 1], 
                    a_ll = combination_matrix[a, 3], 
                    C_leaf = combination_matrix[a, 2], 
                    seasonality = traits_end.seasonality[i, j], 
                    r_s_r = traits_end.r_s_r[i, j], 
                    Tave_optim = traits_end.Tave_optim[i, j], 
                    Tmax_optim = traits_end.Tmax_optim[i, j], 
                    Tmin_optim = traits_end.Tmin_optim[i, j], 
                    Pave_optim = traits_end.Pave_optim[i, j])

                    tr_complete = plant_allometry(tr = tr, par = pars2)
                    GPP_output = GPP_function_for_optimization(env = env, tr = tr_complete, par = pars2) # For now I use the same procedure for lambda as in optimization - is more consistent! 
                    R_maintenance_output = R_maintenance_function(env = env, tr = tr_complete, par = pars2, GPP_out = GPP_output)
                    NPP_output = calc_NPP(GPP_out = GPP_output, R_maintenance = R_maintenance_output, par = pars2)
                    C_turnover_total_output = C_turnover_function(env = env, tr = tr_complete, par = pars2)
                    Net_C_gain_output = calc_net_C_gain(NPP = NPP_output, C_turnover_total = C_turnover_total_output, par=pars2)

                    out = Net_C_gain_output
                    if NPP_output < 0 || Net_C_gain_output < 0
                        out = 0
                    else 
                        out = 1
                    end

                    push!(output, out)

                end

                functional_diversity_record[i, j] = sum(output)/length(combination_matrix[:,1])

            end
        end
    end


    return(functional_diversity_record)

end


# Species richness sampling based on functional diversity and landscape complexity 
run_TREED_richness_assessment = function(functional_diversity_record, vegetation_output, climate, RI_landscape_window)

    # Get topography
    topo = deepcopy(climate.topo)

    # make habitability mask
    habitable = deepcopy(topo)
    habitable[habitable.<0] .= NaN
    habitable[habitable.>=0] .= 1

    # Define habitat classes 
    # Second attempt only considering biological characteristics: 
    # seasonality = 0-0.5, 0.5-1
    # H = 0-5, 5-10, 10-20, 20-30, 30-40, >40 
    # C_leaf = 0-500, 500-1000, 1000-1500, 1500-2000, >2000
    # a_ll = 0-0.5, 0.5-1, 1-2, 2-3, 3-4, >4
    # unproductive = NPP < 0

    H_classes = [range(0, stop=prevfloat(5.0), length=2), range(5, stop=prevfloat(10.0), length=2), range(10, stop=prevfloat(20.0), length=2), range(20, stop=prevfloat(30.0), length=2), range(30, stop=prevfloat(40.0), length=2), range(40, stop=100, length=2)]
    a_ll_classes = [range(0, stop=prevfloat(0.5), length=2), range(0.5, stop=prevfloat(1.0), length=2), range(1, stop=prevfloat(2.0), length=2), range(2, stop=prevfloat(3.0), length=2), range(3, stop=prevfloat(4.0), length=2), range(4, stop=10, length=2)]
    C_leaf_classes = [range(0, stop=prevfloat(500.0), length=2), range(500, stop=prevfloat(1000.0), length=2), range(1000, stop=prevfloat(2000.0), length=2), range(2000, stop=prevfloat(3000.0), length=2), range(3000, stop=prevfloat(4000.0), length=2), range(4000, stop=10000, length=2)]
    seasonality_classes = [range(0, stop=prevfloat(1.0), length=2), range(1, stop=1, length=2)]

    class_array = collect(Base.product(H_classes, a_ll_classes, C_leaf_classes, seasonality_classes))

    ##### 
    # Loop through every grid cell and assign a habitat class 1 to as many as there are classes
    #####

    habitat_class = deepcopy(habitable)
    habitat_class .= NaN

    for i = 1:length(lookup(vegetation_output.H, X))
        for j = 1:length(lookup(vegetation_output.H, Y))
            if habitable[i, j] == 1

                H_current = vegetation_output.H[i, j]
                a_ll_current = vegetation_output.a_ll[i, j]
                C_leaf_current = vegetation_output.C_leaf[i, j]
                seasonality_current = vegetation_output.seasonality[i, j]

                for c = 1:length(class_array)
                    H_range, a_ll_range, C_leaf_range, seasonality_range = class_array[c]
                    if minimum(H_range) <= H_current <= maximum(H_range) && minimum(a_ll_range) <= a_ll_current <= maximum(a_ll_range) && minimum(C_leaf_range) <= C_leaf_current <= maximum(C_leaf_range) && minimum(seasonality_range) <= seasonality_current <= maximum(seasonality_range)
                        habitat_class[i, j] = c
                    end
                end

            end
        end
    end

    # Set non productive regions to habitat class 0 
    habitat_class[vegetation_output.NPP.<=0] .= 0

    # Set ocean to habitat class -1
    habitat_class[isnan.(habitat_class)] .= -1

    # Convert to integers 
    habitat_class = convert.(Int32, habitat_class)

    ##### 
    # Loop through every grid cell and calculated Diversity metrices
    #####

    window_size = RI_landscape_window # km: how large of a window size should be considered in the metric calculation? 
    gamma_EH_collection = deepcopy(habitable)
    gamma_GI_collection = deepcopy(habitable)


    habitat_class_circular = CircularArray(habitat_class)

    for i = 1:length(lookup(habitable, X))
        for j = 1:length(lookup(habitable, Y))
            if habitable[i, j] == 1

                # lon/lat resolution in °
                longitude_resolution = 360 / length(lookup(habitable, X))
                latitude_resolution = 180 / length(lookup(habitable, Y))

                # Get latitude of current location 
                current_latitude = lookup(habitable, Y)[j]

                # cell to cell longitude distance at current latitude
                cell_to_cell_longitude = (6378 * cos(current_latitude * pi / 180) * pi / 180) * longitude_resolution
                cell_to_cell_latitude = 111.11 * latitude_resolution

                window_width_NS = round(Int, (window_size / 2) / cell_to_cell_latitude)
                window_width_WE = round(Int, (window_size / 2) / cell_to_cell_longitude)

                # Check for maximum one surrounding
                if window_width_WE >= round(Int, length(lookup(habitable, X)) / 2)
                    window_width_WE = round(Int, length(lookup(habitable, X)) / 2)
                end

                # Get habitats within window
                habitats_window = habitat_class_circular[i-window_width_WE:i+window_width_WE, j-window_width_NS:j+window_width_NS]


                # Calculate gamma diversity due to environmental variability 
                # Number of potential different vegetation types = number of cells (if every cell was occupied by a different type)
                gamma_EH = length(unique(habitats_window[habitats_window.!=-1])) ./ length(habitats_window[habitats_window.!=-1])
                gamma_EH_collection[i, j] = gamma_EH


                # Caluculate index of habitat fragmentation 
                habitat_patches = label_components(habitats_window, strel_box((3, 3)), bkg=-1) # background is ocean
                gamma_GI = maximum(habitat_patches) ./ length(habitats_window[habitats_window.!=-1])
                gamma_GI_collection[i, j] = gamma_GI

            end
        end
    end

    #### Save outputs 
    richness_output = (gamma_EH=gamma_EH_collection,
        gamma_GI=gamma_GI_collection, 
        diversity_index = functional_diversity_record .* gamma_EH_collection .* gamma_GI_collection)

    return(richness_output)

end
