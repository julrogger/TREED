###########################################################
########## TREED helper functions  ########################
###########################################################

# Get raster of area of each grid cell in km2
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

# Calculate an area weighted average of a raster
area_weighted_average = function(raster)
    weights = raster_area(raster)
    mask = .!isnan.(raster)
    area_weighted_average = sum(raster[mask] .* weights[mask]) / sum(weights[mask])
    return(area_weighted_average)
end

### Daylength calculation 
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


### Get env, tr, par for debugging from output
### Afterwards, it is possible to run physiological functions: 

# e.g.: sample_out = get_env_tr_par(output=TREED_output, tair=tair, precip=precip, clt=clt, rsds=rsds, topo=topo, CO2=360.0, res=0.5, lon=-100, lat=40)
# tr = TREED.plant_allometry(tr=sample_out.tr, par=sample_out.pars)
# @enter TREED.GPP_function_for_optimization(env=sample_out.env, tr=tr, par=sample_out.pars)

get_env_tr_par = function(; output, tair, precip, clt, rsds, topo, CO2, res, lon, lat)

  climate = create_TREED_climate_input(tair, precip, clt, rsds, topo, CO2, res)

  # Extract current traits at location 
  tr = (H=output.H[Near(lon), Near(lat)],
    a_ll=output.a_ll[Near(lon), Near(lat)],
    C_leaf=output.C_leaf[Near(lon), Near(lat)],
    seasonality=output.seasonality[Near(lon), Near(lat)],
    r_s_r=output.r_s_r[Near(lon), Near(lat)],
    Tave_optim=output.Tave_optim[Near(lon), Near(lat)],
    Tmax_optim=output.Tmax_optim[Near(lon), Near(lat)],
    Tmin_optim=output.Tmin_optim[Near(lon), Near(lat)],
    Pave_optim=0
  )

  # Extract climate at location 
  env = (precip_monthly=parent(climate.precip[Near(lon), Near(lat), :]),
    tair_monthly=parent(climate.tair[Near(lon), Near(lat), :]),
    tair_annual=mean(parent(climate.tair[Near(lon), Near(lat), :])),
    rsds_monthly=parent(climate.rsds[Near(lon), Near(lat), :]),
    rss_monthly=parent(climate.rss[Near(lon), Near(lat), :]),
    rls_monthly=parent(climate.rls[Near(lon), Near(lat), :]),
    daylength=parent(climate.daylength[Near(lon), Near(lat), :]),
    CO2_ppm=climate.CO2_ppm,
    precip_annual=mean(parent(climate.precip[Near(lon), Near(lat), :]))
  )

  sample_out = (
  tr = tr,
  env = env,
  pars = pars
  )


end