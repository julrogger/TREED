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