# __TREED__ 
### Trait Ecology and Evolution over Deep time vegetation model 

This repository contains the code of the TREED vegetation model, a trait- and optimality-based vegetation model, designed for eco-evolutionary vegetation modelling over geologic time.  

To run the model, download this repository and open julia from the repository. 

To enable multithreading, start the julia session enabling multiple threads, for example by: 

```bash
$ julia --threads 12
```

(Or edit your "julia.NumThreads" option in your settings.json file in Visual Studio Code)

Next, activate the project environment

```julia
] activate .
```

Next, run 

```julia
] instantiate 
```

which will install all the packages listed in the Project.toml and the versions specified in the Manifest.toml. You should now be good to go. 

### Source code description

You can find the model code in the src directory. 

 - TREED_physiological_functions.jl: contains all the plant physiological functions from photosynthesis and respiration to carbon turnover. 
- TREED_model_functions.jl: contains the functions applied to the global grid at each time step, including the calculation and optimization of photosynthetic fluxes, evolution and dispersal, as well as the end of timestep carbon flux and storage calculation 

- TREED_pars.jl: contains the default values of fixed paramters. 

- TREEDsteadystep.jl: function definition to run a one time step, steady state simulation.

- TREEDsteadycontinuous.jl: function definition to run the model multiple time steps, always in steady state. 

- TREEDnonsteadycontinuous.jl: function definition to run the model over multiple time steps, not in steady state and considering lagged eco-evolutionary dynamics (considering limited rates of trait evolution and dispersal between time steps).

- TREEDnonsteadystep.jl: function definition to run one time step in eco-evolutionary mode, restarting from a previous trait distribution.


### Case studies: 

This repostiroy comes with three case study applications of the model, which can be used to get started with the model, or can be modified according to the user's interest. (A link to access the climate and topography input data will be provided soon..)

- Case study 1: Run the model one time step in steady state using present-day climate inputs for validation of simulated vegetation structures and fluxes.

- Case study 2: Run the model two time steps in steady state for the Paleocene--Eocene Thermal Maximum (PETM) using pre-PETM and peak-PETM climate conditions. 

- Case study 3: Run the model 20 time steps across the PETM and consider lagged eco-evolutionary dynamics across the climate warming in the form of limited trait adaptation and dispersal. Further, several examples of how to alter default parameters of the model using the TREEDnonsteadycontinuous.jl function. 

Essentially, one time step of the model can be run with the following lines of code: 


```julia
include("../src/TREED.jl")
using .TREED
using Rasters, ArchGDAL, NCDatasets
using DimensionalData
using DimensionalData.Lookups

# read climate input rasters 
tair = Raster("./your_monthly_tair.nc")
precip = Raster("./your_monthly_precip.nc") 
clt= Rasster("./your_monthly_cloud_cover.nc")
rsds = Raster("./your_monthly_downwelling_radiation.nc")
topo = Raster("./time_step_topography.nc")

# Additional arguments needed:
res = 0.5 # Target resolution
CO2 = 360.0 # Current atmospheric CO2, transferred to vegetation model
FDsampling = true # Assessment of functional diversity  
RIsampling = true # Assessment of species richness potential ("diversity index")
RI_landscape_window = 300.0 # Width of the landscape window used for the diversity assessment (in km)
outputdir = "./case_study_1_present_day_validation/TREED_present_output" # where should outputs be stored

# Run TREED model
TREED_output = TREEDsteadystep(tair=tair, precip=precip, clt=clt, rsds=rsds, topo=topo, CO2=CO2, res=res, FDsampling=FDsampling, RIsampling=RIsampling, RI_landscape_window=RI_landscape_window, outputdir=outputdir)
```



### Basic usage: 

To run one time step of the model, you need topographic and climatic input data. This includes a raster of monthly average surface air temperatures in K, monthly average precipitation rates in m/s, monthly average cloud cover data as 0-1 fraction, monthly average surface downwelling shortwave radiation in W/m2. Further, you need a topography raster required to define the land-sea mask and a atmospheric CO2 concentration value. 

When defining your inputs, make sure that all rasters are in the same orientation as the model will just iterate over the grid and not check that each grid index is the same latitude and longitude in each of the input rasters. 

For more info, please also see ?TREEDsteadystep and the case studies. 


### Credits: 


This model depends strongly on the Raster.jl package (https://rafaqz.github.io/Rasters.jl/stable/), the DimensionalData.jl package (https://rafaqz.github.io/DimensionalData.jl/stable/). 

Also credits to the GMT.jl (https://www.generic-mapping-tools.org/GMTjl_doc/) package used for plotting of the case study outputs. 


### License: 

This model is published under the Apache License, Version 2.0. 






