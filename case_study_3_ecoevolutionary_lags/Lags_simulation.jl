###########################################################
########## TREED eco-evolutionary lags ####################
###########################################################

# Author: Julian Rogger
# Startdate: 04.2025
# Purpose: Run TREED across the PETM climatic shifts considering different evolution and disperal rates.

# Packages necessary to handle inputs and analyse output
include("../src/TREED.jl")
using .TREED
using Rasters, ArchGDAL, NCDatasets
using DimensionalData
using DimensionalData.Lookups
using Plots
using DataFrames
using GMT


###########################################################
########## Run TREED model 

# Get climate data of the present in right format 
# Derived from Christine Shields, NCAR and Vera Korasidis, University of Melbourne: Korasidis, V. A., Wing, S. L., Shields, C. A., & Kiehl, J. T. (2022). Global changes in terrestrial vegetation and continental climate during the Paleocene-Eocene Thermal Maximum. https://doi.org/10.1029/2021PA004325
# Information about argument format required: ?TREEDnonsteadycontinuous
tair_pre = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/pre_PETM_tair_0.25.nc") # in K
precip_pre = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/pre_PETM_precip_0.25.nc") # in m/s
clt_pre = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/pre_PETM_clt_pseudo_0.25.nc")
rsds_pre = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/pre_PETM_rsds_0.25.nc")

tair_peak = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/peak_PETM_tair_0.25.nc") # in K
precip_peak = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/peak_PETM_precip_0.25.nc") # in m/s
clt_peak = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/peak_PETM_clt_pseudo_0.25.nc")
rsds_peak = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/peak_PETM_rsds_0.25.nc")

topo = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/PETM_topo_0.25.nc")


# Construct landscape/environment to which simulation will be applied
tairvec = [tair_pre, tair_peak]
precipvec = [precip_pre, precip_peak]
cltvec = [clt_pre, clt_peak]
rsdsvec = [rsds_pre, rsds_peak]
CO2vec = [680, 1590]
topovec = [topo, topo]

# Additional arguments needed:
res = 2. # Target resolution
evorate = 1. # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 700. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_immediate_evolution"

# Run TREED for the two PETM states, assuming quasi immediate adaptation 
TREED_output = TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)






###########################################################
########## PETM experiments for step change in CO2

# Experiment 1: slow evolution, fast dispersal
# Experiment 2: fast evolution, slow dispersal
# Experiment 3: Intermediate evolution, intermediate dispersal 

# Timesteps = 20, CO2 change after 5 
# Construct landscape/environment to which simulation will be applied
tairvec = vcat(repeat([tair_pre], 5), repeat([tair_peak], 15));
precipvec = vcat(repeat([precip_pre], 5), repeat([precip_peak], 15));
cltvec = vcat(repeat([clt_pre], 5), repeat([clt_peak], 15));
rsdsvec = vcat(repeat([rsds_pre], 5), repeat([rsds_peak], 15));
CO2vec = vcat(repeat([680], 5), repeat([1590], 15));
topovec = repeat([topo], 20);


############
# Experiment 1: slow evolution, fast dispersal
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp1 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp1, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot(NPP_exp1, label="Experiment 1", xlab="Timestep", ylab="NPP (Pg C per year)")


############
# Experiment 2: fast evolution, slow dispersal
res = 4. # Target resolution
evorate = 0.75 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 200. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_fast_evo_slow_dispersal"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp2 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_fast_evo_slow_dispersal/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp2, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp2, label="Experiment 2")



############
# Experiment 3: intermediate evolution, intermediate dispersal
res = 4. # Target resolution
evorate = 0.1 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 400. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp3 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp3, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp3, label="Experiment 3")







###########################################################
########## temperature niche breadth sensitivity 

# Experiment 4: niche breadth 0
# Experiment 5: niche breadth 0.02
# Experiment 6: niche breadth 0.04
# Experiment 7: niche breadth 0.06
# Experiment 8: niche breadth 0.08

############
# Experiment 4: nichebreadth = 0
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir, temp_niche_breadth_parameter = 0.0)

NPP_exp4 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp4, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot(NPP_exp4, label="Experiment 4")


############
# Experiment 5: nichebreadth = 0.02 (default)
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0.02"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir, temp_niche_breadth_parameter = 0.02)

NPP_exp5 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0.02/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp5, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp5, label="Experiment 5")

############
# Experiment 6: nichebreadth = 0.04
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0.04"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir, temp_niche_breadth_parameter = 0.04)

NPP_exp6 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0.04/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp6, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp6, label="Experiment 6")


############
# Experiment 7: nichebreadth = 0.06
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0.06"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir, temp_niche_breadth_parameter = 0.06)

NPP_exp7 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0.06/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp7, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp7, label="Experiment 7")


############
# Experiment 8: nichebreadth = 0.08
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0.08"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir, temp_niche_breadth_parameter = 0.08)

NPP_exp8 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_interm_dispersal_k0.08/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp8, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp8, label="Experiment 8")







###########################################################
########## dispersal sensitivity 

# Experiment 9: dispersal radius = 0 km
# Experiment 10: dispersal radius = 200 km
# Experiment 11: dispersal radius = 400 km
# Experiment 12: dispersal radius = 600 km
# Experiment 13: dispersal radius = 800 km



############
# Experiment 9: dispersal radius = 0 km
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 0. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_0"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp9 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_0/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp9, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot(NPP_exp9, label="Experiment 9")



############
# Experiment 10: dispersal radius = 200 km
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 200. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_200"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp10 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_200/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp10, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp10, label="Experiment 10")


############
# Experiment 11: dispersal radius = 400 km
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 400. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_400"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp11 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_400/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp11, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp11, label="Experiment 11")


############
# Experiment 12: dispersal radius = 600 km
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_600"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp12 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_600/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp12, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp12, label="Experiment 12")



############
# Experiment 13: dispersal radius = 800 km
res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 800. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_800"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir, temp_niche_breadth_parameter = 0.03)

NPP_exp13 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_interm_evo_dispersal_800/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp13, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp13, label="Experiment 13")





###########################################################
########## warming scenarios

# Experiment 14: land surface warming ~2°C
# Experiment 15: land surface warming ~4°C
# Experiment 16: land surface warming ~6°C
# Experiment 17: land surface warming ~8°C

tair_pre_orig = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/pre_PETM_tair_0.25.nc") # in K
precip_pre_orig = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/pre_PETM_precip_0.25.nc") # in m/s
clt_pre_orig = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/pre_PETM_clt_pseudo_0.25.nc")
rsds_pre_orig = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/pre_PETM_rsds_0.25.nc")

tair_peak_orig = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/peak_PETM_tair_0.25.nc") # in K
precip_peak_orig = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/peak_PETM_precip_0.25.nc") # in m/s
clt_peak_orig = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/peak_PETM_clt_pseudo_0.25.nc")
rsds_peak_orig = Raster("./case_study_3_ecoevolutionary_lags/PETM_climate_topo_inputs/peak_PETM_rsds_0.25.nc")

CO2_orig = [680, 1590]



# Experiment 14: slow evolution, fast dispersal
# Linear interpolation to new peak CO2 level 
CO2_2deg = [680, 907.5]
tair_2deg = (1 - ((907.5 - 680)/(1590 - 680))) .* tair_pre_orig .+ (1 - ((1590 - 907.5)/(1590 - 680))) .* tair_peak_orig
precip_2deg = (1 - ((907.5 - 680)/(1590 - 680))) .* precip_pre_orig .+ (1 - ((1590 - 907.5)/(1590 - 680))) .* precip_peak_orig
clt_2deg = (1 - ((907.5 - 680)/(1590 - 680))) .* clt_pre_orig .+ (1 - ((1590 - 907.5)/(1590 - 680))) .* clt_peak_orig
rsds_2deg = (1 - ((907.5 - 680)/(1590 - 680))) .* rsds_pre_orig .+ (1 - ((1590 - 907.5)/(1590 - 680))) .* rsds_peak_orig


tairvec = vcat(repeat([tair_pre], 5), repeat([tair_2deg], 15));
precipvec = vcat(repeat([precip_pre], 5), repeat([precip_2deg], 15));
cltvec = vcat(repeat([clt_pre], 5), repeat([clt_2deg], 15));
rsdsvec = vcat(repeat([rsds_pre], 5), repeat([rsds_2deg], 15));
CO2vec = vcat(repeat([680], 5), repeat([907.5], 15));
topovec = repeat([topo], 20);

res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal_warming2deg"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp14 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal_warming2deg/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp14, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot(NPP_exp14, label="Experiment 14", xlab="Timestep", ylab="NPP (Pg C per year)")



# Experiment 15: slow evolution, fast dispersal
# Linear interpolation to new peak CO2 level 
CO2_4deg = [680, 1135]
tair_2deg = (1 - ((1135 - 680)/(1590 - 680))) .* tair_pre_orig .+ (1 - ((1590 - 1135)/(1590 - 680))) .* tair_peak_orig
precip_2deg = (1 - ((1135 - 680)/(1590 - 680))) .* precip_pre_orig .+ (1 - ((1590 - 1135)/(1590 - 680))) .* precip_peak_orig
clt_2deg = (1 - ((1135 - 680)/(1590 - 680))) .* clt_pre_orig .+ (1 - ((1590 - 1135)/(1590 - 680))) .* clt_peak_orig
rsds_2deg = (1 - ((1135 - 680)/(1590 - 680))) .* rsds_pre_orig .+ (1 - ((1590 - 1135)/(1590 - 680))) .* rsds_peak_orig


tairvec = vcat(repeat([tair_pre], 5), repeat([tair_2deg], 15));
precipvec = vcat(repeat([precip_pre], 5), repeat([precip_2deg], 15));
cltvec = vcat(repeat([clt_pre], 5), repeat([clt_2deg], 15));
rsdsvec = vcat(repeat([rsds_pre], 5), repeat([rsds_2deg], 15));
CO2vec = vcat(repeat([680], 5), repeat([1135], 15));
topovec = repeat([topo], 20);

res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal_warming4deg"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp15 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal_warming4deg/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp15, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp15, label="Experiment 15", xlab="Timestep", ylab="NPP (Pg C per year)")



# Experiment 16: slow evolution, fast dispersal
# Linear interpolation to new peak CO2 level 
CO2_6deg = [680, 1362.5]
tair_2deg = (1 - ((1362.5 - 680)/(1590 - 680))) .* tair_pre_orig .+ (1 - ((1590 - 1362.5)/(1590 - 680))) .* tair_peak_orig
precip_2deg = (1 - ((1362.5 - 680)/(1590 - 680))) .* precip_pre_orig .+ (1 - ((1590 - 1362.5)/(1590 - 680))) .* precip_peak_orig
clt_2deg = (1 - ((1362.5 - 680)/(1590 - 680))) .* clt_pre_orig .+ (1 - ((1590 - 1362.5)/(1590 - 680))) .* clt_peak_orig
rsds_2deg = (1 - ((1362.5 - 680)/(1590 - 680))) .* rsds_pre_orig .+ (1 - ((1590 - 1362.5)/(1590 - 680))) .* rsds_peak_orig


tairvec = vcat(repeat([tair_pre], 5), repeat([tair_2deg], 15));
precipvec = vcat(repeat([precip_pre], 5), repeat([precip_2deg], 15));
cltvec = vcat(repeat([clt_pre], 5), repeat([clt_2deg], 15));
rsdsvec = vcat(repeat([rsds_pre], 5), repeat([rsds_2deg], 15));
CO2vec = vcat(repeat([680], 5), repeat([1362.5], 15));
topovec = repeat([topo], 20);

res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal_warming6deg"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp16 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal_warming6deg/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp16, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp16, label="Experiment 16", xlab="Timestep", ylab="NPP (Pg C per year)")


# Experiment 17: slow evolution, fast dispersal
# Linear interpolation to new peak CO2 level 
CO2_8deg = [680, 1590]
tair_2deg = (1 - ((1590 - 680)/(1590 - 680))) .* tair_pre_orig .+ (1 - ((1590 - 1590)/(1590 - 680))) .* tair_peak_orig
precip_2deg = (1 - ((1590 - 680)/(1590 - 680))) .* precip_pre_orig .+ (1 - ((1590 - 1590)/(1590 - 680))) .* precip_peak_orig
clt_2deg = (1 - ((1590 - 680)/(1590 - 680))) .* clt_pre_orig .+ (1 - ((1590 - 1590)/(1590 - 680))) .* clt_peak_orig
rsds_2deg = (1 - ((1590 - 680)/(1590 - 680))) .* rsds_pre_orig .+ (1 - ((1590 - 1590)/(1590 - 680))) .* rsds_peak_orig


tairvec = vcat(repeat([tair_pre], 5), repeat([tair_2deg], 15));
precipvec = vcat(repeat([precip_pre], 5), repeat([precip_2deg], 15));
cltvec = vcat(repeat([clt_pre], 5), repeat([clt_2deg], 15));
rsdsvec = vcat(repeat([rsds_pre], 5), repeat([rsds_2deg], 15));
CO2vec = vcat(repeat([680], 5), repeat([1590], 15));
topovec = repeat([topo], 20);

res = 4. # Target resolution
evorate = 0.01 # Rate of trait evolution from current to optimum 
startinsteady = true # Initialize the model at the first time step in optimum state 
dispersal = 600. # Dispersal rate in km (radius)
outputdir = "./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal_warming8deg"

# Run TREED 
TREED_output = TREED.TREEDnonsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, evorate=evorate, startinsteady=startinsteady, dispersal=dispersal, res=res, outputdir=outputdir)

NPP_exp17 = Vector()
for t=1:20
    out = RasterStack(string("./case_study_3_ecoevolutionary_lags/TREED_PETM_lags_output_slow_evo_fast_dispersal_warming8deg/TREED_output_timestep_",t,".nc"))
    push!(NPP_exp17, sum(replace!(out.NPP .* TREED.raster_area(out.NPP) .* 1e+6, NaN => 0)) ./ 1e+15)
end
Plots.plot!(NPP_exp17, label="Experiment 17", xlab="Timestep", ylab="NPP (Pg C per year)")



# Combined plot of all sensitivity simulations 
data_base_scenarios = DataFrame(time = 1:20, NPP1 = NPP_exp1, NPP2 = NPP_exp2, NPP3 = NPP_exp3)
data_base_scenarios = convert.(Float64, data_base_scenarios)

GMT.basemap(projection=:linear, region=(1, 20, 35, 115), figsize=(10, 8), theme=("A2xy"), 
xlabel="", ylabel="Total NPP (Pg C year@+-1@+)", par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=7))
GMT.plot!(data_base_scenarios.time, data_base_scenarios.NPP1, pen=(4,:darkred), legend=(label="Slow evolution, fast dispersal", pos=:TL, box=:none,))
GMT.plot!(data_base_scenarios.time, data_base_scenarios.NPP2, pen=(4,:darkblue), legend=(label="Fast evolution, slow dispersal", pos=:TL, box=:none,))
GMT.plot!(data_base_scenarios.time, data_base_scenarios.NPP3, pen=(4,:black), legend=(label="Intermediate evolution and dispersal", pos=:TL, box=:none,))


data_nichebreadth = DataFrame(time = 1:20, NPP4 = NPP_exp4, NPP5 = NPP_exp5, NPP6 = NPP_exp6, NPP7 = NPP_exp7, NPP8 = NPP_exp8)
data_nichebreadth = convert.(Float64, data_nichebreadth)

GMT.basemap!(projection=:linear, region=(1, 20, 10, 105), figsize=(10, 8), theme=("A2xy"), 
xlabel="", ylabel=" ", par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=7), xshift=11)
GMT.plot!(data_nichebreadth.time, data_nichebreadth.NPP4, pen=(4,"253/231/37"), legend=(label="nichebreadth k=0", pos=:BR, box=:none,))
GMT.plot!(data_nichebreadth.time, data_nichebreadth.NPP5, pen=(4,"94/201/98"), legend=(label="nichebreadth k=0.02", pos=:BR, box=:none,))
GMT.plot!(data_nichebreadth.time, data_nichebreadth.NPP6, pen=(4,"33/145/140"), legend=(label="nichebreadth k=0.04", pos=:BR, box=:none,))
GMT.plot!(data_nichebreadth.time, data_nichebreadth.NPP7, pen=(4, "59/82/139"), legend=(label="nichebreadth k=0.06", pos=:BR, box=:none,))
GMT.plot!(data_nichebreadth.time, data_nichebreadth.NPP8, pen=(4, "68/1/84"), legend=(label="nichebreadth k=0.08", pos=:BR, box=:none,))
text!("@~a@~ = 0.01, dispersal = 600 km",frame=:none,region=(0,10,0,10), proj=:linear, x=2.5, y=0.5, noclip=true ,font=(9,"Helvetica",:black)) 


data_dispersal = DataFrame(time = 1:20, NPP9 = NPP_exp9, NPP10 = NPP_exp10, NPP11 = NPP_exp11, NPP12 = NPP_exp12, NPP13 = NPP_exp13)
data_dispersal = convert.(Float64, data_dispersal)

GMT.basemap!(projection=:linear, region=(1, 20, 0, 105), figsize=(10, 8), theme=("A2xy"), 
xlabel="Timestep", ylabel="Total NPP (Pg C year@+-1@+)", par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=7), xshift=-11, yshift=-9)
GMT.plot!(data_dispersal.time, data_dispersal.NPP9, pen=(4,"252/253/191"), legend=(label="dispersal radius=0 km", pos=:RM, box=:none,))
GMT.plot!(data_dispersal.time, data_dispersal.NPP10, pen=(4,"252/137/97"), legend=(label="dispersal radius=200 km", pos=:RM, box=:none,))
GMT.plot!(data_dispersal.time, data_dispersal.NPP11, pen=(4,"183/55/121"), legend=(label="dispersal radius=400 km", pos=:RM, box=:none,))
GMT.plot!(data_dispersal.time, data_dispersal.NPP12, pen=(4, "81/18/124"), legend=(label="dispersal radius=600 km", pos=:RM, box=:none,))
GMT.plot!(data_dispersal.time, data_dispersal.NPP13, pen=(4, "0/0/4"), legend=(label="dispersal radius=800 km", pos=:RM, box=:none,))
text!("@~a@~ = 0.01, k = 0.02",frame=:none,region=(0,10,0,10), proj=:linear, x=1.5, y=9.5, noclip=true ,font=(9,"Helvetica",:black)) 


data_warming = DataFrame(time = 1:20, NPP14 = NPP_exp14, NPP15 = NPP_exp15, NPP16 = NPP_exp16, NPP17 = NPP_exp17)
data_warming = convert.(Float64, data_warming)

GMT.basemap!(projection=:linear, region=(1, 20, 40, 105), figsize=(10, 8), theme=("A2xy"), 
xlabel="Timestep", ylabel=" ", par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=7), xshift=11)
GMT.plot!(data_warming.time, data_warming.NPP14, pen=(4,"254/230/206"), legend=(label="Land warming ~ 2@~\\260@~C", pos=:BR, box=:none,))
GMT.plot!(data_warming.time, data_warming.NPP15, pen=(4,"253/174/107"), legend=(label="Land warming ~ 4@~\\260@~C", pos=:BR, box=:none,))
GMT.plot!(data_warming.time, data_warming.NPP16, pen=(4,"230/85/13"), legend=(label="Land warming ~ 6@~\\260@~C", pos=:BR, box=:none,))
GMT.plot!(data_warming.time, data_warming.NPP17, pen=(4, "166/54/3"), legend=(label="Land warming ~ 8@~\\260@~C", pos=:BR, box=:none,))
text!("@~a@~ = 0.01, k = 0.02, dispersal = 600 km",frame=:none,region=(0,10,0,10), proj=:linear, x=3, y=9.5, noclip=true ,font=(9,"Helvetica",:black), 
dpi=700, name="./case_study_3_ecoevolutionary_lags/plots/evolutionary_lags_sensitivities.png")
























# Combine with steady state plots from case study 2
convert_raster_to_GMT_grid = function(raster)
    increment = step(lookup(raster, X))
    raster = replace_missing(raster, NaN)
    data_xyz = DataFrame(raster)
    matrix_xyz = Matrix(data_xyz)
    GMT_grd = xyz2grd(matrix_xyz, limits=(minimum(matrix_xyz[:,1]), maximum(matrix_xyz[:,1]), minimum(matrix_xyz[:,2]), maximum(matrix_xyz[:,2])), inc=increment)
    return(GMT_grd)
end

TREED_output_pre = RasterStack("./case_study_2_PETM/TREED_PETM_output/TREED_output_timestep_1.nc")
TREED_output_peak = RasterStack("./case_study_2_PETM/TREED_PETM_output/TREED_output_timestep_2.nc")

# Steady state vegetation height 
H_cpt = makecpt(cmap=:bamako, range=(0, 50), inverse=true, overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="0/59/71", COLOR_FOREGROUND="255/229/172"))
grdimage(convert_raster_to_GMT_grid(TREED_output_pre.H), projection=:Mollweide, theme="A2xy",
    cmap=H_cpt, xaxis=(annot=0,), yaxis=(annot=60,), figsize=10, par=(FONT_ANNOT=7,))
grdimage!(convert_raster_to_GMT_grid(TREED_output_peak.H), projection=:Mollweide, theme="A2xy",
    cmap=H_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=10, par=(FONT_ANNOT=7,), xshift=10.25)
colorbar!(pos=(achor=:RM,), frame=(annot=:auto, ticks=:auto, xlabel="H (m)"),par=(FONT_ANNOT_PRIMARY=12,))
text!("(a)",frame=:none,region=(0,10,0,10), xshift=-10.25, proj=:linear, x=0, y=5.0, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(b)",frame=:none,region=(0,10,0,10), proj=:linear, x=6.75, y=5.0, noclip=true ,font=(10,"Helvetica",:black)) 

NPP_cpt = makecpt(cmap=:plasma, range=(0, 1300), overrule_bg=true, par=(COLOR_NAN=230, COLOR_BACKGROUND="15/7/136", COLOR_FOREGROUND="240/248/35"))
grdimage!(convert_raster_to_GMT_grid(TREED_output_pre.NPP), projection=:Mollweide, theme="A2xy",
    cmap=NPP_cpt, xaxis=(annot=0,), yaxis=(annot=60,), figsize=10, par=(FONT_ANNOT=7,), yshift=-6.25)
grdimage!(convert_raster_to_GMT_grid(TREED_output_peak.NPP), projection=:Mollweide, theme="A2xy",
    cmap=NPP_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=10, par=(FONT_ANNOT=7,), xshift=10.25)
colorbar!(pos=(achor=:RM,), frame=(annot=:auto, ticks=:auto, xlabel="NPP (g C m@+-2@+)"),par=(FONT_ANNOT_PRIMARY=12,))
text!("(c)",frame=:none,region=(0,10,0,10), xshift=-10.25, proj=:linear, x=0, y=5.0, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(d)",frame=:none,region=(0,10,0,10), proj=:linear, x=6.75, y=5.0, noclip=true ,font=(10,"Helvetica",:black)) 

GMT.basemap!(projection=:linear, region=(1, 20, 40, 115), figsize=(12, 7), theme=("A2xy"), 
xlabel="Timestep", ylabel="Total NPP (Pg C year@+-1@+)", par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=7), yshift=-9, xshift=2) # 
GMT.plot!(data.time, data.NPP1, pen=(4,:darkred), legend=(label="Slow evolution, fast dispersal", pos=:TL, box=:none,))
GMT.plot!(data.time, data.NPP2, pen=(4,:darkblue), legend=(label="Fast evolution, slow dispersal", pos=:TL, box=:none,))
GMT.plot!(data.time, data.NPP3, pen=(4,:black), legend=(label="Intermediate evolution and dispersal", pos=:TL, box=:none,))
text!("(e)", x=-1, y=110, noclip=true, font=(10,"Helvetica",:black), 
dpi=700, name="./case_study_3_ecoevolutionary_lags/plots/evolutionary_lags.png")
