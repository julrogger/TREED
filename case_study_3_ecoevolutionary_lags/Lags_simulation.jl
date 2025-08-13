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


# Plot NPP trajectories under the different scenarios 
data = DataFrame(time = 1:20, NPP1 = NPP_exp1, NPP2 = NPP_exp2, NPP3 = NPP_exp3)
data = convert.(Float64, data)

GMT.basemap(projection=:linear, region=(1, 20, 40, 115), figsize=(12, 8), theme=("A2xy"), 
xlabel="Timestep", ylabel="Total NPP (Pg C year@+-1@+)", par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=7))
GMT.plot!(data.time, data.NPP1, pen=(4,:darkred), legend=(label="Slow evolution, fast dispersal", pos=:TL, box=:none,))
GMT.plot!(data.time, data.NPP2, pen=(4,:darkblue), legend=(label="Fast evolution, slow dispersal", pos=:TL, box=:none,))
GMT.plot!(data.time, data.NPP3, pen=(4,:black), legend=(label="Intermediate evolution and dispersal", pos=:TL, box=:none,),
name="./test.png")




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
