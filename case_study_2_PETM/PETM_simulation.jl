###########################################################
########## TREED run for the PETM #########################
###########################################################

# Author: Julian Rogger
# Startdate: 04.2025
# Purpose: Run TREED using PETM climatic and topographic inputs.

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
########## Run TREED model for pre PETM conditions

# Get climate data of the present in right format 
# Derived from Christine Shields, NCAR and Vera Korasidis, University of Melbourne: Korasidis, V. A., Wing, S. L., Shields, C. A., & Kiehl, J. T. (2022). Global changes in terrestrial vegetation and continental climate during the Paleocene-Eocene Thermal Maximum. https://doi.org/10.1029/2021PA004325
# Information about argument format required: ?TREEDsteadystep
tair = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/pre_PETM_tair_0.25.nc") # in K
precip = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/pre_PETM_precip_0.25.nc") # in m/s
clt = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/pre_PETM_clt_pseudo_0.25.nc")
rsds = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/pre_PETM_rsds_0.25.nc")
topo = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/PETM_topo_0.25.nc")

# Additional arguments needed:
res = 0.5 # Target resolution
CO2 = 600.0 # Current atmospheric CO2, transferred to vegetation model
FDsampling = false # Assessment of functional diversity  
RIsampling = false # Assessment of species richness potential ("diversity index")
RI_landscape_window = 300.0 # Width of the landscape window used for the diversity assessment (in km)
outputdir = "./case_study_2_PETM/TREED_PETM_output"

# Run TREED 
TREED_output = TREEDsteadystep(tair=tair, precip=precip, clt=clt, rsds=rsds, topo=topo, CO2=CO2, res=res, FDsampling=FDsampling, RIsampling=RIsampling, RI_landscape_window=RI_landscape_window, outputdir=outputdir)


#########################################################################
########## Run TREED model for pre + peak PETM conditions to steady state 

# Get climate data of the present in right format 
# Derived from Christine Shields, NCAR and Vera Korasidis, University of Melbourne: Korasidis, V. A., Wing, S. L., Shields, C. A., & Kiehl, J. T. (2022). Global changes in terrestrial vegetation and continental climate during the Paleocene-Eocene Thermal Maximum. https://doi.org/10.1029/2021PA004325
# Information about argument format required: ?TREEDsteadystep
tair_pre = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/pre_PETM_tair_0.25.nc") # in K
precip_pre = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/pre_PETM_precip_0.25.nc") # in m/s
clt_pre = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/pre_PETM_clt_pseudo_0.25.nc")
rsds_pre = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/pre_PETM_rsds_0.25.nc")

tair_peak = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/peak_PETM_tair_0.25.nc") # in K
precip_peak = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/peak_PETM_precip_0.25.nc") # in m/s
clt_peak = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/peak_PETM_clt_pseudo_0.25.nc")
rsds_peak = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/peak_PETM_rsds_0.25.nc")

topo = Raster("./case_study_2_PETM/PETM_climate_topo_inputs/PETM_topo_0.25.nc")

# Construct landscape/environment to which simulation will be applied
tairvec = [tair_pre, tair_peak]
precipvec = [precip_pre, precip_peak]
cltvec = [clt_pre, clt_peak]
rsdsvec = [rsds_pre, rsds_peak]
CO2vec = [680, 1590]
topovec = [topo, topo]


# Additional arguments needed:
res = 0.5 # Target resolution
FDsampling = false # Assessment of functional diversity  
RIsampling = false # Assessment of species richness potential ("diversity index")
RI_landscape_window = 300.0 # Width of the landscape window used for the diversity assessment (in km)
outputdir = "./case_study_2_PETM/TREED_PETM_output"

# Run TREED for both timesteps 
TREEDsteadycontinuous(tairvec=tairvec, precipvec=precipvec, cltvec=cltvec, rsdsvec=rsdsvec, topovec=topovec, CO2vec=CO2vec, res=res, FDsampling=FDsampling, RIsampling=RIsampling, RI_landscape_window=RI_landscape_window, outputdir=outputdir)

TREED_output_pre = RasterStack("./case_study_2_PETM/TREED_PETM_output/TREED_output_timestep_1.nc")
TREED_output_peak = RasterStack("./case_study_2_PETM/TREED_PETM_output/TREED_output_timestep_2.nc")

###########################################################
########## Plotting (with GMT.jl)

# Helper function to make transfer from Raster to GMT grid
convert_raster_to_GMT_grid = function(raster)
    increment = step(lookup(raster, X))
    raster = replace_missing(raster, NaN)
    data_xyz = DataFrame(raster)
    matrix_xyz = Matrix(data_xyz)
    GMT_grd = xyz2grd(matrix_xyz, limits=(minimum(matrix_xyz[:,1]), maximum(matrix_xyz[:,1]), minimum(matrix_xyz[:,2]), maximum(matrix_xyz[:,2])), inc=increment)
    return(GMT_grd)
end


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
text!("(d)",frame=:none,region=(0,10,0,10), proj=:linear, x=6.75, y=5.0, noclip=true ,font=(10,"Helvetica",:black), 
dpi=700, name="./case_study_2_PETM/plots/PETM_pre_peak_steady_state.png") 

