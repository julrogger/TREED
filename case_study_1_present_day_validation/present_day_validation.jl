###########################################################
########## TREED validation present #######################
###########################################################

# Author: Julian Rogger
# Startdate: 04.2025
# Purpose: Run TREED using present day climatic and topographic inputs and validate the model with observational data.

# Packages necessary to handle inputs and analyse output
include("../src/TREED.jl")
using .TREED
using Rasters, ArchGDAL, NCDatasets
using DimensionalData
using DimensionalData.Lookups
using Plots
using DataFrames
using CSV
using Statistics
using GLM
using GMT

###########################################################
########## Run TREED model 

# Get climate data of the present in right format 
# From CHELSA: https://chelsa-climate.org/; Karger et al. (2017): Climatologies at high resolution for the Earth land surface areas. Scientific Data. 4 170122. https://doi.org/10.1038/sdata.2017.122
# Information about argument format required: ?TREEDsteadystep
tair = Raster("./case_study_1_present_day_validation/present_day_climate_topo_inputs/monmean_tas_climatology_1981-2010.nc") .+ 273.15 # in K
precip = Raster("./case_study_1_present_day_validation/present_day_climate_topo_inputs/monmean_pr_climatology_1981-2010.nc") 
month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
for m in 1:12
    precip[:,:,m] .= (precip[:,:,m] ./ 1000) ./ (month_days[m] * 24 * 60 * 60)
end # convert to m/s
clt = Raster("./case_study_1_present_day_validation/present_day_climate_topo_inputs/monmean_clt_climatology_1981-2010.nc")
rsds = Raster("./case_study_1_present_day_validation/present_day_climate_topo_inputs/monmean_rsds_climatology_1981-2010.nc")
topo = Raster("./case_study_1_present_day_validation/present_day_climate_topo_inputs/present_day_topography.nc")

# Additional arguments needed:
res = 0.5 # Target resolution
CO2 = 360.0 # Current atmospheric CO2, transferred to vegetation model
FDsampling = true # Assessment of functional diversity  
RIsampling = true # Assessment of species richness potential ("diversity index")
RI_landscape_window = 300.0 # Width of the landscape window used for the diversity assessment (in km)
outputdir = "./case_study_1_present_day_validation/TREED_present_output"

# Run TREED for the present 
TREED_output = TREEDsteadystep(tair=tair, precip=precip, clt=clt, rsds=rsds, topo=topo, CO2=CO2, res=res, FDsampling=FDsampling, RIsampling=RIsampling, RI_landscape_window=RI_landscape_window, outputdir=outputdir)

# Or read output from previously completed run 
TREED_output = RasterStack("./case_study_1_present_day_validation/TREED_present_output/TREED_output.nc")



###########################################################
########## Get validation data

topography = TREED_output.topography

NPP_ref = Raster("./case_study_1_present_day_validation/present_day_validation_data/NPP_ref.nc")
NPP_ref = replace_missing(NPP_ref, NaN)
NPP_ref = NPP_ref[:,:,1]
NPP_ref = resample(NPP_ref; to=TREED_output.H, method="bilinear")
NPP_ref[isnan.(TREED_output.H)] .= NaN
NPP_ref[isnan.(NPP_ref) .&& .!isnan.(TREED_output.H)] .= 0

GPP_ref = Raster("./case_study_1_present_day_validation/present_day_validation_data/GPP_ref.nc")
GPP_ref = replace_missing(GPP_ref, NaN)
GPP_ref = GPP_ref[:,:,1]
GPP_ref = resample(GPP_ref; to=TREED_output.H, method="bilinear")
GPP_ref[isnan.(TREED_output.H)] .= NaN
GPP_ref[isnan.(GPP_ref) .&& .!isnan.(TREED_output.H)] .= 0

AGB_ref = Raster("./case_study_1_present_day_validation/present_day_validation_data/AGB_ref.nc")
AGB_ref = replace_missing(AGB_ref, NaN)
AGB_ref = resample(AGB_ref; to=TREED_output.H, method="bilinear")
AGB_ref[isnan.(TREED_output.H)] .= NaN
AGB_ref[isnan.(AGB_ref) .&& .!isnan.(TREED_output.H)] .= 0
AGB_ref = AGB_ref .* 0.47 # multiply with assumed model carobn density in dry biomass -> Mg C/ha

BGB_ref = Raster("./case_study_1_present_day_validation/present_day_validation_data/BGB_ref.nc")
BGB_ref = replace_missing(BGB_ref, NaN)
BGB_ref = resample(BGB_ref; to=TREED_output.H, method="bilinear")
BGB_ref[isnan.(TREED_output.H)] .= NaN
BGB_ref[isnan.(BGB_ref) .&& .!isnan.(TREED_output.H)] .= 0
BGB_ref = BGB_ref .* 0.47 # multiply with assumed model carobn density in dry biomass -> Mg C/ha

H_ref = Raster("./case_study_1_present_day_validation/present_day_validation_data/H_ref.nc")
H_ref = replace_missing(H_ref, NaN)
H_ref = resample(H_ref; to=TREED_output.H, method="bilinear")
H_ref[isnan.(TREED_output.H)] .= NaN
H_ref[isnan.(H_ref) .&& .!isnan.(TREED_output.H)] .= 0

AET_ref = Raster("./case_study_1_present_day_validation/present_day_validation_data/AET_ref.nc")
AET_ref = replace_missing(AET_ref[:,:,1], NaN)
AET_ref = resample(AET_ref; to=TREED_output.H, method="bilinear")
AET_ref[isnan.(TREED_output.H)] .= NaN
AET_ref[isnan.(AET_ref) .&& .!isnan.(TREED_output.H)] .= 0


SR_ref = Raster("./case_study_1_present_day_validation/present_day_validation_data/SR_ref.nc")
SR_ref= replace_missing(SR_ref, NaN)
SR_ref = resample(SR_ref; to=TREED_output.H, method="bilinear")
SR_ref[isnan.(TREED_output.H)] .= NaN
SR_ref[isnan.(SR_ref) .&& .!isnan.(TREED_output.H)] .= 0

# Power model diversity index data 
valid_index = .!isnan.(vec(TREED_output.diversity_index)) .&& .!isnan.(vec(SR_ref)) .&& (vec(TREED_output.diversity_index) .> 0 .&& vec(SR_ref .> 0))
# Fit a linear model in log space
data = DataFrame(y = log10.(vec(SR_ref)[valid_index]), x = log10.(vec(TREED_output.diversity_index)[valid_index]))
linear_fit = lm(@formula(y ~ x), data)
cor(data.y, data.x)

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

########## Comparison fluxes 

GPP_cpt = makecpt(cmap=:viridis, range=(0, 3000), continuous=true, overrule_bg=true, par=(COLOR_NAN=230, COLOR_BACKGROUND="68/1/84", COLOR_FOREGROUND="251/231/35"))
grdimage(convert_raster_to_GMT_grid(TREED_output.GPP), projection=:Mollweide, theme="A2xy",
    cmap=GPP_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"))
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
grdimage!(convert_raster_to_GMT_grid(GPP_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=GPP_cpt, figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15))), frame=(annot=:auto, ticks=:auto, xlabel="GPP (g C m@+-2@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
diff_cpt = makecpt(cmap=:bam, range=(-2000, 2000),overrule_bg=true, continuous=true, par=(COLOR_NAN=230,COLOR_BACKGROUND="101/2/75", COLOR_FOREGROUND="15/77/1"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.GPP .- GPP_ref), projection=:Mollweide, theme="A2xy",
    cmap=diff_cpt, xaxis=(annot=0, ), yaxis=(annot=0,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15), triangles=true)), frame=(annot=:auto, ticks=:auto, xlabel="GPP difference (g C m@+-2@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
text!("(a)",frame=:none,region=(0,10,0,10),x=-21, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("Model",frame=:none,region=(0,10,0,10),x=-19, y=4.5, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(b)",frame=:none,region=(0,10,0,10),x=-10, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("Observation",frame=:none,region=(0,10,0,10),x=-9, y=4.5, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(c)",frame=:none,region=(0,10,0,10),x=0.5, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 


NPP_cpt = makecpt(cmap=:plasma, range=(0, 1300), continuous=true, overrule_bg=true, par=(COLOR_NAN=230, COLOR_BACKGROUND="15/7/136", COLOR_FOREGROUND="240/248/35"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.NPP), projection=:Mollweide, theme="A2xy",
    cmap=NPP_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), yshift=-4.5, xshift=-11.4)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
grdimage!(convert_raster_to_GMT_grid(NPP_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=NPP_cpt, figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15))), frame=(annot=:auto, ticks=:auto, xlabel="NPP (g C m@+-2@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
diff_cpt = makecpt(cmap=:bam, range=(-1000, 1000),overrule_bg=true, continuous=true, par=(COLOR_NAN=230,COLOR_BACKGROUND="101/2/75", COLOR_FOREGROUND="15/77/1"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.NPP .- NPP_ref), projection=:Mollweide, theme="A2xy",
    cmap=diff_cpt, xaxis=(annot=0, ), yaxis=(annot=0,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15), triangles=true)), frame=(annot=:auto, ticks=:auto, xlabel="NPP difference (g C m@+-2@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
text!("(d)",frame=:none,region=(0,10,0,10),x=-21, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(e)",frame=:none,region=(0,10,0,10),x=-10, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(f)",frame=:none,region=(0,10,0,10),x=0.5, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 



AET_cpt = makecpt(cmap=:devon, range=(0, 1500), continuous=true, overrule_bg=true, par=(COLOR_NAN=230, COLOR_BACKGROUND="44/26/76", COLOR_FOREGROUND="254/254/255"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.AET), projection=:Mollweide, theme="A2xy",
    cmap=AET_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), yshift=-4.5, xshift=-11.4)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
grdimage!(convert_raster_to_GMT_grid(AET_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=AET_cpt, figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15))), frame=(annot=:auto, ticks=:auto, xlabel="AET (mm year@+-1@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
diff_cpt = makecpt(cmap=:bam, range=(-1000, 1000),overrule_bg=true, continuous=true, par=(COLOR_NAN=230,COLOR_BACKGROUND="101/2/75", COLOR_FOREGROUND="15/77/1"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.AET .- AET_ref), projection=:Mollweide, theme="A2xy",
    cmap=diff_cpt, xaxis=(annot=0, ), yaxis=(annot=0,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15), triangles=true)), frame=(annot=:auto, ticks=:auto, xlabel="AET difference (mm year@+-1@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
text!("(g)",frame=:none,region=(0,10,0,10),x=-21, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(h)",frame=:none,region=(0,10,0,10),x=-10, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(i)",frame=:none,region=(0,10,0,10),x=0.5, y=3.75, noclip=true ,font=(10,"Helvetica",:black), 
show=true, dpi=330, name="./case_study_1_present_day_validation/plots/flux_comparison_spatial.png") 



########## Comparison structures
H_cpt = makecpt(cmap=:bamako, range=(-1, 55), inverse=true, continuous=true, overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="255/254/254", COLOR_FOREGROUND="255/254/254"))
grdimage(convert_raster_to_GMT_grid(TREED_output.H), projection=:Mollweide, theme="A2xy",
    cmap=H_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"))
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
grdimage!(convert_raster_to_GMT_grid(H_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=H_cpt, figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15))), frame=(annot=:auto, ticks=:auto, xlabel="H (m)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
diff_cpt = makecpt(cmap=:bam, range=(-20, 20),overrule_bg=true, continuous=true, par=(COLOR_NAN=230,COLOR_BACKGROUND="101/2/75", COLOR_FOREGROUND="15/77/1"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.H .- H_ref), projection=:Mollweide, theme="A2xy",
    cmap=diff_cpt, xaxis=(annot=0, ), yaxis=(annot=0,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15), triangles=true)), frame=(annot=:auto, ticks=:auto, xlabel="H difference (m)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
text!("(a)",frame=:none,region=(0,10,0,10),x=-21, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("Model",frame=:none,region=(0,10,0,10),x=-19, y=4.5, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(b)",frame=:none,region=(0,10,0,10),x=-10, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("Observation",frame=:none,region=(0,10,0,10),x=-9, y=4.5, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(c)",frame=:none,region=(0,10,0,10),x=0.5, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 

AGB_model = ((TREED_output.C_leaf .+ TREED_output.C_heartwood .+ TREED_output.C_sapwood) ./ TREED_output.CA) .* (10000/1e+6) # in Mg C ha-1
AGB_model[isnan.(AGB_model) .&& .!isnan.(TREED_output.H)] .= 0

AGB_cpt = makecpt(cmap=:roma, reverse=true, continuous=true, range=(-1, 150), overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="4/50/152", COLOR_FOREGROUND="126/24/0"))
grdimage!(convert_raster_to_GMT_grid(AGB_model), projection=:Mollweide, theme="A2xy",
    cmap=AGB_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), yshift=-4.5, xshift=-11.4)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
grdimage!(convert_raster_to_GMT_grid(AGB_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=AGB_cpt, figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15))), frame=(annot=:auto, ticks=:auto, xlabel="AGB (Mg C ha@+-1@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
diff_cpt = makecpt(cmap=:bam, range=(-100, 100),overrule_bg=true, continuous=true, par=(COLOR_NAN=230,COLOR_BACKGROUND="101/2/75", COLOR_FOREGROUND="15/77/1"))
grdimage!(convert_raster_to_GMT_grid(AGB_model .- AGB_ref), projection=:Mollweide, theme="A2xy",
    cmap=diff_cpt, xaxis=(annot=0, ), yaxis=(annot=0,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15), triangles=true)), frame=(annot=:auto, ticks=:auto, xlabel="AGB difference (Mg C ha@+-1@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
text!("(d)",frame=:none,region=(0,10,0,10),x=-21, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(e)",frame=:none,region=(0,10,0,10),x=-10, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(f)",frame=:none,region=(0,10,0,10),x=0.5, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 


BGB_model = ((TREED_output.C_coarseroot .+ TREED_output.C_fineroot) ./ TREED_output.CA) .* (10000/1e+6) # in Mg C ha-1
BGB_model[isinf.(BGB_model) .&& .!isnan.(TREED_output.H)] .= 0

BGB_cpt = makecpt(cmap=:roma, reverse=true, range=(-1, 40), overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="4/50/152", COLOR_FOREGROUND="126/24/0"))
grdimage!(convert_raster_to_GMT_grid(BGB_model), projection=:Mollweide, theme="A2xy",
    cmap=BGB_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), yshift=-4.5, xshift=-11.4)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
grdimage!(convert_raster_to_GMT_grid(BGB_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=BGB_cpt, figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15))), frame=(annot=:auto, ticks=:auto, xlabel="BGB (Mg C ha@+-1@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
diff_cpt = makecpt(cmap=:bam, range=(-30, 30),overrule_bg=true, continuous=true, par=(COLOR_NAN=230,COLOR_BACKGROUND="101/2/75", COLOR_FOREGROUND="15/77/1"))
grdimage!(convert_raster_to_GMT_grid(BGB_model .- BGB_ref), projection=:Mollweide, theme="A2xy",
    cmap=diff_cpt, xaxis=(annot=0, ), yaxis=(annot=0,), figsize=5.6, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"), xshift=5.7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    colorbar!(pos=(position=(outside=true, anchor=:BC, size=(2, 0.15), triangles=true)), frame=(annot=:auto, ticks=:auto, xlabel="BGB difference (Mg C ha@+-1@+)"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=20))
text!("(g)",frame=:none,region=(0,10,0,10),x=-21, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(h)",frame=:none,region=(0,10,0,10),x=-10, y=3.75, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(i)",frame=:none,region=(0,10,0,10),x=0.5, y=3.75, noclip=true ,font=(10,"Helvetica",:black), 
show=true, dpi=330, name="./case_study_1_present_day_validation/plots/structure_comparison_spatial.png") 



########## 1:1 plots structures

# Height density plot $
RMSE = sqrt(mean(filter(!isnan, vec(TREED_output.H .- H_ref).^2)))
D = binstats([vec(H_ref)[valid_index] convert.(Float64, vec(TREED_output.H)[valid_index])], inc=1.5, tiling=:hex, stats=:number)
density_cpt=makecpt(cmap=:nuuk, reverse=true, continuous=true, range=(0, 200), par=(COLOR_BACKGROUND="5/89/140", COLOR_FOREGROUND="253/253/177"))
GMT.plot(D, region=(0, 40, 0, 40), hexbin=true,colorbar=false, cmap=density_cpt, figsize=(6, 6), theme=("A2xy"),
    xlabel="Observed height (m)", ylabel="Modelled height (m)",
    par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6, MAP_FRAME_PEN="0.3p"))
GMT.plot!([0, 46], [0, 46], linecolor=:darkred, lw=0.8)
colorbar!(pos=(position=(outside=true, anchor=:MR, offset=(0.25,-2), size=(2, 0.15), triangles=:f)), frame=(annot=:auto, ticks=:auto, xlabel="Count"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=16))
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," m"), x = 32, y = 2, font=7)


AGB_model = ((TREED_output.C_leaf .+ TREED_output.C_heartwood .+ TREED_output.C_sapwood) ./ TREED_output.CA) .* (10000/1e+6) # in Mg C ha-1
AGB_model[isnan.(AGB_model) .&& .!isnan.(TREED_output.H)] .= 0
RMSE = sqrt(mean(filter(!isnan, vec(AGB_model .- AGB_ref).^2)))
D = binstats([vec(AGB_ref)[valid_index] convert.(Float64, vec(AGB_model)[valid_index])], inc=2.5, tiling=:hex, stats=:number)
density_cpt=makecpt(cmap=:nuuk, reverse=true, continuous=true, range=(0, 70), par=(COLOR_BACKGROUND="5/89/140", COLOR_FOREGROUND="253/253/177"))
GMT.plot!(D, region=(0, 150, 0, 150), hexbin=true,colorbar=false, cmap=density_cpt, figsize=(6, 6), theme=("A2xy"),
    xlabel="Observed AGB (Mg C ha@+-1@+)", ylabel="Modelled AGB (Mg C ha@+-1@+)",
    par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6, MAP_FRAME_PEN="0.3p") , yshift=-7)
GMT.plot!([0, 160], [0, 160], linecolor=:darkred, lw=0.8)
colorbar!(pos=(position=(outside=true, anchor=:MR, offset=(0.25,-2), size=(2, 0.15), triangles=:f)), frame=(annot=:auto, ticks=:auto, xlabel="Count"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=16))
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," Mg C ha@+-1@+"), x = 110, y = 8, font=7)


BGB_model = ((TREED_output.C_coarseroot .+ TREED_output.C_fineroot) ./ TREED_output.CA) .* (10000/1e+6) # in Mg C ha-1
BGB_model[isinf.(BGB_model) .&& .!isnan.(TREED_output.H)] .= 0
RMSE = sqrt(mean(filter(!isnan, vec(BGB_model .- BGB_ref).^2)))
D = binstats([vec(BGB_ref)[valid_index] convert.(Float64, vec(BGB_model)[valid_index])], inc=1, tiling=:hex, stats=:number)
density_cpt=makecpt(cmap=:nuuk, reverse=true, continuous=true, range=(0, 90), par=(COLOR_BACKGROUND="5/89/140", COLOR_FOREGROUND="253/253/177"))
GMT.plot!(D, region=(0, 42, 0, 42), hexbin=true,colorbar=false, cmap=density_cpt, figsize=(6, 6), theme=("A2xy"),
    xlabel="Observed BGB (Mg C ha@+-1@+)", ylabel="Modelled BGB (Mg C ha@+-1@+)",
    par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6, MAP_FRAME_PEN="0.3p") , yshift=-7)
GMT.plot!([0, 160], [0, 160], linecolor=:darkred, lw=0.8)
colorbar!(pos=(position=(outside=true, anchor=:MR, offset=(0.25,-2), size=(2, 0.15), triangles=:f)), frame=(annot=:auto, ticks=:auto, xlabel="Count"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=16))
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," Mg C ha@+-1@+"), x = 32, y = 2, font=7)
GMT.text!("(f)", x = -8, y = 41, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(e)", x = -8, y = 91, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(d)", x = -8, y = 141, noclip=true, font=(10,"Helvetica",:black))



########## 1:1 plots fluxes
RMSE = sqrt(mean(filter(!isnan, vec(TREED_output.GPP .- GPP_ref).^2)))
D = binstats([vec(GPP_ref)[valid_index] convert.(Float64, vec(TREED_output.GPP)[valid_index])], inc=50, tiling=:hex, stats=:number)
density_cpt=makecpt(cmap=:nuuk, reverse=true, continuous=true, range=(0, 50), par=(COLOR_BACKGROUND="5/89/140", COLOR_FOREGROUND="253/253/177"))
GMT.plot!(D, region=(0, 3000, 0, 3000), hexbin=true,colorbar=false, cmap=density_cpt, figsize=(6, 6), theme=("A2xy"),
    xlabel="Observed GPP (g C m@+-2@+)", ylabel="Modelled GPP (g C m@+-2@+)",
    par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6, MAP_FRAME_PEN="0.3p"), xshift=-8, yshift=14)
GMT.plot!([0, 3100], [0, 3100], linecolor=:darkred, lw=0.8)
colorbar!(pos=(position=(outside=true, anchor=:MR, offset=(0.25,-2), size=(2, 0.15), triangles=:f)), frame=(annot=:auto, ticks=:auto, xlabel="Count"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=16))
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," g C m@+-2@+"), x = 2200, y = 180, font=7)


RMSE = sqrt(mean(filter(!isnan, vec(TREED_output.NPP .- NPP_ref).^2)))
D = binstats([vec(NPP_ref)[valid_index] convert.(Float64, vec(TREED_output.NPP)[valid_index])], inc=25, tiling=:hex, stats=:number)
density_cpt=makecpt(cmap=:nuuk, reverse=true, continuous=true, range=(0, 50), par=(COLOR_BACKGROUND="5/89/140", COLOR_FOREGROUND="253/253/177"))
GMT.plot!(D, region=(0, 1800, 0, 1800), hexbin=true,colorbar=false, cmap=density_cpt, figsize=(6, 6), theme=("A2xy"),
    xlabel="Observed NPP (g C m@+-2@+)", ylabel="Modelled NPP (g C m@+-2@+)",
    par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6, MAP_FRAME_PEN="0.3p") , yshift=-7)
GMT.plot!([0, 2000], [0, 2000], linecolor=:darkred, lw=0.8)
colorbar!(pos=(position=(outside=true, anchor=:MR, offset=(0.25,-2), size=(2, 0.15), triangles=:f)), frame=(annot=:auto, ticks=:auto, xlabel="Count"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=16))
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," g C m@+-2@+"), x = 1380, y = 100, font=7)

RMSE = sqrt(mean(filter(!isnan, vec(TREED_output.AET .- AET_ref).^2)))
D = binstats([vec(AET_ref)[valid_index] convert.(Float64, vec(TREED_output.AET)[valid_index])], inc=20, tiling=:hex, stats=:number)
density_cpt=makecpt(cmap=:nuuk, reverse=true, continuous=true, range=(0, 40), par=(COLOR_BACKGROUND="5/89/140", COLOR_FOREGROUND="253/253/177"))
GMT.plot!(D, region=(0, 1500, 0, 1500), hexbin=true,colorbar=false, cmap=density_cpt, figsize=(6, 6), theme=("A2xy"),
    xlabel="Observed AET (mm year@+-1@+)", ylabel="Modelled AET (mm year@+-1@+)",
    par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6, MAP_FRAME_PEN="0.3p") , yshift=-7)
GMT.plot!([0, 1500], [0, 1500], linecolor=:darkred, lw=0.8)
colorbar!(pos=(position=(outside=true, anchor=:MR, offset=(0.25,-2), size=(2, 0.15), triangles=:f)), frame=(annot=:auto, ticks=:auto, xlabel="Count"),par=(FONT_ANNOT_PRIMARY=14,FONT_LABEL=16))
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," mm year@+-1@+"), x = 1090, y = 70, font=7)
GMT.text!("(c)", x = -230, y = 1450, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(b)", x = -230, y = 3210, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(a)", x = -240, y = 5000, noclip=true, font=(10,"Helvetica",:black),
show=true, dpi=330, name="./case_study_1_present_day_validation/plots/flux_structure_relationships_onetoone.png")


########## Flux~structure relationships 
# See R script for final version of density plots (in folder plots)
valid_index = .!isnan.(vec(H_ref)) .&& .!isnan.(vec(NPP_ref)) .&& (vec(H_ref) .> 0 .&& vec(NPP_ref .> 0))
GMT.basemap(region=(1, 40, 0, 2000), figsize=(6, 5), xlabel="Height (m)", ylabel="NPP (g C m@+-2@+)",
theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6, MAP_FRAME_PEN="0.3p"))
GMT.scatter!(vec(H_ref)[valid_index], vec(NPP_ref)[valid_index], fill=:black, alpha=95, legend=(label="Data", pos=:TL, box=:none,))
valid_index = .!isnan.(vec(TREED_output.H)) .&& .!isnan.(vec(TREED_output.NPP)) .&& (vec(TREED_output.H) .> 0 .&& vec(TREED_output.NPP .> 0))
GMT.scatter!(vec(TREED_output.H)[valid_index], vec(TREED_output.NPP)[valid_index], fill=:lightblue, alpha=98, legend="Model", show=true)

valid_index = .!isnan.(vec(H_ref)) .&& .!isnan.(vec(AGB_ref)) .&& (vec(H_ref) .> 0 .&& vec(AGB_ref .> 0))
GMT.basemap!(region=(1, 50, 0, 150), figsize=(6, 5), xlabel="Height (m)", ylabel="AGB (Mg C ha@+-1@+)",
    theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6, MAP_FRAME_PEN="0.3p"), xshift = 7.25)
GMT.scatter!(vec(H_ref)[valid_index], vec(AGB_ref)[valid_index], fill=:black, alpha=95, legend=(label="Data", pos=:TL, box=:none,))
valid_index = .!isnan.(vec(TREED_output.H)) .&& .!isnan.(vec(AGB_model)) .&& (vec(TREED_output.H) .> 0 .&& vec(AGB_model .> 0))
GMT.scatter!(convert.(Float64, vec(TREED_output.H)[valid_index]), vec(AGB_model)[valid_index], fill=:lightblue, alpha=98, legend="Model") 
GMT.text!("(b)", x = -5.5, y = 152, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(a)", x = -68, y = 152, noclip=true, font=(10,"Helvetica",:black),show=true)


########## Richness estimation 
FD_cpt = makecpt(cmap=:batlowK, continuous=true, range=(0, 0.25), overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="4/6/11", COLOR_FOREGROUND="250/204/249"))
grdimage(convert_raster_to_GMT_grid(TREED_output.functional_diversity), projection=:Mollweide, theme="A2xy",
    cmap=FD_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=6.5, par=(FONT_ANNOT=6, MAP_FRAME_PEN="0.2p"))
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
colorbar!(pos=(paper=true, anchor=(3.5,-0.1), size=(4,0.2), justify=:TC, horizontal=true, triangles=:f), frame=(annot=:auto, ), xlabel="Functional diversity index (-)", 
    par=(FONT_LABEL=16, FONT_ANNOT_PRIMARY=12))

EH_cpt = makecpt(cmap=:batlowK, continuous=true, range=(0, 0.5), overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="4/6/11", COLOR_FOREGROUND="250/204/249"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.gamma_EH), projection=:Mollweide, theme="A2xy",
    cmap=EH_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=6.5, par=(FONT_ANNOT=7, MAP_FRAME_PEN="0.2p"), xshift=7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
colorbar!(pos=(paper=true, anchor=(3.5,-0.1), size=(4,0.2), justify=:TC, horizontal=true, triangles=:f), frame=(annot=:auto, ), xlabel="Landscape heterogeneity (-)", 
    par=(FONT_LABEL=16, FONT_ANNOT_PRIMARY=12))

GI_cpt = makecpt(cmap=:batlowK, continuous=true, range=(0, 1.0), overrule_bg=true, par=(COLOR_NAN=235, COLOR_BACKGROUND="4/6/11", COLOR_FOREGROUND="250/204/249"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.gamma_GI), projection=:Mollweide, theme="A2xy",
    cmap=GI_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=6.5, par=(FONT_ANNOT=7, MAP_FRAME_PEN="0.2p"), xshift=7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
colorbar!(pos=(paper=true, anchor=(3.5,-0.1), size=(4,0.2), justify=:TC, horizontal=true, triangles=:f), frame=(annot=:auto, ), xlabel="Landscape fragmentation (-)", 
    par=(FONT_LABEL=16, FONT_ANNOT_PRIMARY=12))

SR_cpt = makecpt(cmap=:batlowK, continuous = true, 
    range=(0, 2000), overrule_bg=true, par=(COLOR_NAN=235, COLOR_BACKGROUND="4/6/11", COLOR_FOREGROUND="163.0/19.0/1.0"))

grdimage!(convert_raster_to_GMT_grid((10^3.8355) .* (TREED_output.diversity_index .^ 0.302779)), projection=:Mollweide, theme="A2xy",
    cmap=SR_cpt, xaxis=(annot=0,), yaxis=(annot=60,), figsize=10, par=(FONT_ANNOT=7, MAP_FRAME_PEN="0.2p"), xshift=-14, yshift=-7)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
    text!("Z = 6847 × DP@+0.30@+",frame=:none,region=(0,10,0,10), proj=:linear, x=5, y=5.3, noclip=true ,font=(8,"Helvetica",:black)) 
    text!("Model",frame=:none,region=(0,10,0,10), proj=:linear, x=4.5, y=5.7, noclip=true ,font=(9,"Helvetica",:black)) 
    text!("Data",frame=:none,region=(0,10,0,10), proj=:linear, x=12, y=5.3, noclip=true ,font=(9,"Helvetica",:black)) 
    colorbar!(pos=(paper=true, anchor=(10,-0.2), size=(4,0.2), justify=:TC, horizontal=true, triangles=:f), frame=(annot=:auto, ), xlabel="Species richness density", 
    par=(FONT_LABEL=16, FONT_ANNOT_PRIMARY=12))

grdimage!(convert_raster_to_GMT_grid(SR_ref), projection=:Mollweide, theme="A2xy",
    cmap=SR_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=10, par=(FONT_ANNOT=7, MAP_FRAME_PEN="0.2p"), xshift=10.5)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")

image!("./case_study_1_present_day_validation/plots/richness_density_plot.png", frame=:none, xshift=-16, yshift=-8.25)

GMT.text!("(d)",frame=:none,region=(0,10,0,10), x=5, y=11,  noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(e)",frame=:none,region=(0,10,0,10), x=17, y=11,  noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(f)",frame=:none,region=(0,10,0,10), x=5, y=6,  noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(a)",frame=:none,region=(0,10,0,10), x=5, y=15, noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(b)",frame=:none,region=(0,10,0,10), x=13, y=15, noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(c)",frame=:none,region=(0,10,0,10), x=20, y=15, noclip=true, font=(10, "Helvetica", :black), 
show=true, dpi=330, name="./case_study_1_present_day_validation/plots/richness_estimation_validation.png")


##########  Plot additional structures 
C = makecpt(cmap=((33,113,181), (254,228,171)), T=[-0.1,0.5,1.1],overrule_bg=true, par=(COLOR_NAN=235, COLOR_BACKGROUND="56/102/149", COLOR_FOREGROUND="15/42/3"))
C.label = ["Deciduous","Evergreen"]
grdimage(convert_raster_to_GMT_grid(TREED_output.seasonality), projection=:Mollweide, theme="A2xy",
    cmap=C, xaxis=(annot=0,), yaxis=(annot=60,), figsize=7.25, par=(FONT_ANNOT=7, MAP_FRAME_PEN="0.2p"))
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
colorbar!(pos=(paper=true, anchor=(4,-0.2), size=(4,0.2), justify=:TC, horizontal=true),
    B=:none, equal_size=(gap=0.1,),par=(FONT_ANNOT=12,))

a_ll_cpt = makecpt(cmap=:cork, hinge=1, range=(0, 5), overrule_bg=true, par=(COLOR_NAN=235, COLOR_BACKGROUND="56/102/149", COLOR_FOREGROUND="15/42/3"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.a_ll), projection=:Mollweide, theme="A2xy",
    cmap=a_ll_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=7.25, par=(FONT_ANNOT=7, MAP_FRAME_PEN="0.2p"), xshift=9.5)
    grdcontour!(convert_raster_to_GMT_grid(TREED_output.topography), projection=:Mollweide, levels=[0], pen="0.08p,black")
colorbar!(pos=(paper=true, anchor=(4,-0.2), size=(4,0.2), justify=:TC, horizontal=true), xlabel="Leaf longevity (years)",
par=(FONT_ANNOT=12,))
text!("(a)",frame=:none,region=(0,10,0,10), proj=:linear, x=-6, y=4, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(b)",frame=:none,region=(0,10,0,10), proj=:linear, x=-0.25, y=4, noclip=true ,font=(10,"Helvetica",:black))

image!("./case_study_1_present_day_validation/plots/lma_comparison_model_data.png", yshift=-8, xshift=-19)
text!("(c)",frame=:none,region=(0,10,0,10), proj=:linear, x=6.5, y=7, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(d)",frame=:none,region=(0,10,0,10), proj=:linear, x=12.5, y=7, noclip=true ,font=(10,"Helvetica",:black), 
dpi=330, name="./case_study_1_present_day_validation/plots/phenology_all.png", show=true)



# Additional plot for methods: illstruate SLA~a_ll and k~a_ll relationship 
a_ll = 0.001:0.01:5
SLA_model = (2e-4) .* (1/0.4763) .* 10 .^(2.25 .- 0.5 .* log10.(a_ll .* 12))
SLA_LPJ_needle = (2e-4) .* (1/0.4763) .* 10 .^(2.08 .- 0.4 .* log10.(a_ll .* 12))
SLA_LPJ_broad = (2e-4) .* (1/0.4763) .* 10 .^(2.22 .- 0.4 .* log10.(a_ll .* 12))

GMT.basemap(region=(0.001, 4.5, 0.009, 0.042), figsize=(6, 5), xlabel="a@-ll@- (years)", ylabel="SLA (m@+2@+ / g C)",
theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6, MAP_FRAME_PEN="0.2p"))
GMT.plot!(a_ll, SLA_LPJ_needle, linecolor=:darkgreen, lw=1.5, legend=(label="LPJ needle", pos=:TR, box=:none,))
GMT.plot!(a_ll, SLA_LPJ_broad, linecolor=:lightred, lw=1.5, legend=(label="LPJ broad", pos=:TR, box=:none,))
GMT.plot!(a_ll, SLA_model, linecolor=:black, lw=1.5, legend=(label="TREED", pos=:TR, box=:none,), 
dpi=700, name="./case_study_1_present_day_validation/plots/all_SLA_relationship.png", show=true)
