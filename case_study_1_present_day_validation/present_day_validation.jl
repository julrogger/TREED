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
# Derived from CHELSA: https://chelsa-climate.org/; Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, P., Kessler, M. (2017): Climatologies at high resolution for the Earth land surface areas. Scientific Data. 4 170122. https://doi.org/10.1038/sdata.2017.122
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
NPP_cpt = makecpt(cmap=:plasma, range=(0, 1300), overrule_bg=true, par=(COLOR_NAN=230, COLOR_BACKGROUND="15/7/136", COLOR_FOREGROUND="240/248/35"))
grdimage(convert_raster_to_GMT_grid(TREED_output.NPP), projection=:Mollweide, theme="A2xy",
    cmap=NPP_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=6.5, par=(FONT_ANNOT=7,))
grdimage!(convert_raster_to_GMT_grid(NPP_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=NPP_cpt, figsize=6.5, xshift=6.6)
colorbar!(pos=(achor=:RM,), frame=(annot=:auto, ticks=:auto, xlabel="NPP (g C m@+-2@+)"),par=(FONT_ANNOT_PRIMARY=12,))
text!("(a)",frame=:none,region=(0,10,0,10), xshift=-6.6, proj=:linear, x=-0.5, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(b)",frame=:none,region=(0,10,0,10), proj=:linear, x=4.25, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 

GPP_cpt = makecpt(cmap=:viridis, range=(0, 3000), overrule_bg=true, par=(COLOR_NAN=230, COLOR_BACKGROUND="68/1/84", COLOR_FOREGROUND="251/231/35"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.GPP), projection=:Mollweide, theme="A2xy",
    cmap=GPP_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=6.5, par=(FONT_ANNOT=7,), yshift=-4)
grdimage!(convert_raster_to_GMT_grid(GPP_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=GPP_cpt, figsize=6.5, xshift=6.6)
colorbar!(pos=(achor=:RM,), frame=(annot=:auto, ticks=:auto, xlabel="GPP (g C m@+-2@+)"),par=(FONT_ANNOT_PRIMARY=12,))
text!("(c)",frame=:none,region=(0,10,0,10), xshift=-6.6, proj=:linear, x=-0.5, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(d)",frame=:none,region=(0,10,0,10), proj=:linear, x=4.25, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 

AET_cpt = makecpt(cmap=:devon, range=(0, 1500), overrule_bg=true, par=(COLOR_NAN=230, COLOR_BACKGROUND="44/26/76", COLOR_FOREGROUND="254/254/255"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.AET), projection=:Mollweide, theme="A2xy",
    cmap=AET_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=6.5, par=(FONT_ANNOT=7,), yshift=-4)
grdimage!(convert_raster_to_GMT_grid(AET_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=AET_cpt, figsize=6.5, xshift=6.6)
colorbar!(pos=(achor=:RM,), frame=(annot=:auto, ticks=:auto, xlabel="AET (mm year@+-1@+)"),par=(FONT_ANNOT_PRIMARY=12,))
text!("(e)",frame=:none,region=(0,10,0,10), xshift=-6.6, proj=:linear, x=-0.5, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(f)",frame=:none,region=(0,10,0,10), proj=:linear, x=4.25, y=3.25, noclip=true ,font=(10,"Helvetica",:black), 
dpi=700, name="./case_study_1_present_day_validation/plots/fluxes_comparison.png") 


########## Comparison structures
H_cpt = makecpt(cmap=:bamako, range=(-1, 55), inverse=true, overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="255/254/254", COLOR_FOREGROUND="255/254/254"))
grdimage(convert_raster_to_GMT_grid(TREED_output.H), projection=:Mollweide, theme="A2xy",
    cmap=H_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=6.5, par=(FONT_ANNOT=7,))
grdimage!(convert_raster_to_GMT_grid(H_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=H_cpt, figsize=6.5, xshift=6.6)
colorbar!(pos=(achor=:RM,), frame=(annot=:auto, ticks=:auto, xlabel="H (m)"),par=(FONT_ANNOT_PRIMARY=12,))
text!("(a)",frame=:none,region=(0,10,0,10), xshift=-6.6, proj=:linear, x=-0.5, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(b)",frame=:none,region=(0,10,0,10), proj=:linear, x=4.25, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 


AGB_model = ((TREED_output.C_leaf .+ TREED_output.C_heartwood .+ TREED_output.C_sapwood) ./ TREED_output.CA) .* (10000/1e+6) # in Mg C ha-1
AGB_model[isnan.(AGB_model) .&& .!isnan.(TREED_output.H)] .= 0

AGB_cpt = makecpt(cmap=:roma, reverse=true, range=(-1, 150), overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="4/50/152", COLOR_FOREGROUND="126/24/0"))
grdimage!(convert_raster_to_GMT_grid(AGB_model), projection=:Mollweide, theme="A2xy",
    cmap=AGB_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=6.5, par=(FONT_ANNOT=7,), yshift=-4)
grdimage!(convert_raster_to_GMT_grid(AGB_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=AGB_cpt, figsize=6.5, xshift=6.6)
colorbar!(pos=(achor=:RM,), frame=(annot=:auto, ticks=:auto, xlabel="AGB (Mg C ha@+-1@+)"),par=(FONT_ANNOT_PRIMARY=12,))
text!("(c)",frame=:none,region=(0,10,0,10), xshift=-6.6, proj=:linear, x=-0.5, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(d)",frame=:none,region=(0,10,0,10), proj=:linear, x=4.25, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 


BGB_model = ((TREED_output.C_coarseroot .+ TREED_output.C_fineroot) ./ TREED_output.CA) .* (10000/1e+6) # in Mg C ha-1
BGB_model[isinf.(BGB_model) .&& .!isnan.(TREED_output.H)] .= 0

BGB_cpt = makecpt(cmap=:roma, reverse=true, range=(-1, 40), overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="4/50/152", COLOR_FOREGROUND="126/24/0"))
grdimage!(convert_raster_to_GMT_grid(BGB_model), projection=:Mollweide, theme="A2xy",
    cmap=BGB_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=6.5, par=(FONT_ANNOT=7,), yshift=-4)
grdimage!(convert_raster_to_GMT_grid(BGB_ref), yaxis=(annot=0, ), xaxis=(annot=0,), projection=:Mollweide, theme="A2xy",
    cmap=BGB_cpt, figsize=6.5, xshift=6.6)
colorbar!(pos=(achor=:RM,), frame=(annot=:auto, ticks=:auto, xlabel="BGB (Mg C ha@+-1@+)"),par=(FONT_ANNOT_PRIMARY=12,))
text!("(e)",frame=:none,region=(0,10,0,10), xshift=-6.6, proj=:linear, x=-0.5, y=3.25, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(f)",frame=:none,region=(0,10,0,10), proj=:linear, x=4.25, y=3.25, noclip=true ,font=(10,"Helvetica",:black), 
dpi=700, name="./case_study_1_present_day_validation/plots/structures_comparison.png") 


########## 1:1 plots structures
valid_index = .!isnan.(vec(TREED_output.H)) .&& .!isnan.(vec(H_ref))
RMSE = sqrt(sum(filter(!isnan, vec(TREED_output.H .- H_ref).^2)) / length(vec(TREED_output.H .- H_ref)))
GMT.basemap(region=(0, 46, 0, 46), figsize=(6, 5), xlabel="Observed height (m)", ylabel="Modelled height (m)", theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,))
GMT.scatter!(vec(H_ref)[valid_index], convert.(Float64, vec(TREED_output.H)[valid_index]), fill=:black, alpha=98)
GMT.plot!([0, 46], [0, 46], linecolor=:darkred, lw=0.8)
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," m"), x = 38.5, y = 2, font=7)

AGB_model = ((TREED_output.C_leaf .+ TREED_output.C_heartwood .+ TREED_output.C_sapwood) ./ TREED_output.CA) .* (10000/1e+6) # in Mg C ha-1
AGB_model[isnan.(AGB_model) .&& .!isnan.(TREED_output.H)] .= 0
valid_index = .!isnan.(vec(AGB_model)) .&& .!isnan.(vec(AGB_ref))
RMSE = sqrt(sum(filter(!isnan, vec(AGB_model .- AGB_ref).^2)) / length(vec(AGB_model .- AGB_ref)))
GMT.basemap!(region=(-1, 160, -1, 160), figsize=(6, 5), xlabel="Observed AGB (Mg C ha@+-1@+)", ylabel="Modelled AGB (Mg C ha@+-1@+)",
theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,), yshift=-6)
GMT.scatter!(vec(AGB_ref)[valid_index], convert.(Float64, vec(AGB_model)[valid_index]), fill=:black, alpha=98)
GMT.plot!([0, 160], [0, 160], linecolor=:darkred, lw=0.8)
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," Mg C ha@+-1@+"), x = 120, y = 10, font=7)

BGB_model = ((TREED_output.C_coarseroot .+ TREED_output.C_fineroot) ./ TREED_output.CA) .* (10000/1e+6) # in Mg C ha-1
BGB_model[isinf.(BGB_model) .&& .!isnan.(TREED_output.H)] .= 0
valid_index = .!isnan.(vec(BGB_model)) .&& .!isnan.(vec(BGB_ref))
RMSE = sqrt(sum(filter(!isnan, vec(BGB_model .- BGB_ref).^2)) / length(vec(BGB_model .- BGB_ref)))
GMT.basemap!(region=(0, 42, 0, 42), figsize=(6, 5), xlabel="Observed BGB (Mg C ha@+-1@+)", ylabel="Modelled BGB (Mg C ha@+-1@+)",
theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,), yshift=-6)
GMT.scatter!(vec(BGB_ref)[valid_index], convert.(Float64, vec(BGB_model)[valid_index]), fill=:black, alpha=98)
GMT.plot!([0, 42], [0, 42], linecolor=:darkred, lw=0.8)
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," Mg C ha@+-1@+"), x = 32, y = 2, font=7)
GMT.text!("(f)", x = -8, y = 44, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(e)", x = -8, y = 94, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(d)", x = -8, y = 144, noclip=true, font=(10,"Helvetica",:black))

########## 1:1 plots fluxes
valid_index = .!isnan.(vec(TREED_output.GPP)) .&& .!isnan.(vec(GPP_ref)) .&& (vec(GPP_ref) .> 0 .&& vec(TREED_output.GPP .> 0))
RMSE = sqrt(sum(filter(!isnan, vec(TREED_output.GPP .- GPP_ref).^2)) / length(vec(TREED_output.GPP .- GPP_ref)))
GMT.basemap!(region=(0, 3000, 0, 3000), figsize=(6, 5), xlabel="Observed GPP (g C m@+-2@+)", ylabel="Modelled GPP (g C m@+-2@+)",
theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,), xshift=-7.5, yshift=12)
GMT.scatter!(vec(GPP_ref)[valid_index], convert.(Float64, vec(TREED_output.GPP)[valid_index]), fill=:black, alpha=98)
GMT.plot!([0, 3100], [0, 3100], linecolor=:darkred, lw=0.8)
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," g C m@+-2@+"), x = 2200, y = 200, font=7)

valid_index = .!isnan.(vec(TREED_output.NPP)) .&& .!isnan.(vec(NPP_ref)) .&& (vec(NPP_ref) .> 0 .&& vec(TREED_output.NPP .> 0))
RMSE = sqrt(sum(filter(!isnan, vec(TREED_output.NPP .- NPP_ref).^2)) / length(vec(TREED_output.NPP .- NPP_ref)))
GMT.basemap!(region=(0, 2000, 0, 2000), figsize=(6, 5), xlabel="Observed NPP (g C m@+-2@+)", ylabel="Modelled NPP (g C m@+-2@+)", theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,), yshift=-6)
GMT.scatter!(vec(NPP_ref)[valid_index], convert.(Float64, vec(TREED_output.NPP)[valid_index]), fill=:black, alpha=98)
GMT.plot!([0, 2000], [0, 2000], linecolor=:darkred, lw=0.8)
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," g C m@+-2@+"), x = 500, y = 1900, font=7)

valid_index = .!isnan.(vec(TREED_output.AET)) .&& .!isnan.(vec(AET_ref)) .&& (vec(AET_ref) .> 0 .&& vec(TREED_output.AET .> 0))
RMSE = sqrt(sum(filter(!isnan, vec(TREED_output.AET .- AET_ref).^2)) / length(vec(TREED_output.AET .- AET_ref)))
GMT.basemap!(region=(0, 1500, 0, 1500), figsize=(6, 5), xlabel="Observed AET (mm year@+-1@+)", ylabel="Modelled AET (mm year@+-1@+)", theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,), yshift=-6)
GMT.scatter!(vec(AET_ref)[valid_index], convert.(Float64, vec(TREED_output.AET)[valid_index]), fill=:black, alpha=98)
GMT.plot!([0, 1500], [0, 1500], linecolor=:darkred, lw=0.8)
GMT.text!(string("RMSE = ",round(RMSE,digits=2)," mm year@+-1@+"), x = 1120, y = 70, font=7)
GMT.text!("(c)", x = -230, y = 1580, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(b)", x = -230, y = 3350, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(a)", x = -230, y = 5200, noclip=true, font=(10,"Helvetica",:black),
dpi=700, name="./case_study_1_present_day_validation/plots/combined_scatter_fluxes_structures.png")


########## Flux~structure relationships 
valid_index = .!isnan.(vec(H_ref)) .&& .!isnan.(vec(NPP_ref)) .&& (vec(H_ref) .> 0 .&& vec(NPP_ref .> 0))
GMT.basemap(region=(1, 40, 0, 2000), figsize=(6, 5), xlabel="Height (m)", ylabel="NPP (g C m@+-2@+)",
theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,))
GMT.scatter!(vec(H_ref)[valid_index], vec(NPP_ref)[valid_index], fill=:black, alpha=95, legend=(label="Data", pos=:TL, box=:none,))
valid_index = .!isnan.(vec(TREED_output.H)) .&& .!isnan.(vec(TREED_output.NPP)) .&& (vec(TREED_output.H) .> 0 .&& vec(TREED_output.NPP .> 0))
GMT.scatter!(vec(TREED_output.H)[valid_index], vec(TREED_output.NPP)[valid_index], fill=:lightblue, alpha=98, legend="Model")

valid_index = .!isnan.(vec(H_ref)) .&& .!isnan.(vec(AGB_ref)) .&& (vec(H_ref) .> 0 .&& vec(AGB_ref .> 0))
GMT.basemap!(region=(1, 50, 0, 150), figsize=(6, 5), xlabel="Height (m)", ylabel="AGB (Mg C ha@+-1@+)",
    theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,), xshift = 7.25)
GMT.scatter!(vec(H_ref)[valid_index], vec(AGB_ref)[valid_index], fill=:black, alpha=95, legend=(label="Data", pos=:TL, box=:none,))
valid_index = .!isnan.(vec(TREED_output.H)) .&& .!isnan.(vec(AGB_model)) .&& (vec(TREED_output.H) .> 0 .&& vec(AGB_model .> 0))
GMT.scatter!(convert.(Float64, vec(TREED_output.H)[valid_index]), vec(AGB_model)[valid_index], fill=:lightblue, alpha=98, legend="Model") 
GMT.text!("(b)", x = -5.5, y = 152, noclip=true, font=(10,"Helvetica",:black))
GMT.text!("(a)", x = -68, y = 152, noclip=true, font=(10,"Helvetica",:black),
    dpi=700, name="./case_study_1_present_day_validation/plots/flux_structure_relationship.png")



########## Richness estimation 
FD_cpt = makecpt(cmap=:batlowK, continuous=true, range=(0, 0.25), overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="4/6/11", COLOR_FOREGROUND="250/204/249"))
grdimage(convert_raster_to_GMT_grid(TREED_output.functional_diversity), projection=:Mollweide, theme="A2xy",
    cmap=FD_cpt, xaxis=(annot=0, ), yaxis=(annot=60,), figsize=6.5, par=(FONT_ANNOT=7,))
colorbar!(pos=(paper=true, anchor=(3.5,-0.1), size=(4,0.2), justify=:TC, horizontal=true, triangles=:f), frame=(annot=:auto, ), xlabel="Functional diversity index (0-1)", 
    par=(FONT_LABEL=12, FONT_ANNOT_PRIMARY=10))

EH_cpt = makecpt(cmap=:batlowK, continuous=true, range=(0, 0.5), overrule_bg=true, par=(COLOR_NAN=235,COLOR_BACKGROUND="4/6/11", COLOR_FOREGROUND="250/204/249"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.gamma_EH), projection=:Mollweide, theme="A2xy",
    cmap=EH_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=6.5, par=(FONT_ANNOT=7,), xshift=7)
colorbar!(pos=(paper=true, anchor=(3.5,-0.1), size=(4,0.2), justify=:TC, horizontal=true, triangles=:f), frame=(annot=:auto, ), xlabel="Landscape heterogeneity (0-1)", 
    par=(FONT_LABEL=12, FONT_ANNOT_PRIMARY=10))

GI_cpt = makecpt(cmap=:batlowK, continuous=true, range=(0, 1.0), overrule_bg=true, par=(COLOR_NAN=235, COLOR_BACKGROUND="4/6/11", COLOR_FOREGROUND="250/204/249"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.gamma_GI), projection=:Mollweide, theme="A2xy",
    cmap=GI_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=6.5, par=(FONT_ANNOT=7,), xshift=7)
colorbar!(pos=(paper=true, anchor=(3.5,-0.1), size=(4,0.2), justify=:TC, horizontal=true, triangles=:f), frame=(annot=:auto, ), xlabel="Landscape fragmentation (0-1)", 
    par=(FONT_LABEL=12, FONT_ANNOT_PRIMARY=10))

SR_cpt = makecpt(cmap=:batlowK, continuous = true, 
    range=(0, 2000), overrule_bg=true, par=(COLOR_NAN=235, COLOR_BACKGROUND="4/6/11", COLOR_FOREGROUND="163.0/19.0/1.0"))

grdimage!(convert_raster_to_GMT_grid((10^3.8355) .* (TREED_output.diversity_index .^ 0.302779)), projection=:Mollweide, theme="A2xy",
    cmap=SR_cpt, xaxis=(annot=0,), yaxis=(annot=60,), figsize=10, par=(FONT_ANNOT=7,), xshift=-14, yshift=-7)
    text!("Y = 6847 × Diversity index@+0.30@+",frame=:none,region=(0,10,0,10), proj=:linear, x=5, y=5.3, noclip=true ,font=(8,"Helvetica",:black)) 
    colorbar!(pos=(paper=true, anchor=(10,-0.2), size=(4,0.2), justify=:TC, horizontal=true, triangles=:f), frame=(annot=:auto, ), xlabel="Species richness", 
    par=(FONT_LABEL=12, FONT_ANNOT_PRIMARY=10))

grdimage!(convert_raster_to_GMT_grid(SR_ref), projection=:Mollweide, theme="A2xy",
    cmap=SR_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=10, par=(FONT_ANNOT=7,), xshift=10.5)

basemap!(projection=:linear, frame=(annot=:auto, ), region=(-5.5, -0.3, 2.0, 3.9), xaxis=(label="log@-10@- modelled diversity potential (0-1)",),
    yaxis=(label="Observed richness (log@-10@-(N) per 10@+4@+km@+2@+)",), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=7), figsize=(10, 5),
    xshift=-10, yshift=-6.5)
GMT.scatter!(log10.(vec(TREED_output.diversity_index)), log10.(vec(SR_ref)), fill=:black, alpha=95)
x = -6:0.01:-0.3
GMT.plot!(x, 3.8355 .+ 0.302779  .* x, linecolor=:darkred, lw=2)
GMT.text!("Y = 3.84 + 0.30 × X", x=-1.25, y=2.25, font=(7, "Helvetica", :darkred))
GMT.text!("Correlation coefficient = 0.81", x=-1.2, y=2.12, font=(7, "Helvetica", :black))
GMT.text!("(f)", x=-5.5, y=4.1, noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(d)", x=-5.5, y=6.3, noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(e)", x=0.5, y=6.3, noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(a)", x=-5.5, y=8.5, noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(b)", x=-1.75, y=8.5, noclip=true, font=(10, "Helvetica", :black))
GMT.text!("(c)", x=1.5, y=8.5, noclip=true, font=(10, "Helvetica", :black), 
dpi=700, name="./case_study_1_present_day_validation/plots/richness_estimation_validation.png")



# Plot additional structures 
C = makecpt(cmap=((0,59,71), (254,228,171)), T=[-0.1,0.5,1.1],overrule_bg=true, par=(COLOR_NAN=235, COLOR_BACKGROUND="56/102/149", COLOR_FOREGROUND="15/42/3"))
C.label = ["Deciduous","Evergreen"]
grdimage(convert_raster_to_GMT_grid(TREED_output.seasonality), projection=:Mollweide, theme="A2xy",
    cmap=C, xaxis=(annot=0,), yaxis=(annot=60,), figsize=10, par=(FONT_ANNOT=7,))
colorbar!(pos=(paper=true, anchor=(5,-0.2), size=(4,0.2), justify=:TC, horizontal=true),
    B=:none, equal_size=(gap=0.1,),par=(FONT_ANNOT=12,))

a_ll_cpt = makecpt(cmap=:cork, hinge=1, range=(0, 5), overrule_bg=true, par=(COLOR_NAN=235, COLOR_BACKGROUND="56/102/149", COLOR_FOREGROUND="15/42/3"))
grdimage!(convert_raster_to_GMT_grid(TREED_output.a_ll), projection=:Mollweide, theme="A2xy",
    cmap=a_ll_cpt, xaxis=(annot=0,), yaxis=(annot=0,), figsize=10, par=(FONT_ANNOT=7,), xshift=10.5)
colorbar!(pos=(paper=true, anchor=(5,-0.2), size=(4,0.2), justify=:TC, horizontal=true), xlabel="Leaf longevity (years)",
par=(FONT_ANNOT=12,))
text!("(a)",frame=:none,region=(0,10,0,10), proj=:linear, x=-7, y=5, noclip=true ,font=(10,"Helvetica",:black)) 
text!("(b)",frame=:none,region=(0,10,0,10), proj=:linear, x=-0.5, y=5, noclip=true ,font=(10,"Helvetica",:black))

data = DataFrame(TREED_output)
data = convert.(Float64, data)
lat = data.Y
a_ll = data.a_ll
phenology = data.seasonality
H = data.H
LMA  = 1 ./ (data.SLA .* 0.47)

#cpt = makecpt(range = (0, 5), continuous=true, color=:cork, hinge=1)
cpt = makecpt(cmap="darkblue,darkgreen", range=[0, 1, 5])
GMT.basemap!(region=(0, 86, 45, 225), figsize=(8, 6), xlabel="Absolute latitude (degrees)", ylabel="LMA (g / m@+2@+)",
    theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,), yshift = -8, xshift=-9.5)
GMT.scatter!(abs.(lat[.!isnan.(a_ll) .&& H .> 0]), LMA[.!isnan.(a_ll) .&& H .> 0],
    zcolor=a_ll[.!isnan.(a_ll) .&& H .> 0], color=cpt, alpha = 97)
#GMT.colorbar!(pos=(paper=true, anchor=(8.5,5), size=(4,0.2), justify=:TC), xlabel="Leaf longevity (years)",
#    par=(FONT_ANNOT=12, FONT_ANNOT_PRIMARY = 12))
GMT.text!("(c)",x=-10, y=220, noclip=true ,font=(10,"Helvetica",:black))


glopnet_lma = CSV.read("case_study_1_present_day_validation/present_day_validation_data/Glopnet-subset.csv", DataFrame)
lat = convert.(Float64, glopnet_lma.Latitude)
lat = abs.(lat)
LMA = convert.(Float64, glopnet_lma.LMA)
DE = glopnet_lma.DecEv

# Linear fits, complete data 
data = DataFrame(y = LMA[DE .== "E"], x = lat[DE .== "E"])
linear_fit_E = lm(@formula(y ~ x), data)
data = DataFrame(y = LMA[DE .== "D"], x = lat[DE .== "D"])
linear_fit_D = lm(@formula(y ~ x), data)

Plots.scatter(abs.(glopnet_lma.Latitude), glopnet_lma.LMA, ylim=(20, 300))

GMT.basemap!(region=(0, 86, 10, 400), figsize=(8, 6), xlabel="Absolute latitude (degrees)", ylabel="LMA (g / m@+2@+)",
    theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,), xshift=11)
GMT.scatter!(lat[DE .== "E"], LMA[DE .== "E"], alpha = 50, color=:darkgreen)
GMT.scatter!(lat[DE .== "D"], LMA[DE .== "D"], alpha = 50, color=:darkblue)
x = 0:1:90
y_evergreen = 150.712 .+ 1.26 .* x
y_deciduous = 101.75 .- 0.5 .* x
GMT.plot!(x, y_evergreen, linecolor=:darkgreen, lw=1.5, legend=(label="Longevity > 1", pos=:TR, box=:none))
GMT.plot!(x, y_deciduous, linecolor=:darkblue, lw=1.5, legend=(label="Longevity <= 1", pos=:TR, box=:none))
GMT.text!("(d)",x=-10.25, y=400, noclip=true ,font=(10,"Helvetica",:black), 
dpi=700, name="./case_study_1_present_day_validation/plots/phenology_all.png")



# Additional plot for methods: illstruate SLA~a_ll and k~a_ll relationship 
a_ll = 0.001:0.01:5
SLA_model = (2e-4) .* (1/0.4763) .* 10 .^(2.25 .- 0.5 .* log10.(a_ll .* 12))
SLA_LPJ_needle = (2e-4) .* (1/0.4763) .* 10 .^(2.08 .- 0.4 .* log10.(a_ll .* 12))
SLA_LPJ_broad = (2e-4) .* (1/0.4763) .* 10 .^(2.22 .- 0.4 .* log10.(a_ll .* 12))

GMT.basemap(region=(0.001, 4.5, 0.009, 0.042), figsize=(6, 5), xlabel="a@-ll@- (years)", ylabel="SLA (m@+2@+ / g C)",
theme=("A2xy"), par=(FONT_LABEL=7, FONT_ANNOT_PRIMARY=6,))
GMT.plot!(a_ll, SLA_LPJ_needle, linecolor=:darkgreen, lw=1.5, legend=(label="LPJ needle", pos=:TR, box=:none,))
GMT.plot!(a_ll, SLA_LPJ_broad, linecolor=:lightred, lw=1.5, legend=(label="LPJ broad", pos=:TR, box=:none,))
GMT.plot!(a_ll, SLA_model, linecolor=:black, lw=1.5, legend=(label="Model", pos=:TR, box=:none,), 
dpi=700, name="./case_study_1_present_day_validation/plots/all_SLA_relationship.png")


