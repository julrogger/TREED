module TREED

using Rasters, ArchGDAL, NCDatasets
using DimensionalData
using DimensionalData.Lookups
using Optimization
using OptimizationMetaheuristics
using CircularArrays
using Statistics
using Distributions
using Roots
using JLD2
using Images

include("./TREED_model_functions.jl")
include("./TREED_physiological_functions.jl")
include("./TREED_pars.jl")

export TREEDnonsteadycontinuous
export TREEDsteadystep
export TREEDsteadycontinuous
export TREEDnonsteadystep
include("./TREEDnonsteadycontinuous.jl")
include("./TREEDnonsteadystep.jl")
include("./TREEDsteadycontinuous.jl")
include("./TREEDsteadystep.jl")

export raster_area
export area_weighted_average
include("./TREED_auxiliary.jl")


end