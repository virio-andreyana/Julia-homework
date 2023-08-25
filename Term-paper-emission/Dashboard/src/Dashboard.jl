module Dashboard

using PlotlyJS
using Dash
using DataFrames
using DataFramesMeta
using CSV
using DataStructures
using UnPack
using JSON3
using Blink
using Base.Threads

include("dash_functions.jl")
include("plotting_funcs.jl")
include("sankey.jl")
include("helper_functions.jl")

export dashboard, Result, readcsv, readin

end # module Dashboard
