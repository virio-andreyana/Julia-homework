using Pkg
Pkg.activate(joinpath(@__DIR__, "."))
# Pkg.develop(path=joinpath(@__DIR__, "Dashboard"))

using Dashboard

result_path = joinpath(@__DIR__, "results")
result = Result(result_path)
dashboard(result)