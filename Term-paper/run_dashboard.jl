using Pkg
Pkg.activate(joinpath(@__DIR__, "."))
# Pkg.develop(path=joinpath(@__DIR__, "Dashboard"))

using Dashboard

file_path = "test"

result_path = joinpath(@__DIR__, "results\\$file_path")
result = Result(result_path)
dashboard(result)