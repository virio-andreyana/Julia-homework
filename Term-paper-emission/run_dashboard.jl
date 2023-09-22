using Pkg
Pkg.activate(joinpath(@__DIR__, "."))
# Pkg.develop(path=joinpath(@__DIR__, "Dashboard"))

result_path = joinpath(@__DIR__, "results CO2L1.0")
result = Result(result_path)
dashboard(result)



#Virio Blink Error solution
#readdir(joinpath(Sys.BINDIR, Base.LIBEXECDIR))
#cp(joinpath(Sys.BINDIR, Base.LIBEXECDIR, "julia", "7z.exe"), joinpath(Sys.BINDIR, Base.LIBEXECDIR, "7z.exe"))
#Blink.AtomShell.install()