using Pkg
Pkg.activate(joinpath(@__DIR__, "."))
Pkg.develop(path=joinpath(@__DIR__, "Dashboard"))
Pkg.instantiate()


using JuMP # building models
using DataStructures # using dictionaries with a default value
using Clp # solver for the JuMP model
using CSV # readin of CSV files
using DataFrames # data tables
using JSON3
using Dashboard
include(joinpath(@__DIR__, "colors.jl")) # colors for the plots

data_dir = joinpath(@__DIR__, "data")

### Read in of parameters ###
# We define our sets from the csv files
technologies = readcsv("technologies.csv", dir=data_dir).technology
fuels = readcsv("fuels.csv", dir=data_dir).fuel
hour = 1:120
n_hour = length(hour)
storages = readcsv("storages.csv", dir=data_dir).storage

### define readin for regions
regions = readcsv("regions.csv", dir=data_dir).region

# Also, we read our input parameters via csv files
Demand = readin("demand_regions.csv", default=0, dims=2, dir=data_dir)
OutputRatio = readin("outputratio.csv", dims=2, dir=data_dir)
InputRatio = readin("inputratio.csv", dims=2, dir=data_dir)
VariableCost = readin("variablecost.csv", dims=1, dir=data_dir)
InvestmentCost = readin("investmentcost.csv", dims=1, dir=data_dir)
EmissionRatio = readin("emissionratio.csv", dims=1, dir=data_dir)
DemandProfile = readin("demand_timeseries_regions.csv", default=1/n_hour, dims=3, dir=data_dir)
MaxCapacity = readin("maxcapacity.csv",default=999,dims=2, dir=data_dir)
TagDispatchableTechnology = readin("tag_dispatchabletechnology.csv",dims=1, dir=data_dir)
CapacityFactor = readin("capacity_factors_regions.csv",default=0, dims=3, dir=data_dir)
for t in technologies
    if TagDispatchableTechnology[t] > 0
        for h in hour
            for r in regions
                CapacityFactor[r,t,h] = 1
            end
        end
    end
end

InvestmentCostStorage = readin("investmentcoststorage.csv",dims=1, dir=data_dir)
E2PRatio = readin("e2pratio.csv",dims=1, dir=data_dir)
StorageChargeEfficiency = readin("storagechargeefficiency.csv",dims=2, dir=data_dir)
StorageDisChargeEfficiency = readin("storagedischargeefficiency.csv",dims=2, dir=data_dir)
MaxStorageCapacity = readin("maxstoragecapacity.csv",default=999,dims=2, dir=data_dir)
StorageLosses = readin("storagelosses.csv",default=1,dims=2, dir=data_dir)

### Define your readin for MaxTradeCapacity
MaxTradeCapacity= readin("maxtradecapacity.csv",default=0,dims=3, dir=data_dir)
# our emission limit
EmissionLimit = 20000

# instantiate a model with an optimizer
ESM = Model(Clp.Optimizer)

# this creates our variables
@variable(ESM,TotalCost[technologies]>=0)
@variable(ESM,Production[regions,hour,technologies, fuels] >= 0)
@variable(ESM,Capacity[regions,technologies] >=0)
@variable(ESM,Use[regions,hour,technologies, fuels] >=0)
@variable(ESM,Emissions[regions,technologies] >=0)
@variable(ESM,Curtailment[regions,hour,fuels] >=0)

@variable(ESM,StorageEnergyCapacity[regions,s=storages,f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageCharge[regions,s=storages, hour, f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageDischarge[regions,s=storages, hour, f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageLevel[regions,s=storages, hour, f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,TotalStorageCost[storages] >= 0)


### Define your Export and Import variables
@variable(ESM,Export[hour,regions,regions,fuels] >= 0)
@variable(ESM,Import[hour,regions,regions,fuels] >= 0)

## constraints ##
# Generation must meet demand
@constraint(ESM, DemandAdequacy[r in regions,h in hour,f in fuels],
    sum(Production[r,h,t,f] for t in technologies) + sum(StorageDischarge[r,s,h,f] for s in storages if StorageDisChargeEfficiency[s,f]>0) + sum(Import[h,r,rr,f] for rr in regions)== 
        Demand[r,f]*DemandProfile[r,f,h] + sum(Use[r,h,t,f] for t in technologies)+Curtailment[r,h,f] + sum(StorageCharge[r,s,h,f] for s in storages if StorageChargeEfficiency[s,f] > 0) + sum(Export[h,r,rr,f] for rr in regions)
)

# calculate the total cost
@constraint(ESM, ProductionCost[t in technologies],
    sum(Production[r,h,t,f] * VariableCost[t] for f in fuels, h in hour, r in regions) + sum(Capacity[r,t] * InvestmentCost[t] for r in regions) == TotalCost[t]
)

# limit the production by the installed capacity
@constraint(ESM, ProductionFuntion_disp[r in regions,h in hour, t in technologies, f in fuels;TagDispatchableTechnology[t]>0],
    OutputRatio[t,f] * Capacity[r,t] * CapacityFactor[r,t,h] >= Production[r,h,t,f]
)
# for variable renewables, the production needs to be always at maximum
@constraint(ESM, ProductionFunction_res[r in regions,h in hour, t in technologies, f in fuels; TagDispatchableTechnology[t]==0], 
    OutputRatio[t,f] * Capacity[r,t] * CapacityFactor[r,t,h] == Production[r,h,t,f]
)

# define the use by the production
@constraint(ESM, UseFunction[r in regions,h in hour,t in technologies, f in fuels],
    InputRatio[t,f] * sum(Production[r,h,t,ff]/OutputRatio[t,ff] for ff in fuels if OutputRatio[t,ff]>0) == Use[r,h,t,f]
)

# define the emissions
@constraint(ESM, TechnologyEmissions[t in technologies, r in regions],
    sum(Production[r,h,t,f] for f in fuels, h in hour) * EmissionRatio[t] == Emissions[r,t]
)

# limit the emissions
@constraint(ESM, TotalEmissionsFunction,
    sum(Emissions[r,t] for t in technologies, r in regions) <= EmissionLimit
)

# installed capacity is limited by the maximum capacity
@constraint(ESM, MaxCapacityFunction[r in regions,t in technologies],
     Capacity[r,t] <= MaxCapacity[r,t]
)

# storage charge is limited by storage energy capacity and E2PRatio
@constraint(ESM, StorageChargeFunction[r in regions, s in storages, h in hour, f in fuels; StorageDisChargeEfficiency[s,f]>0], 
    StorageCharge[r,s,h,f] <= StorageEnergyCapacity[r,s,f]/E2PRatio[s]
)

# storage discharge is limited by storage energy capacity and E2PRatio
@constraint(ESM, StorageDischargeFunction[r in regions,s in storages, h in hour, f in fuels; StorageDisChargeEfficiency[s,f]>0], 
    StorageDischarge[r,s,h,f] <= StorageEnergyCapacity[r,s,f]/E2PRatio[s]
)

# storage level depends on previous period's storage level and current period charge/discharge
@constraint(ESM, StorageLevelFunction[r in regions,s in storages, h in hour, f in fuels; h>1 && StorageDisChargeEfficiency[s,f]>0], 
    StorageLevel[r,s,h,f] == StorageLevel[r,s,h-1,f]*StorageLosses[s,f] + StorageCharge[r,s,h,f]*StorageChargeEfficiency[s,f] - StorageDischarge[r,s,h,f]/StorageDisChargeEfficiency[s,f]
)

# storage level for first period does not depend on previous level but we set it to 50% energy capacity
@constraint(ESM, StorageLevelStartFunction[r in regions,s in storages, h in hour, f in fuels; h==1 && StorageDisChargeEfficiency[s,f]>0], 
    StorageLevel[r,s,h,f] == 0.5*StorageEnergyCapacity[r,s,f]*StorageLosses[s,f] + StorageCharge[r,s,h,f]*StorageChargeEfficiency[s,f] - StorageDischarge[r,s,h,f]/StorageDisChargeEfficiency[s,f]
)

# storage level is limited by storage capacity
@constraint(ESM, MaxStorageLevelFunction[r in regions,s in storages, h in hour, f in fuels; StorageDisChargeEfficiency[s,f]>0], 
    StorageLevel[r,s,h,f] <= StorageEnergyCapacity[r,s,f]
)

# storage cost are the sum of all storage technology costs
@constraint(ESM, StorageCostFunction[s in storages], 
    TotalStorageCost[s] == sum(StorageEnergyCapacity[r,s,f]*InvestmentCostStorage[s] for f in fuels, r in regions if StorageDisChargeEfficiency[s,f]>0)
)

# storage level at the end of a year has to equal storage level at the beginning of year
@constraint(ESM, StorageAnnualBalanceFunction[r in regions,s in storages, f in fuels; StorageDisChargeEfficiency[s,f]>0], 
    StorageLevel[r,s,n_hour,f] == 0.5*StorageEnergyCapacity[r,s,f]
)

# storage capacity is limited by max storage capacity
@constraint(ESM, StorageMaxCapacityConstraint[r in regions,s in storages], 
    sum(StorageEnergyCapacity[r,s,f] for f in fuels if StorageDisChargeEfficiency[s,f]>0) <= MaxStorageCapacity[r,s]
)

### write your equations for import/export
@constraint(ESM, ImportExportBalance[r in regions,rr in regions, f in fuels, h in hour], 
    Export[h,r,rr,f] == Import[h,rr,r,f]
)

@constraint(ESM, MaxTradeCapacityFunction[r in regions,rr in regions, f in fuels, h in hour], 
    Export[h,r,rr,f] <= MaxTradeCapacity[r,rr,f]
)

# the objective function
# total costs should be minimized
@objective(ESM, Min,
    sum(TotalCost[t] for t in technologies)
    + sum(TotalStorageCost[s] for s in storages)
    + sum(Export[h,r,rr,f]*0.1 for h in hour, r in regions, rr in regions, f in fuels)
)

# this starts the optimization
# the assigned solver (here Clp) will takes care of the solution algorithm
optimize!(ESM)
# reading our objective value
objective_value(ESM)

# some result analysis
value.(Production)
value.(Capacity)
value.(StorageEnergyCapacity)
value.(StorageDischarge)
value.(StorageLevel)
value.(StorageCharge)
value.(TotalStorageCost)

df_production = DataFrame(Containers.rowtable(value,Production; header = [:Region, :Hour, :Technology, :Fuel, :value]))
df_use = DataFrame(Containers.rowtable(value,Use; header = [:Region, :Hour, :Technology, :Fuel, :value]))
df_capacity = DataFrame(Containers.rowtable(value,Capacity; header = [:Region, :Technology, :value]))

df_storage_production = DataFrame(Containers.rowtable(value,StorageDischarge; header = [:Region, :Technology, :Hour, :Fuel, :value]))
df_storage_charge = DataFrame(Containers.rowtable(value,StorageCharge; header = [:Region, :Technology, :Hour, :Fuel, :value]))
df_storage_level = DataFrame(Containers.rowtable(value,StorageLevel; header = [:Region, :Technology, :Hour, :Fuel, :value]))

df_demand = DataFrame(
    (Region=r, Hour=h, Fuel=f, value=Demand[r,f]*DemandProfile[r,f,h]) for r in regions, f in fuels, h in hour
)

df_export = DataFrame(Containers.rowtable(value,Export; header = [:Hour, :From, :To, :Fuel, :value]))
df_import = DataFrame(Containers.rowtable(value,Import; header = [:Hour, :To, :From, :Fuel, :value]))

append!(df_use, df_storage_charge)
append!(df_production, df_storage_production)


# Define the path to the results directory
result_path = mkpath(joinpath(@__DIR__, "results"))
CSV.write(joinpath(result_path, "production.csv"), df_production)
CSV.write(joinpath(result_path, "use.csv"), df_use)
CSV.write(joinpath(result_path, "demand.csv"), df_demand)
CSV.write(joinpath(result_path, "capacity.csv"), df_capacity)
CSV.write(joinpath(result_path, "level.csv"), df_storage_level)
CSV.write(joinpath(result_path, "ex_import.csv"), df_import)
CSV.write(joinpath(result_path, "ex_export.csv"), df_export)

open(joinpath(result_path, "colors.json"), "w") do f
    JSON3.pretty(f, JSON3.write(colors))
    println(f)
end