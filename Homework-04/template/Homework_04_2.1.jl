using Pkg
Pkg.activate(joinpath(@__DIR__, "."))
Pkg.develop(path=joinpath(@__DIR__, "Dashboard"))
Pkg.instantiate()


using JuMP # building models
using DataStructures # using dictionaries with a default value
using HiGHS # solver for the JuMP model
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
year = 2020:10:2050

### define readin for regions
regions = readcsv("regions.csv", dir=data_dir).region

# Also, we read our input parameters via csv files
Demand = readin("demand.csv", default=0, dims=3, dir=data_dir)
OutputRatio = readin("outputratio.csv", dims=2, dir=data_dir)
InputRatio = readin("inputratio.csv", dims=2, dir=data_dir)
VariableCost = readin("variablecost.csv", dims=2, dir=data_dir)
InvestmentCost = readin("investmentcost.csv", dims=2, dir=data_dir)
EmissionRatio = readin("emissionratio.csv", dims=1, dir=data_dir)
DemandProfile = readin("demand_timeseries_regions.csv", default=1/n_hour, dims=3, dir=data_dir)
MaxCapacity = readin("maxcapacity.csv",default=999,dims=3, dir=data_dir)
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

InvestmentCostStorage = readin("investmentcoststorage.csv",dims=2, dir=data_dir)
E2PRatio = readin("e2pratio.csv",dims=1, dir=data_dir)
StorageChargeEfficiency = readin("storagechargeefficiency.csv",dims=2, dir=data_dir)
StorageDisChargeEfficiency = readin("storagedischargeefficiency.csv",dims=2, dir=data_dir)
MaxStorageCapacity = readin("maxstoragecapacity.csv",default=9999,dims=3, dir=data_dir)
StorageLosses = readin("storagelosses.csv",default=1,dims=2, dir=data_dir)
ResidualCapacity = readin("residualcapacity.csv",default=0,dims=3, dir=data_dir)
TechnologyLifetime = readin("technologylifetime.csv", default=10,dims=1, dir=data_dir)

### Define your readin for MaxTradeCapacity
MaxTradeCapacity = readin("maxtradecapacity.csv",default=0,dims=4, dir=data_dir)
TradeDistance = readin("tradedistance.csv",default=0,dims=2, dir=data_dir)
TradeCostFactor = readin("tradecostfactor.csv",default=0,dims=1, dir=data_dir)
TradeLossFactor = readin("tradelossfactor.csv",default=0,dims=1, dir=data_dir)

# our emission limit
AnnualEmissionLimit = readin("annualemissionlimit.csv",default=99999,dims=1, dir=data_dir)
#AnnualEmissionLimit = DefaultDict(99999,)
TotalEmissionLimit = 400000
DiscountRate = 0.05



# create a multiplier to weight the different years correctly
YearlyDifferenceMultiplier = Dict()
for i in 1:length(year)-1
    difference = year[i+1] - year[i]
    # Store the difference in the dictionary
    YearlyDifferenceMultiplier[year[i]] = difference
end
YearlyDifferenceMultiplier[year[end]] = 1
# this gives us the distance between each year for all years
YearlyDifferenceMultiplier

# instantiate a model with an optimizer
ESM = Model(HiGHS.Optimizer)

# this creates our variables
@variable(ESM,TotalCost[technologies]>=0)
@variable(ESM,Production[year,regions,hour,technologies, fuels] >= 0)
@variable(ESM,NewCapacity[year,regions,technologies] >=0)
@variable(ESM,TotalCapacity[year,regions,technologies] >=0)
@variable(ESM,Use[year,regions,hour,technologies, fuels] >=0)
@variable(ESM,AnnualEmissions[year,regions,technologies])
@variable(ESM,Curtailment[year,regions,hour,fuels] >=0)

@variable(ESM,TotalStorageEnergyCapacity[year,regions,s=storages,f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,NewStorageEnergyCapacity[year,regions,s=storages,f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageCharge[year,regions,s=storages, hour, f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageDischarge[year,regions,s=storages, hour, f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageLevel[year,regions,s=storages, hour, f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,TotalStorageCost[storages] >= 0)

@variable(ESM,Import[year,hour,regions,regions,fuels] >= 0)
@variable(ESM,Export[year,hour,regions,regions,fuels] >= 0)

@variable(ESM,TotalEmissions >= 0)


## constraints ##
# Generation must meet demand
@constraint(ESM, DemandAdequacy[y in year, r in regions,h in hour,f in fuels],
    sum(Production[y,r,h,t,f] for t in technologies) + sum(StorageDischarge[y,r,s,h,f] for s in storages if StorageDisChargeEfficiency[s,f]>0) + sum(Import[y,h,r,rr,f] for rr in regions) == 
        Demand[y,r,f]*DemandProfile[r,f,h] + sum(Use[y,r,h,t,f] for t in technologies)+Curtailment[y,r,h,f] + sum(StorageCharge[y,r,s,h,f] for s in storages if StorageChargeEfficiency[s,f] > 0) + sum(Export[y,h,r,rr,f] for rr in regions)
)

# calculate the total cost
@constraint(ESM, ProductionCost[t in technologies],
    sum(Production[y,r,h,t,f] * VariableCost[y,t] * YearlyDifferenceMultiplier[y] / (1+DiscountRate)^(y - minimum(year)) for f in fuels, h in hour, r in regions, y in year) 
    + sum(NewCapacity[y,r,t] * InvestmentCost[y,t] / (1+DiscountRate)^(y - minimum(year)) for r in regions, y in year)
    == TotalCost[t]
)

# calculate the total installed capacity in each year
@constraint(ESM, TotalCapacityFunction[y in year, t in technologies, r in regions], 
    sum(NewCapacity[yy,r,t] for yy in year if yy <= y &&  yy + TechnologyLifetime[t] >= y) + ResidualCapacity[y,r,t] == TotalCapacity[y,r,t] ####################### &&  yy + TechnologyLifetime[t] > y // ; (y-1) + TechnologyLifetime[t] >= y
)

# limit the production by the installed capacity
@constraint(ESM, ProductionFuntion_disp[y in year,r in regions,h in hour, t in technologies, f in fuels;TagDispatchableTechnology[t]>0],
    OutputRatio[t,f] * TotalCapacity[y,r,t] * CapacityFactor[r,t,h] >= Production[y,r,h,t,f]
)
# for variable renewables, the production needs to be always at maximum
@constraint(ESM, ProductionFunction_res[y in year,r in regions,h in hour, t in technologies, f in fuels; TagDispatchableTechnology[t]==0], 
    OutputRatio[t,f] * TotalCapacity[y,r,t] * CapacityFactor[r,t,h] == Production[y,r,h,t,f]
)

# define the use by the production
@constraint(ESM, UseFunction[y in year,r in regions,h in hour,t in technologies, f in fuels],
    InputRatio[t,f] * sum(Production[y,r,h,t,ff] for ff in fuels) == Use[y,r,h,t,f]
)

# define the technology emissions
@constraint(ESM, AnnualTechnologyEmissions[y in year,t in technologies, r in regions],
    sum(Production[y,r,h,t,f] for f in fuels, h in hour) * EmissionRatio[t] == AnnualEmissions[y,r,t]
)

# limit the emissions per year
@constraint(ESM, AnnualEmissionsLimitFunction[y in year],
    sum(AnnualEmissions[y,r,t] for t in technologies, r in regions) <= AnnualEmissionLimit[y]
)

# account for the total emissions
@constraint(ESM, TotalEmissionsAccounting,
    sum(Production[y,r,h,t,f] * EmissionRatio[t] * YearlyDifferenceMultiplier[y] for f in fuels, h in hour,y in year,t in technologies,r in regions)  == TotalEmissions
)
# limit the total emissions
@constraint(ESM, TotalEmissionsLimitFunction,
    TotalEmissions <= TotalEmissionLimit
)


# installed capacity is limited by the maximum capacity
@constraint(ESM, MaxCapacityFunction[y in year, r in regions,t in technologies],
     TotalCapacity[y,r,t]  <= MaxCapacity[y,r,t]
)

# storage charge is limited by storage energy capacity and E2PRatio
@constraint(ESM, StorageChargeFunction[y in year,r in regions, s in storages, h in hour, f in fuels; StorageDisChargeEfficiency[s,f]>0], 
    StorageCharge[y,r,s,h,f] <= TotalStorageEnergyCapacity[y,r,s,f]/E2PRatio[s]
)

# account for currently installed storage capacities
@constraint(ESM, TotalStorageCapacityFunction[y in year, s in storages, r in regions, f in fuels; StorageDisChargeEfficiency[s,f]>0],
    sum(NewStorageEnergyCapacity[yy,r,s,f] for yy in year if yy<=y) == TotalStorageEnergyCapacity[y,r,s,f]
)


# storage discharge is limited by storage energy capacity and E2PRatio
@constraint(ESM, StorageDischargeFunction[y in year,r in regions,s in storages, h in hour, f in fuels; StorageDisChargeEfficiency[s,f]>0], 
    StorageDischarge[y,r,s,h,f] <= TotalStorageEnergyCapacity[y,r,s,f]/E2PRatio[s]
)

# storage level depends on previous period's storage level and current period charge/discharge
@constraint(ESM, StorageLevelFunction[y in year,r in regions,s in storages, h in hour, f in fuels; h>1 && StorageDisChargeEfficiency[s,f]>0], 
    StorageLevel[y,r,s,h,f] == StorageLevel[y,r,s,h-1,f]*StorageLosses[s,f] + StorageCharge[y,r,s,h,f]*StorageChargeEfficiency[s,f] - StorageDischarge[y,r,s,h,f]/StorageDisChargeEfficiency[s,f]
)

# storage level for first period does not depend on previous level but we set it to 50% energy capacity
@constraint(ESM, StorageLevelStartFunction[y in year,r in regions,s in storages, h in hour, f in fuels; h==1 && StorageDisChargeEfficiency[s,f]>0], 
    StorageLevel[y,r,s,h,f] == 0.5*TotalStorageEnergyCapacity[y,r,s,f]*StorageLosses[s,f] + StorageCharge[y,r,s,h,f]*StorageChargeEfficiency[s,f] - StorageDischarge[y,r,s,h,f]/StorageDisChargeEfficiency[s,f]
)

# storage level is limited by storage capacity
@constraint(ESM, MaxStorageLevelFunction[y in year,r in regions,s in storages, h in hour, f in fuels; StorageDisChargeEfficiency[s,f]>0], 
    StorageLevel[y,r,s,h,f] <= TotalStorageEnergyCapacity[y,r,s,f]
)

# storage cost are the sum of all storage technology costs
@constraint(ESM, StorageCostFunction[s in storages], 
    TotalStorageCost[s] == 
    (sum(NewStorageEnergyCapacity[y,r,s,f]*InvestmentCostStorage[y,s] for f in fuels, r in regions, y in year if StorageDisChargeEfficiency[s,f]>0))
    / (1+DiscountRate)^(y - minimum(year))
)

# storage level at the end of a year has to equal storage level at the beginning of year
@constraint(ESM, StorageAnnualBalanceFunction[y in year,r in regions,s in storages, f in fuels; StorageDisChargeEfficiency[s,f]>0], 
    StorageLevel[y,r,s,n_hour,f] == 0.5*TotalStorageEnergyCapacity[y,r,s,f]
)

# storage capacity is limited by max storage capacity
@constraint(ESM, StorageMaxCapacityConstraint[y in year,r in regions,s in storages], 
    sum(TotalStorageEnergyCapacity[y,r,s,f] for f in fuels if StorageDisChargeEfficiency[s,f]>0) <= MaxStorageCapacity[y,r,s]
)

### write your equations for import/export
@constraint(ESM, ImportExportBalance[y in year,h in hour,r in regions, rr in regions, f in fuels],
    Export[y,h,r,rr,f]*(1-TradeLossFactor[f]*TradeDistance[r,rr]) == Import[y,h,rr,r,f]
)

@constraint(ESM, MaxImportFunction[y in year,h in hour,r in regions, rr in regions, f in fuels],
    Import[y,h,r,rr,f] <= MaxTradeCapacity[y,r,rr,f]
)

# the objective function
# total costs should be minimized
@objective(ESM, Min,
    sum(TotalCost[t] for t in technologies)
    + sum(TotalStorageCost[s] for s in storages)
    + sum(Export[y,h,r,rr,f]*TradeCostFactor[f]*TradeDistance[r,rr]  * YearlyDifferenceMultiplier[y] / (1+DiscountRate)^(y - minimum(year)) for h in hour, r in regions, rr in regions, f in fuels, y in year)
)

# this starts the optimization
# the assigned solver (here Clp) will takes care of the solution algorithm
optimize!(ESM)
# reading our objective value
objective_value(ESM)

# some result analysis
value.(Production)
value.(TotalCapacity)
value.(NewStorageEnergyCapacity)
value.(TotalStorageEnergyCapacity)
value.(StorageDischarge)
value.(StorageLevel)
value.(StorageCharge)
value.(TotalStorageCost)
sum(value.(SalvageValue))

df_production = DataFrame(Containers.rowtable(value, Production; header = [:Year, :Region, :Hour, :Technology, :Fuel, :value]))
df_use = DataFrame(Containers.rowtable(value, Use; header = [:Year, :Region, :Hour, :Technology, :Fuel, :value]))
df_capacity = DataFrame(Containers.rowtable(value, TotalCapacity; header = [:Year, :Region, :Technology, :value]))
df_newcapacity = DataFrame(Containers.rowtable(value, NewCapacity; header = [:Year, :Region, :Technology, :value]))
df_annualemissions = filter(row -> row.value != 0, DataFrame(Containers.rowtable(value,AnnualEmissions; header = [:Year, :Region,:Technology, :value])))

df_storage_production = DataFrame(Containers.rowtable(value,StorageDischarge; header = [:Year, :Region, :Technology, :Hour, :Fuel, :value]))
df_storage_charge = DataFrame(Containers.rowtable(value,StorageCharge; header = [:Year, :Region, :Technology, :Hour, :Fuel, :value]))
df_storage_level = DataFrame(Containers.rowtable(value,StorageLevel; header = [:Year, :Region, :Technology, :Hour, :Fuel, :value]))

df_demand = DataFrame(
    (Year=y, Region=r, Hour=h, Fuel=f, value=Demand[y,r,f]*DemandProfile[r,f,h]) for y in year, r in regions, f in fuels, h in hour
)

df_export = DataFrame(Containers.rowtable(value,Export; header = [:Year, :Hour, :From, :To, :Fuel, :value]))
df_import = DataFrame(Containers.rowtable(value,Import; header = [:Year, :Hour, :To, :From, :Fuel, :value]))

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
CSV.write(joinpath(result_path, "emission.csv"), df_annualemissions)
CSV.write(joinpath(result_path, "emission.csv"), df_annualemissions)
CSV.write(joinpath(result_path, "newcapacity.csv"), df_newcapacity)

open(joinpath(result_path, "colors.json"), "w") do f
    JSON3.pretty(f, JSON3.write(colors))
    println(f)
end