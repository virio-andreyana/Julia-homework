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
residual_year = 1990:10:2050 # NEW: residual capacity have its own lifetime. See residual lifetime
first_year = year[1]

### define readin for regions
regions = readcsv("regions.csv", dir=data_dir).region

# Also, we read our input parameters via csv files
Demand = readin("demand.csv", default=0, dims=3, dir=data_dir)
OutputRatio = readin("outputratio.csv", dims=2, dir=data_dir)
InputRatio = readin("inputratio.csv", dims=2, dir=data_dir)
VariableCost = readin("variablecost.csv", dims=2, dir=data_dir)
InvestmentCost = readin("investmentcost.csv", dims=2, dir=data_dir)
EmissionRatio = readin("emissionratio.csv", dims=1, dir=data_dir)
SequesterRatio = readin("sequesterratio.csv", dims=1, dir=data_dir) # NEW: add the amount of co2 a technology can sequester
DemandProfile = readin("demand_timeseries_regions.csv", default=1/n_hour, dims=3, dir=data_dir)
MaxCapacity = readin("maxcapacity.csv",default=999,dims=3, dir=data_dir)
TagDispatchableTechnology = readin("tag_dispatchabletechnology.csv",dims=1, dir=data_dir)
TagCCTechnology = readin("tag_cctechnology.csv",dims=1, dir=data_dir) # NEW: tag for CC technologies is added
TagCCAbleTechnology = readin("tag_cc_abletechnology.csv",dims=1, dir=data_dir) # NEW: tag for CC technologies is added
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
ResidualCapacityData = readin("residualcapacitylifetime.csv",default=0.0,dims=3, dir=data_dir) # Edited: now ResidualCapacity has its own lifetime
TechnologyLifetime = readin("technologylifetime.csv", default=10,dims=1, dir=data_dir)

### Define your readin for MaxTradeCapacity
MaxTradeCapacity = readin("maxtradecapacity.csv",default=0,dims=4, dir=data_dir)
TradeDistance = readin("tradedistance.csv",default=0,dims=2, dir=data_dir)
TradeCostFactor = readin("tradecostfactor.csv",default=0,dims=1, dir=data_dir)
TradeLossFactor = readin("tradelossfactor.csv",default=0,dims=1, dir=data_dir)

# our emission limit
AnnualEmissionLimit = readin("annualemissionlimit.csv",default=99999,dims=1, dir=data_dir)
#AnnualEmissionLimit = DefaultDict(99999,)
CO2L = 1.0
TotalEmissionLimit = 400000 * CO2L
DiscountRate = 0.05

# our sequester limit
TotalSequesterLimit = readin("totalsequesterlimit.csv",default=99999,dims=1, dir=data_dir)

techCC = Dict("GasPowerPlantCC"=>"GasPowerPlant", "CoalPowerPlantCC"=>"CoalPowerPlant", "GasCHPPlantCC"=>"GasCHPPlant", "CoalCHPPlantCC"=>"CoalCHPPlant")

# Adjust the share of CO2 captured by the CC powerplants. Value ranges between 0 to 1. The higher the CC_share, the higher CO2 captured as 'fuel' and the lower is it's emission.
function adjust_CO2_capture(CC_share)
    for t in keys(techCC)
        SequesterRatio[t] = round(CC_share*EmissionRatio[techCC[t]]; digits = 2)
        EmissionRatio[t] = round((1-CC_share)*EmissionRatio[techCC[t]]; digits = 2)
    end
    return EmissionRatio, SequesterRatio
end

# Adjust the investement cost of CC powerplants in relation to normal powerplants. Logically it will be more expensive, hence cost_diff > 1.
function adjust_investmentcost(cost_diff)
    for y in year
        for t in keys(techCC)
            InvestmentCost[y,t] = cost_diff*InvestmentCost[y,techCC[t]]
        end
    end
    return InvestmentCost
end

function extend_year_scenario(year_ref, year_add)
    for y in year_add, r in regions
        AnnualEmissionLimit[y] = AnnualEmissionLimit[year_ref]
        for t in technologies
            InvestmentCost[y,t] = InvestmentCost[year_ref,t]
            VariableCost[y,t] = VariableCost[year_ref,t]
            MaxCapacity[y,r,t] = MaxCapacity[year_ref,r,t]
        end
        for s in storages
            InvestmentCostStorage[y,s] = InvestmentCostStorage[year_ref,s]
            MaxStorageCapacity[y,r,s] = MaxStorageCapacity[year_ref,r,s]
        end
        for f in fuels
            Demand[y,r,f] = Demand[year_ref,r,f]
            for rr in regions
                MaxTradeCapacity[y,r,rr,f] = MaxTradeCapacity[year_ref,r,rr,f]
            end
        end
    end
end

function limit_regional_sequestering(seq_share)
    for r in regions
        TotalSequesterLimit[r] = TotalSequesterLimit[r] * seq_share
    end
    return TotalSequesterLimit 
end

# CC Scenario options
EmissionRatio, SequesterRatio = adjust_CO2_capture(0.9)
InvestmentCost = adjust_investmentcost(1.4)
SequesterCost = 0.26 # Cost of sequestering per Ton of CO2
# TotalSequesterLimit = limit_regional_sequestering(0.5)

#year_add = 2060:10:2100
#extend_year_scenario(2050,year_add)
#year = minimum(year):10:maximum(year_add)

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
@variable(ESM,ResidualCapacity[residual_year,regions,t=technologies] >=0) #NEW: Residual capacities are now a variable
@variable(ESM,TotalCapacity[year,regions,technologies] >=0)
@variable(ESM,Use[year,regions,hour,technologies, fuels] >=0)
@variable(ESM,AnnualEmissions[year,regions,technologies])
@variable(ESM,AnnualSequester[year,regions,technologies])
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
@variable(ESM,TotalSequester[regions] >= 0)
@variable(ESM,TotalSequesterCost >= 0)


## constraints ##
# Generation must meet demand
@constraint(ESM, DemandAdequacy[y in year, r in regions,h in hour,f in fuels],
    sum(Production[y,r,h,t,f] for t in technologies) + sum(StorageDischarge[y,r,s,h,f] for s in storages if StorageDisChargeEfficiency[s,f]>0) + sum(Import[y,h,r,rr,f] for rr in regions) == 
        Demand[y,r,f]*DemandProfile[r,f,h] + sum(Use[y,r,h,t,f] for t in technologies)+Curtailment[y,r,h,f] + sum(StorageCharge[y,r,s,h,f] for s in storages if StorageChargeEfficiency[s,f] > 0) + sum(Export[y,h,r,rr,f] for rr in regions)
)


# EDITED: calculate the total cost NOTE: Retrofiting residual capacity for CC technology has a simplified investment cost. Calculating when does the retroffiting occur is too complex.
@constraint(ESM, ProductionCost[t in technologies],
    sum(Production[y,r,h,t,f] * VariableCost[y,t] * YearlyDifferenceMultiplier[y] / (1+DiscountRate)^(y - minimum(year)) for f in fuels, h in hour, r in regions, y in year) 
    + sum(NewCapacity[y,r,t] * InvestmentCost[y,t] / (1+DiscountRate)^(y - minimum(year)) for r in regions, y in year)
    + sum(ResidualCapacity[y_res,r,t] * (InvestmentCost[minimum(year),t] - InvestmentCost[minimum(year),techCC[t]]) for r in regions, y_res in residual_year if TagCCTechnology[t]>0)
    == TotalCost[t]
)

# calculate the total installed capacity in each year
@constraint(ESM, TotalCapacityFunction[y in year, t in technologies, r in regions], 
    sum(NewCapacity[yy,r,t] for yy in year if yy <= y &&  yy + TechnologyLifetime[t] >= y)
    + sum(ResidualCapacity[y_res,r,t] for y_res in residual_year if y_res <= y &&  y_res + TechnologyLifetime[t] >= y)
    == TotalCapacity[y,r,t] ####################### &&  yy + TechnologyLifetime[t] > y // ; (y-1) + TechnologyLifetime[t] >= y
)

# NEW: add the possibility to convert residual capacity of fossil fuel powerplant into CC powerplant
@constraint(ESM, RetrofitResidualCapacityFunction[y_res in residual_year, t in technologies, r in regions; TagCCTechnology[t]>0], 
    0.9*ResidualCapacityData[y_res,r,techCC[t]] == ResidualCapacity[y_res,r,t] + ResidualCapacity[y_res,r,techCC[t]]
)

# NEW: other residual capacity that are not connected with techCC remain the same
@constraint(ESM, ResidualCapacityFunction[y_res in residual_year, t in technologies, r in regions; TagCCAbleTechnology[t]==0], 
    0.9*ResidualCapacityData[y_res,r,t] == ResidualCapacity[y_res,r,t]
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

# NEW: define the technology sequester emission
@constraint(ESM, AnnualTechnologySequester[y in year,t in technologies, r in regions],
    sum(Production[y,r,h,t,f] for f in fuels, h in hour) * SequesterRatio[t] == AnnualSequester[y,r,t]
)
# limit the emissions per year
@constraint(ESM, AnnualEmissionsLimitFunction[y in year],
    sum(AnnualEmissions[y,r,t] for t in technologies, r in regions) <= AnnualEmissionLimit[y]
)

# account for the total emissions
@constraint(ESM, TotalEmissionsAccounting,
    sum(Production[y,r,h,t,f] * EmissionRatio[t] * YearlyDifferenceMultiplier[y] for f in fuels, h in hour,y in year,t in technologies,r in regions)  == TotalEmissions
)

# NEW: account for the total sequester emission
@constraint(ESM, TotalSequesterAccounting[r in regions],
    sum(Production[y,r,h,t,f] * SequesterRatio[t] * YearlyDifferenceMultiplier[y] for f in fuels, h in hour,y in year,t in technologies)  == TotalSequester[r]
)

# limit the total emissions
@constraint(ESM, TotalEmissionsLimitFunction,
    TotalEmissions <= TotalEmissionLimit
)

# NEW: limit the total sequester by region
@constraint(ESM, TotalSequesterLimitFunction[r in regions],
    TotalSequester[r] <= TotalSequesterLimit[r]
)

# NEW: sequestering co2 has cost per ton. the value 0f 0.26 is estimated by Fabio
@constraint(ESM, TotalSequesterCostFunction,
    sum(TotalSequester[r] for r in regions) * SequesterCost == TotalSequesterCost
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
    + TotalSequesterCost
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
#sum(value.(SalvageValue))
value.(TotalSequester)

df_production = DataFrame(Containers.rowtable(value, Production; header = [:Year, :Region, :Hour, :Technology, :Fuel, :value]))
df_use = DataFrame(Containers.rowtable(value, Use; header = [:Year, :Region, :Hour, :Technology, :Fuel, :value]))
df_capacity = DataFrame(Containers.rowtable(value, TotalCapacity; header = [:Year, :Region, :Technology, :value]))
df_newcapacity = DataFrame(Containers.rowtable(value, NewCapacity; header = [:Year, :Region, :Technology, :value]))
df_annualemissions = filter(row -> row.value != 0, DataFrame(Containers.rowtable(value,AnnualEmissions; header = [:Year, :Region,:Technology, :value])))
df_annualsequester = filter(row -> row.value != 0, DataFrame(Containers.rowtable(value,AnnualSequester; header = [:Year, :Region,:Technology, :value]))) # NEW: add sequester
df_residualcapacity = DataFrame(Containers.rowtable(value, ResidualCapacity; header = [:Year, :Region, :Technology, :value]))

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
result_path = mkpath(joinpath(@__DIR__, "results NoCC"))
CSV.write(joinpath(result_path, "production.csv"), df_production)
CSV.write(joinpath(result_path, "use.csv"), df_use)
CSV.write(joinpath(result_path, "demand.csv"), df_demand)
CSV.write(joinpath(result_path, "capacity.csv"), df_capacity)
CSV.write(joinpath(result_path, "level.csv"), df_storage_level)
CSV.write(joinpath(result_path, "ex_import.csv"), df_import)
CSV.write(joinpath(result_path, "ex_export.csv"), df_export)
CSV.write(joinpath(result_path, "emission.csv"), df_annualemissions)
CSV.write(joinpath(result_path, "sequester.csv"), df_annualsequester)
CSV.write(joinpath(result_path, "newcapacity.csv"), df_newcapacity)
CSV.write(joinpath(result_path, "residualcapacity.csv"), df_residualcapacity)

open(joinpath(result_path, "colors.json"), "w") do f
    JSON3.pretty(f, JSON3.write(colors))
    println(f)
end

function extract_string_before_bracket(s::AbstractString)
    parts = split(s, "[")
    return parts[1]
end

function extract_string_before_spaces(s::AbstractString)
    parts = split(s, " ")
    return parts[1]
end

function binding_constraints(model, threshold=1e-8)
    df = DataFrame(con=String[], binding=String[])
    for (F, S) in list_of_constraint_types(model)
        for con in all_constraints(model, F, S)
            if abs(dual(con)) > threshold
                push!(df,[string(con)," Binding"])
            else
                push!(df,[string(con)," Non-binding"])
            end
        end
    end

    df.con = map(extract_string_before_bracket, df.con)
    df.con = map(extract_string_before_spaces, df.con)
    df[!, :con_binding] = df.con .* df.binding
    grouped_df = groupby(df, :con_binding)

    df_counts = combine(grouped_df, nrow => :count)
    split_columns = split.(df_counts.con_binding, " ")
    df_counts.constraint = getindex.(split_columns, 1)
    df_counts.binding = getindex.(split_columns, 2)
    select!(df_counts, Not(:con_binding))

    return df_counts[:, [:constraint, :binding, :count]]
end

df_binding = binding_constraints(ESM)

CSV.write(joinpath(result_path, "binding.csv"), df_binding)

file = open(joinpath(result_path, "objective.txt"), "w")
println(file,objective_value(ESM))
close(file)

