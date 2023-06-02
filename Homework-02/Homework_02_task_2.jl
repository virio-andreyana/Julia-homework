using JuMP # building models
using DataStructures # using dictionaries with a default value
using Clp # solver for the JuMP model
using CSV # readin of CSV files
using DataFrames # data tables
using Statistics # mean function
using Plots  # generate graphs
using StatsPlots # additional features for plots
using Plots.Measures
include(joinpath(@__DIR__, "colors.jl")) 


# read the csv files
df_technologies = CSV.read("technologies.csv", DataFrame,stringtype=String)
df_fuels = CSV.read("fuels.csv", DataFrame,stringtype=String)
df_demand = CSV.read("demand.csv", DataFrame,stringtype=String)
df_inputratio = CSV.read("inputratio.csv", DataFrame,stringtype=String)
df_outputratio = CSV.read("outputratio.csv", DataFrame,stringtype=String)
df_investmentcost = CSV.read("investmentcost.csv", DataFrame,stringtype=String)
df_variablecost = CSV.read("variablecost.csv", DataFrame,stringtype=String)
df_emissionratio = CSV.read("emissionratio.csv", DataFrame,stringtype=String)
df_maxcapacity = CSV.read("maxcapacity.csv", DataFrame,stringtype=String)
df_tagdispatchabletechnology = CSV.read("tag_dispatchabletechnology.csv", DataFrame,stringtype=String)
df_demand_timeseries = CSV.read("demand_timeseries.csv", DataFrame,stringtype=String)
df_capacity_factor = CSV.read("capacity_factors_2018.csv", DataFrame,stringtype=String)

# storage csv files
df_storages = CSV.read("storages.csv", DataFrame,stringtype=String)
df_investmentcoststorage = CSV.read("investmentcoststorage.csv", DataFrame,stringtype=String)
df_e2pratio = CSV.read("e2pratio.csv", DataFrame,stringtype=String)
df_storagechargeefficiency = CSV.read("storagechargeefficiency.csv", DataFrame,stringtype=String)
df_storagedischargeefficiency = CSV.read("storagedischargeefficiency.csv", DataFrame,stringtype=String)
df_maxstoragecapacity = CSV.read("maxstoragecapacity.csv", DataFrame,stringtype=String)
df_storagelosses = CSV.read("storagelosses.csv", DataFrame,stringtype=String)


# readin function for parameters; this makes handling easier
readin(x; default=0,dims=1) = DefaultDict(default,Dict((dims > 1 ? Tuple(row[y] for y in 1:dims) : row[1]) => row[dims+1] for row in eachrow(x)))

# We define our sets from the csv files
technologies = df_technologies.technology
fuels = df_fuels.fuel
hour = 1:120
n_hour = length(hour)
storages = df_storages.storage


# Also, we read our input parameters via csv files
Demand = readin(df_demand,dims=1)
OutputRatio = readin(df_outputratio,dims=2)
InputRatio = readin(df_inputratio,dims=2)
VariableCost = readin(df_variablecost,dims=1)
InvestmentCost = readin(df_investmentcost,dims=1)
EmissionRatio = readin(df_emissionratio,dims=1)
TagDispatchableTechnology = readin(df_tagdispatchabletechnology,dims=1)
DemandProfile = readin(df_demand_timeseries,default=1/120,dims=2)
CapacityFactor = readin(df_capacity_factor,default=0,dims=2)
for t in technologies
    if TagDispatchableTechnology[t] > 0
        for h in hour
            CapacityFactor[t,h] = 1
        end
    end
end
MaxCapacity = readin(df_maxcapacity,default=999,dims=1)

# storage parameters
InvestmentCostStorage = readin(df_investmentcoststorage,dims=1)
E2PRatio = readin(df_e2pratio,dims=1)
StorageChargeEfficiency = readin(df_storagechargeefficiency,dims=2)
StorageDisChargeEfficiency = readin(df_storagedischargeefficiency,dims=2)
MaxStorageCapacity = readin(df_maxstoragecapacity,default=999,dims=1)
StorageLosses = readin(df_storagelosses, dims=1)

# our emission limit
EmissionLimit = 6000

# instantiate a model with an optimizer

ESM = Model(Clp.Optimizer)

# this creates our variables
@variable(ESM,TotalCost[technologies]>=0)
@variable(ESM,Production[hour,technologies, fuels] >= 0)
@variable(ESM,Capacity[technologies] >=0)
@variable(ESM,Use[hour,technologies, fuels] >=0)
@variable(ESM,Emissions[technologies] >=0)
@variable(ESM,Curtailment[hour,fuels] >=0)

### add variables
@variable(ESM,StorageEnergyCapacity[s=storages,f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageCharge[s=storages, hour, f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageDischarge[s=storages, hour, f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,StorageLevel[s=storages, hour, f=fuels; StorageDisChargeEfficiency[s,f]>0]>=0)
@variable(ESM,TotalStorageCost[storages] >= 0)

# constraints
@constraint(ESM, DemandAdequacy[h in hour,f in fuels], sum(Production[h,t,f] for t in technologies) + sum(StorageDischarge[s,h,f] for s in storages if StorageDisChargeEfficiency[s,f]>0) == Demand[f]*DemandProfile[f,h] + sum(Use[h,t,f] for t in technologies)+Curtailment[h,f] + sum(StorageCharge[s,h,f] for s in storages if StorageChargeEfficiency[s,f] > 0))
@constraint(ESM, ProductionCost[t in technologies], sum(Production[h,t,f] for f in fuels, h in hour)*VariableCost[t]+Capacity[t]*InvestmentCost[t] == TotalCost[t])
@constraint(ESM, ProductionFunction_disp[h in hour, t in technologies, f in fuels; TagDispatchableTechnology[t]>0], OutputRatio[t,f]*Capacity[t]*CapacityFactor[t,h] >= Production[h,t,f])

@constraint(ESM, ProductionFunction_res[h in hour, t in technologies, f in fuels; TagDispatchableTechnology[t]==0], OutputRatio[t,f]*Capacity[t]*CapacityFactor[t,h] == Production[h,t,f])

@constraint(ESM, UseFunction[h in hour,t in technologies, f in fuels], InputRatio[t,f]*sum(Production[h,t,ff]/OutputRatio[t,ff] for ff in fuels if OutputRatio[t,ff]>0) == Use[h,t,f])
@constraint(ESM, TechnologyEmissions[t in technologies], sum(Production[h,t,f] for f in fuels, h in hour)*EmissionRatio[t] == Emissions[t])
@constraint(ESM, TotalEmissionsFunction, sum(Emissions[t] for t in technologies) <= EmissionLimit)
@constraint(ESM, MaxCapacityFunction[t in technologies], Capacity[t] <= MaxCapacity[t])

### Add storage constraints
@constraint(ESM, StorageChargeFunction[s in storages, h in hour, f in fuels; StorageDisChargeEfficiency[s,f]>0], StorageCharge[s,h,f] <= StorageEnergyCapacity[s,f]/E2PRatio[s])
@constraint(ESM, StorageDischargeFunction[s in storages, h in hour, f in fuels; StorageDisChargeEfficiency[s,f]>0], StorageDischarge[s,h,f] <= StorageEnergyCapacity[s,f]/E2PRatio[s])
@constraint(ESM, StorageLevelFunction[s in storages, h in hour, f in fuels; h>1 && StorageDisChargeEfficiency[s,f]>0], StorageLevel[s,h,f] == StorageLevel[s,h-1,f] * StorageLosses[s] + StorageCharge[s,h,f]*StorageChargeEfficiency[s,f] - StorageDischarge[s,h,f]/StorageDisChargeEfficiency[s,f])
@constraint(ESM, StorageLevelStartFunction[s in storages, h in hour, f in fuels; h==1 && StorageDisChargeEfficiency[s,f]>0], StorageLevel[s,h,f] == 0.5*StorageEnergyCapacity[s,f]+ StorageCharge[s,h,f]*StorageChargeEfficiency[s,f] - StorageDischarge[s,h,f]/StorageDisChargeEfficiency[s,f])
@constraint(ESM, MaxStorageLevelFunction[s in storages, h in hour, f in fuels; StorageDisChargeEfficiency[s,f]>0], StorageLevel[s,h,f] <= StorageEnergyCapacity[s,f])
@constraint(ESM, StorageCostFunction[s in storages], TotalStorageCost[s] == sum(StorageEnergyCapacity[s,f]*InvestmentCostStorage[s] for f in fuels if StorageDisChargeEfficiency[s,f]>0))
@constraint(ESM, StorageAnnualBalanceFunction[s in storages, f in fuels; StorageDisChargeEfficiency[s,f]>0], StorageLevel[s,n_hour,f] == 0.5*StorageEnergyCapacity[s,f])
@constraint(ESM, MaxStorageCapacityFunction[s in storages], MaxStorageCapacity[s] >= sum(StorageEnergyCapacity[s,f] for f in fuels if StorageDisChargeEfficiency[s,f]>0))


# the objective function
@objective(ESM, Min, sum(TotalCost[t] for t in technologies) + sum(TotalStorageCost[s] for s in storages))

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

df_res_production = DataFrame(Containers.rowtable(value,Production; header = [:Hour, :Technology, :Fuel, :value]))
df_res_capacity = DataFrame(Containers.rowtable(value,Capacity; header = [:Technology, :value]))

df_storage_production = DataFrame(Containers.rowtable(value,StorageDischarge; header = [:Technology, :Hour, :Fuel, :value]))
df_storage_charge = DataFrame(Containers.rowtable(value,StorageCharge; header = [:Technology, :Hour, :Fuel, :value]))
df_storage_level = DataFrame(Containers.rowtable(value,StorageLevel; header = [:Technology, :Hour, :Fuel, :value]))

append!(df_res_production, df_storage_production)

transform!(df_res_production, "Technology" => ByRow(x-> colors[x]) => "Color")
transform!(df_res_capacity, "Technology" => ByRow(x-> colors[x]) => "Color")

# and some plots
df_total_production = combine(
groupby(df_res_production, ["Fuel", "Technology"]),
"value" => sum => "value"
)
transform!(df_total_production, "Technology" => ByRow(x-> colors[x]) => "Color")
groupedbar(
df_total_production.Fuel,
df_total_production.value,
group=df_total_production.Technology,
bar_position=:stack,
title="Production by Technology",
linewidth=0,
color=df_total_production.Color,
legend=false
)

bar(
    df_res_capacity.Technology,
    df_res_capacity.value,
    title="Installed Capacity by Technology",
    color=df_res_capacity.Color,
    linewidth=0,
    rotation=90
)

gdf_production_by_fuel = groupby(df_res_production, :Fuel)
sto_charge = Dict((row.Technology, row.Fuel, row.Hour) => row.value for row in eachrow(df_storage_charge))

n_fuels = length(gdf_production_by_fuel)
plts = map(enumerate(pairs(gdf_production_by_fuel))) do (i,(k,v))
    p = groupedbar(
        v.Hour,
        v.value,
        group=v.Technology,
        bar_position=:stack,
        title="$(k[1])",
        linewidth=0,
        color=v.Color,
        legend=i == n_fuels ? (0.15,-0.5) : false,
        bottom_margin=i == n_fuels ? 20mm : 2mm,
        legend_column=5
    )

    d = [Demand[k[1]]*DemandProfile[k[1],h] for h in hour]
    u = sum(value.(Use)[:,t, k[1]] for t in technologies)
    c = [sum(get(sto_charge, (s, k[1], h), 0) for s in storages) for h in hour]
    du = d .+ u.data
    dus = d .+ u.data .+ c
    plot!(p, hour, d, color=:black, linewidth=2, label="Demand")
    plot!(p, hour, du, color=:black, linestyle=:dash, linewidth=2, label="Demand + Use")
    plot!(p, hour, dus, color=:black, linestyle=:dot, linewidth=2, label="Demand + Use + Storage")

    return p
end

plot(plts..., layout=(n_fuels,1), size=(1200,1200))


sto_prod = Dict((row.Technology, row.Fuel, row.Hour) => row.value for row in eachrow(df_storage_production))
sto_lvl = Dict((row.Technology, row.Fuel, row.Hour) => row.value for row in eachrow(df_storage_level))

plt_storage_lvl = map(storages) do s
    
   val = [[get(sto_prod, (s, f, h), 0) for h in hour ] for f in fuels if StorageDisChargeEfficiency[s,f]>0]

   color = colors[s]
    p = bar(
        hour,
        val,
        title="$(s)",
        label="Production",
        color=:green,
        ylabel="GW",
        linewidth=0
    )

    val_charg = [[get(sto_charge, (s, f, h), 0) for h in hour ] for f in fuels if StorageDisChargeEfficiency[s,f]>0]

    bar!(p,
        hour,
        -val_charg,
        title="$(s)",
        label="Charge",
        color=:red,
        linewidth=0,
        legend=:bottomleft
    )


    val_lvl = [[get(sto_lvl, (s, f, h), 0) for h in hour ] for f in fuels if StorageDisChargeEfficiency[s,f]>0]


    p2 = twinx(p)

    plot!(
        p2,
        hour,
        val_lvl,
        label = "Level",
        color=:black,
        linewidth=3,
        ylabel="GWh",
        legend=:bottomright
    )

    

    return p
end
plot(plt_storage_lvl...)