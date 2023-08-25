

add_color!(df, c) = transform!(df, :Technology => ByRow(x -> get(c,x,default_color)) => :Color)

function filterByYearFuelRegion(df, year, fuel, region; filter_non_zero_techs=false)
    df = @chain df begin
        @rsubset :Year == year
        @rsubset :Fuel == fuel
        @rsubset :Region == region
    end

    if filter_non_zero_techs
        techs_in_df_not_zero = Dict(k.Technology => any(x-> x != 0,v.value) for (k,v) in pairs(groupby(df, :Technology)))
        filter!(row -> techs_in_df_not_zero[row.Technology], df)
    end
    sort!(df, :Hour)
end

function filterExchangeYearFuelRegion(df, year,fuel, region, which, colors)
    which == "Import" ? col = :To : col = :From
    
    df = @chain df begin
        @rsubset :Year == year
        @rsubset :Fuel == fuel
        filter!(row-> row[col] == region, _)
        groupby([:Year, :Hour, :Fuel, col])
        combine(:value => sum => :value)
        sort!(:Hour)
        rename!(col => :Region)
        @rtransform!  :Technology = which :Color = colors[which]
    end

    return df
end

function plot_dispatch(r::Result, year, fuel, region; colors=colors)
        
    selected_prod = filterByYearFuelRegion(r.production, year, fuel, region, filter_non_zero_techs=true)
    selected_use = filterByYearFuelRegion(r.use, year, fuel, region, filter_non_zero_techs=true)
    selected_demand = filterByYearFuelRegion(r.demand, year, fuel, region)
    selected_demand[!, :Technology] .= "Demand"
    selected_demand[!, :Color] .= "black"

    selected_import = filterExchangeYearFuelRegion(r.ex_import, year, fuel, region, "Import", r.colors)
    selected_export = filterExchangeYearFuelRegion(r.ex_export, year, fuel, region, "Export", r.colors)

    add_color!(selected_prod, r.colors)
    add_color!(selected_use, r.colors)

    selected_use = vcat(selected_demand, selected_use, selected_export)
    selected_prod = vcat(selected_prod, selected_import)
    transform!(selected_use, :value => ByRow(x -> -x) => :value)

    traces_prod = stackedarea(selected_prod, x=:Hour, y=:value, group=:Technology, color=:Color)
    traces_use = stackedarea(selected_use, x=:Hour, y=:value, group=:Technology, color=:Color)
    traces = vcat(traces_prod, traces_use)
    return Plot(traces)
end

function plot_total_by_fuel(r::Result, y, region=""; style="relative")
        
    filt(df) = filter(row -> row.Region == region && row.Year == y, df)
    filt_year(df) = filter(row -> row.Year == y, df)

    if region != ""
       prod =  filt(r.production)
       use =  filt(r.use)
       demand =  filt(r.demand)
       ex_im = filter(row -> row.To == region, r.ex_import)
       ex_im = combine(groupby(ex_im, "Fuel"), "value" => sum => "value")
       filter!(row-> row.value != 0, ex_im)
       ex_ex = filter(row -> row.From == region, r.ex_export)
       ex_ex = combine(groupby(ex_ex, "Fuel"), "value" => sum => "value")
       filter!(row-> row.value != 0, ex_ex)
       transform!(ex_ex, :value => ByRow(x -> -x) => :value)
    else
        prod = filt_year(r.production)
        use = filt_year(r.use)
        demand = filt_year(r.demand)
    end

    agg_prod = combine(
        groupby(prod, ["Technology", "Fuel"]),
        "value" => sum=> "value"
    )
    filter!(row -> row.value != 0, agg_prod)

    agg_use = combine(
        groupby(use, ["Technology", "Fuel"]),
        "value" => sum=> "value"
    )
    filter!(row -> row.value != 0, agg_use)

    agg_demand = combine(
        groupby(demand, "Fuel"),
        "value" => sum=> "value"
    )


    transform!(agg_prod, :Technology => ByRow(x -> get(r.colors,x,default_color)) => :Color)
    transform!(agg_use, :Technology => ByRow(x -> get(r.colors,x,default_color)) => :Color)
    if style == "relative"
        transform!(agg_use, :value => ByRow(x -> -x) => :value)
        transform!(agg_demand, :value => ByRow(x -> -x) => :value)
    end

    traces_prod = bar(agg_prod, x=:Fuel, y=:value, group=:Technology, marker_color=:Color)
    traces_use = bar(agg_use, x=:Fuel, y=:value, group=:Technology, marker_color=:Color)
    traces_demand = bar(agg_demand, x=:Fuel, y=:value, marker_color="black", name="Demand")
    layout = Layout(;barmode=style, yaxis=attr(title="GWh"))

    if region != ""
        traces_im = bar(ex_im, x=:Fuel, y=:value, marker_color=r.colors["Import"], name="Import")
        traces_ex = bar(ex_ex, x=:Fuel, y=:value, marker_color=r.colors["Export"], name="Export")
        traces = vcat(traces_prod, traces_im, traces_use, traces_demand, traces_ex)
    else
        traces = vcat(traces_prod, traces_use, traces_demand)
    end

    return Plot(traces, layout)
end


function plot_capacity(r::Result, region=""; variable=:capacity, showlegend=true)
        
    df = getfield(r, variable)
    if region != ""
        cap = filter(row -> row.Region == region, df)
    else
        cap = combine(groupby(df, ["Year", "Technology"]), "value" => sum => "value")
    end

    transform!(cap, :Technology => ByRow(x -> get(r.colors,x,default_color)) => :Color)
    traces = bar(cap, x=:Year, y=:value, group=:Technology, marker_color=:Color)
    layout = Layout(;barmode="stack", yaxis=attr(title="GW"), showlegend=showlegend)

    return Plot(traces, layout)
end

function plot_emission(r::Result, g::String)
    g = Symbol(g)
    df = @chain r.emission begin
        groupby([:Year, g])
        combine("value" => sum => "value")
    end

    if g == :Technology
        transform!(df, :Technology => ByRow(x -> get(r.colors,x,default_color)) => :Color)
        traces = bar(df, x=:Year, y=:value, group=:Technology, marker_color=:Color)
    elseif g == :Region
        traces = bar(df, x=:Year, y=:value, group=:Region, marker_color=:Color)
    end
    layout = Layout(;barmode="relative", yaxis=attr(title="Co2 emissions (kt)"))

    return Plot(traces, layout)
end

function plot_sequester(r::Result, g::String)
    g = Symbol(g)
    df = @chain r.sequester begin
        groupby([:Year, g])
        combine("value" => sum => "value")
    end

    if g == :Technology
        transform!(df, :Technology => ByRow(x -> get(r.colors,x,default_color)) => :Color)
        traces = bar(df, x=:Year, y=:value, group=:Technology, marker_color=:Color)
    elseif g == :Region
        traces = bar(df, x=:Year, y=:value, group=:Region, marker_color=:Color)
    end
    layout = Layout(;barmode="relative", yaxis=attr(title="Co2 sequester (kt)"))

    return Plot(traces, layout)
end