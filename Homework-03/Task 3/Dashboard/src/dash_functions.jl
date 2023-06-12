
function stack_traces(t)

    isempty(t) && error("traces must be non-empty")
    traces = deepcopy(t)
    traces[1]["fill"] = "tozeroy"
    traces[1]["mode"] ="none"
    traces[1]["fillcolor"] = first(traces[1]["color"])


    if length(traces) > 1
        for (i, tr) in enumerate(traces[2:end])
            tr["fill"] = "tonexty"
            tr["mode"] ="none"
            tr["fillcolor"] = first(tr["color"])

            for j in 1:min(length(traces[i]["y"]), length(tr["y"]))
                tr["y"][j] += traces[i]["y"][j]
            end
        end
    end
    traces
end

function stackedarea(args...; kwargs...)
    traces = scatter(args...; kwargs...)
    return stack_traces(traces)
end  # function stackedarea

struct Result
    Hour::Vector{Int}
    Technology::Vector{String}
    Fuel::Vector{String}
    Region::Vector{String}

    production::DataFrame
    use::DataFrame
    demand::DataFrame
    capacity::DataFrame
    ex_import::DataFrame
    ex_export::DataFrame

    colors::Dict{String, String}

    function Result(path)
        results = Dict(filename(file) => readcsv(file) for file in readdir(path, join=true) if endswith(file, ".csv"))
        colors = JSON3.read(read(joinpath(path, "colors.json"), String), Dict)

        @unpack production, use, demand, capacity, ex_import, ex_export = results

        Hour = unique(vcat(production.Hour, use.Hour, demand.Hour))
        Technology = unique(vcat(production.Technology, use.Technology))
        Fuel = unique(vcat(production.Fuel, use.Fuel, demand.Fuel))
        Region = unique(vcat(production.Region, use.Region, demand.Region))

        new(Hour, Technology, Fuel, Region, production, use, demand, capacity, ex_import, ex_export, colors)
    end
end


add_color!(df, c) = transform!(df, :Technology => ByRow(x -> c[x]) => :Color)

function filterByFuelRegion(df, fuel, region; filter_non_zero_techs=false)
    df = filter(row-> row.Fuel == fuel && row.Region == region, df)
    if filter_non_zero_techs
        techs_in_df_not_zero = Dict(k.Technology => any(x-> x != 0,v.value) for (k,v) in pairs(groupby(df, :Technology)))
        filter!(row -> techs_in_df_not_zero[row.Technology], df)
    end
    sort!(df, :Hour)
end

function filterExchangeFuelRegion(df, fuel, region, which, colors)
    which == "Import" ? col = :To : col = :From
    df = filter(row-> row.Fuel == fuel && row[col] == region, df)
    df = combine(groupby(df, [:Hour, :Fuel, col]), :value => sum => :value)
    sort!(df, :Hour)
    rename!(df, col => :Region)
    df[!, :Technology] .= which
    df[!, :Color] .= colors[which]
    return df
end

function plot_dispatch(r::Result, fuel, region; colors=colors)
        
    selected_prod = filterByFuelRegion(r.production, fuel, region, filter_non_zero_techs=true)
    selected_use = filterByFuelRegion(r.use, fuel, region, filter_non_zero_techs=true)
    selected_demand = filterByFuelRegion(r.demand, fuel, region)
    selected_demand[!, :Technology] .= "Demand"
    selected_demand[!, :Color] .= "black"

    selected_import = filterExchangeFuelRegion(r.ex_import, fuel, region, "Import", r.colors)
    selected_export = filterExchangeFuelRegion(r.ex_export, fuel, region, "Export", r.colors)

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

function plot_total_by_fuel(r::Result, region=""; style="relative")
        
    filt(df) = filter(row -> row.Region == region, df)

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
        prod = r.production
        use = r.use
        demand = r.demand
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


    transform!(agg_prod, :Technology => ByRow(x -> r.colors[x]) => :Color)
    transform!(agg_use, :Technology => ByRow(x -> r.colors[x]) => :Color)
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


function plot_capacity(r::Result, region="")
        
    if region != ""
        cap = filter(row -> row.Region == region, r.capacity)
    else
        cap = combine(groupby(r.capacity, "Technology"), "value" => sum => "value")
    end

    transform!(cap, :Technology => ByRow(x -> r.colors[x]) => :Color)
    traces = bar(cap, x=:Technology, y=:value, marker_color=:Color)
    layout = Layout(;barmode="stack", yaxis=attr(title="GW"))

    return Plot(traces, layout)
end



function dashboard(result::Result; window=true, debug=false)
        
    dropdown_options_region = [Dict("label" => r, "value" => r) for r in result.Region]
    dropdown_options_fuel = [Dict("label" => f, "value" => f) for f in result.Fuel]

    app = dash(external_stylesheets = ["bWLwgP.css"], assets_folder="assets")

    app.layout = html_div() do
        dcc_tabs(id="tabs", value="tab-1-example-graph", children=[
            dcc_tab(label="Dispatch", value="tab-1-example-graph", children=[
                html_div(className = "row") do
                    html_div(
                        className = "one column",
                        style = Dict("align-items" => "center", "justify-content" => "center")
                    ) do
                        dcc_dropdown(
                            id = "dropdown1",
                            options = dropdown_options_fuel,
                            value = "Power",
                            clearable=false
                        ),
                        dcc_dropdown(
                            id = "dropdown2",
                            options = dropdown_options_region,
                            value = "DE",
                            clearable=false
                        )
                    end,
                    html_div(className = "eleven columns") do
                        dcc_graph(
                            id = "dispatch plot 1",
                            figure = plot_dispatch(result, "Power", "DE")
                        )
                    end
                end,
            
                html_div(className = "row") do
                    html_div(
                            className = "one column",
                            style = Dict("align-items" => "center", "justify-content" => "center")
                        ) do
                        dcc_dropdown(
                            id = "dropdown3",
                            options = dropdown_options_fuel,
                            value = "Power",
                            clearable=false
                        ),
                        dcc_dropdown(
                            id = "dropdown4",
                            options = dropdown_options_region,
                            value = "FR",
                            clearable=false
                        )
                    end,
                    html_div(className = "eleven columns") do
                        dcc_graph(
                            id = "dispatch plot 2",
                            figure = plot_dispatch(result, "Power", "FR")
                        )
                    end
                end
            ]),
            dcc_tab(label="Total", value="Total", children=[
                html_div(style = Dict("align-items" => "center", "justify-content" => "center", "display" => "flex")) do
                    dcc_dropdown(
                        id = "dropdown5",
                        options = dropdown_options_region,
                        value = nothing
                    )
                end,

                html_div(className = "row") do
                    html_div(className = "six columns") do
                        html_h2("Total Generation"),
                        dcc_graph(
                            id = "generation plot",
                            figure = plot_total_by_fuel(result, "")
                        )
                    end,
                    html_div(className = "six columns") do
                        html_h2("Installed Capacity"),
                        dcc_graph(
                            id = "capacity plot",
                            figure = plot_capacity(result, "")
                        )
                    end
                end,
            ]),
            dcc_tab(label="Sankey", value="Sankey", children=[
                html_div(style = Dict("align-items" => "center", "justify-content" => "center", "display" => "flex")) do
                    dcc_dropdown(
                        id = "dropdown6",
                        options = dropdown_options_region,
                        value = nothing
                    )
                end,
                html_div(className = "twelve columns") do
                    dcc_graph(
                        id = "sankey plot",
                        figure = plot_sankey(result, "")
                    )
                end
            ]),
        ])
        
    end

    callback!(
        app,
        Output("dispatch plot 1", "figure"),
        Input("dropdown1", "value"),
        Input("dropdown2", "value")) do v1, v2
        return plot_dispatch(result, v1, v2)
    end

    callback!(
        app,
        Output("dispatch plot 2", "figure"),
        Input("dropdown3", "value"),
        Input("dropdown4", "value")) do v1, v2
        return plot_dispatch(result, v1, v2)
    end

    callback!(
        app,
        Output("generation plot", "figure"),
        Input("dropdown5", "value")) do v
        isnothing(v) ? v = "" : v = v
        return plot_total_by_fuel(result, v)
    end

    callback!(
        app,
        Output("capacity plot", "figure"),
        Input("dropdown5", "value")) do v
        isnothing(v) ? v = "" : v = v
        return plot_capacity(result, v)
    end

    callback!(
        app,
        Output("sankey plot", "figure"),
        Input("dropdown6", "value")) do v
        isnothing(v) ? v = "" : v = v
        return plot_sankey(result, v)
    end

    if window
        w = Window()
        Threads.@spawn run_server(app, debug=true)
        loadurl(w, "http://127.0.0.1:8050")
    else
        run_server(app, debug=debug)
    end
end
