const default_color = "#000000"


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
    Year::Vector{Int}
    Hour::Vector{Int}
    Technology::Vector{String}
    Fuel::Vector{String}
    Region::Vector{String}

    production::DataFrame
    use::DataFrame
    demand::DataFrame
    capacity::DataFrame
    newcapacity::DataFrame
    ex_import::DataFrame
    ex_export::DataFrame
    emission::DataFrame

    colors::Dict{String, String}

    function Result(path)
        results = Dict(filename(file) => readcsv(file) for file in readdir(path, join=true) if endswith(file, ".csv"))
        colors = JSON3.read(read(joinpath(path, "colors.json"), String), Dict)

        @unpack production, use, demand, capacity, newcapacity, ex_import, ex_export, emission = results

        Year = unique(demand.Year)
        Hour = unique(vcat(production.Hour, use.Hour, demand.Hour))
        Technology = unique(vcat(production.Technology, use.Technology))
        Fuel = unique(vcat(production.Fuel, use.Fuel, demand.Fuel))
        Region = unique(vcat(production.Region, use.Region, demand.Region))

        new(Year, Hour, Technology, Fuel, Region, production, use, demand, capacity, newcapacity, ex_import, ex_export, emission, colors)
    end
end

function dashboard(result::Result; window=true, debug=false)
        
    dropdown_options_region = [Dict("label" => r, "value" => r) for r in result.Region]
    dropdown_options_fuel = [Dict("label" => f, "value" => f) for f in result.Fuel]
    dropdown_options_year = [Dict("label" => f, "value" => f) for f in result.Year]

    asset_path = joinpath(@__DIR__, "assets")
    
    app = dash(assets_folder=asset_path)

    app.layout = html_div() do
        dcc_tabs(id="tabs", value="tab-1-example-graph", children=[
            dcc_tab(label="Dispatch", value="tab-1-example-graph", children=[
                html_div(
                    className = "two columns",
                    style = Dict("margin-top" => 60, "align-items" => "center", "justify-content" => "center")
                ) do
                    html_h5("Filter options"),
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
                    ),
                    dcc_dropdown(
                        id = "dropdown_dispatch_year",
                        options = dropdown_options_year,
                        value = 2020,
                        clearable=false
                    )
                end,
                html_div(className = "ten columns") do
                    dcc_graph(
                        id = "dispatch plot 1",
                        figure = plot_dispatch(result, 2020, "Power", "DE")
                    )
                end,
                html_div(className = "row") do
                    html_div(
                            className = "two columns",
                            style=Dict("margin-top" => 60)
                        ) do
                            html_h5("Filter options"),
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
                            ),
                            dcc_dropdown(
                                id = "dropdown_dispatch2_year",
                                options = dropdown_options_year,
                                value = 2020,
                                clearable=false
                            )
                    end,
                    html_div(className = "ten columns") do
                        dcc_graph(
                            id = "dispatch plot 2",
                            figure = plot_dispatch(result, 2020, "Power", "FR")
                        )
                    end
                end
            ]),
            dcc_tab(label="Total", value="Total", children=[
                html_h2(
                    "Total Generation",
                    style=Dict("display"=> "flex","align-items" => "center", "justify-content" => "center")
                ),
                html_div(children=[
                    html_div(
                        className="two columns",
                        style = style = Dict("margin-top" => 60, "align-items" => "center", "justify-content" => "center"),
                        children=[

                        html_h5("Filter options"),
                        dcc_dropdown(
                            id = "dropdown5",
                            options = dropdown_options_region,
                            value = nothing
                        ),
                        dcc_dropdown(
                            id = "dropdown_total_year",
                            options = dropdown_options_year,
                            value = 2020,
                            clearable=false
                        )
                    ]),
                    html_div(className="ten columns", children=[
                        
                        dcc_graph(
                            id = "generation plot",
                            figure = plot_total_by_fuel(result, "")
                        )
                    ])
                ]),
                html_div(
                    children=[
                        html_div(className="one-half column", children=[
                            html_h2("Installed Capacity", style=Dict("display"=> "flex","align-items" => "center", "justify-content" => "center")),
                            dcc_graph(
                                id = "capacity plot",
                                figure = plot_capacity(result, "", showlegend=false)
                            )
                        ]),
                        html_div(className="one-half column", children=[
                            html_h2("Newly installed Capacity", style=Dict("display"=> "flex","align-items" => "center", "justify-content" => "center")),
                            dcc_graph(
                                id = "newcapacity plot",
                                figure = plot_capacity(result, ""; variable=:newcapacity)
                            )
                        ])
                    ]
                )
            ]),
            dcc_tab(label="Emission", value="tab Emission", children=[
                html_div(children=[
                    html_div(className="six columns", children=[
                        dcc_graph(
                            id = "emission plot by technology",
                            figure = plot_emission(result, "Technology")
                        )
                    ]),
                    html_div(className="six columns", children=[
                        dcc_graph(
                            id = "emission plot by region",
                            figure = plot_emission(result, "Region")
                        )
                    ])
                ])
            ]),
            dcc_tab(label="Sankey", value="Sankey", children=[
                html_h2("Sankey Diagramm"),
                html_div(children=[
                    html_div(className="two columns", children=[
                        html_h5("Filter options"),
                        dcc_dropdown(
                            id = "dropdown6",
                            options = dropdown_options_region,
                            value = nothing
                        ),
                        dcc_dropdown(
                            id = "dropdown_sankey_year",
                            options = dropdown_options_year,
                            value = 2020,
                            clearable=false
                        )
                    ]),
                    html_div(className = "ten columns", children=[
                        dcc_graph(
                            id = "sankey plot",
                            figure = plot_sankey(result, 2020, "")
                        )
                    ])
                ])
            ])
        ])
    end

    callback!(
        app,
        Output("dispatch plot 1", "figure"),
        Input("dropdown1", "value"),
        Input("dropdown2", "value"),
        Input("dropdown_dispatch_year", "value")) do v1, v2, v3
        return plot_dispatch(result, v3, v1, v2)
    end

    callback!(
        app,
        Output("dispatch plot 2", "figure"),
        Input("dropdown3", "value"),
        Input("dropdown4", "value"),
        Input("dropdown_dispatch2_year", "value")) do v1, v2, v3
        return plot_dispatch(result, v3, v1, v2)
    end

    callback!(
        app,
        Output("generation plot", "figure"),
        Input("dropdown5", "value"),
        Input("dropdown_total_year", "value")) do v1, v2
        isnothing(v1) ? v1 = "" : v1 = v1
        return plot_total_by_fuel(result, v2, v1)
    end

    callback!(
        app,
        Output("capacity plot", "figure"),
        Input("dropdown5", "value")) do v
        isnothing(v) ? v = "" : v = v
        return plot_capacity(result, v, showlegend=false)
    end
    
    callback!(
        app,
        Output("newcapacity plot", "figure"),
        Input("dropdown5", "value")) do v
        isnothing(v) ? v = "" : v = v
        return plot_capacity(result, v; variable=:newcapacity)
    end
    
    callback!(
        app,
        Output("sankey plot", "figure"),
        Input("dropdown6", "value"),
        Input("dropdown_sankey_year", "value")) do v1, v2
        isnothing(v1) ? v1 = "" : v1 = v1
        return plot_sankey(result, v2, v1)
    end

    if window
        w = Window(Dict("width"=>1000, "height"=>1200))
        Threads.@spawn run_server(app, debug=true)
        loadurl(w, "http://127.0.0.1:8050")
    else
        run_server(app, "127.0.0.1:8050", debug=debug)
    end
end
