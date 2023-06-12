
function plot_sankey(r::Result, region="")

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

    source = Int[]
    target = Int[]
    value = Float64[]
    link_label = String[]

    node_dict = Dict(f => i for (i,f) in enumerate(r.Fuel))
    region != "" && (node_dict["Import"] = length(node_dict) + 1)

    for tech in r.Technology
        node_dict[tech] = length(node_dict) + 1
    end
 
    node_dict["Demand"] = length(node_dict) + 1
    region != "" && (node_dict["Export"] = length(node_dict) + 1) 

    for f in r.Fuel
        prod = filter(row -> row.Fuel == f, agg_prod)
        use = filter(row -> row.Fuel == f, agg_use)

        for row in eachrow(use)
            push!(target, node_dict[row.Technology])
            push!(link_label, f)
            push!(source, node_dict[f])
            push!(value, row.value)
        end

        for row in eachrow(prod)
            push!(source, node_dict[row.Technology])
            push!(link_label, f)
            push!(target, node_dict[f])
            push!(value, row.value)
        end
    end

    if region != ""
        for row in eachrow(ex_im)
            push!(source, node_dict["Import"])
            push!(link_label, row.Fuel)
            push!(target, node_dict[row.Fuel])
            push!(value, row.value)
        end
    end

    for row in eachrow(agg_demand)
        push!(source, node_dict[row.Fuel])
        push!(link_label, row.Fuel)
        push!(target, node_dict["Demand"])
        push!(value, row.value)
    end

    if region != ""
        for row in eachrow(ex_ex)
            push!(target, node_dict["Export"])
            push!(link_label, row.Fuel)
            push!(source, node_dict[row.Fuel])
            push!(value, row.value)
        end
    end

        # Define node labels
    labels = fill("", length(node_dict)+1)
    for (k,v) in pairs(node_dict)
        labels[v+1] = k
    end

    # Define link data
    link_data = attr(
        source = source,
        target = target,
        value = value,
        label = link_label,
        color = [r.colors[l] for l in link_label]
    )

    # Define node data
    node_data = attr(
        pad = 30,
        thickness = 20,
        line = attr(
            color = "black",
            width = 0.5
        ),
        label = labels,
        color = [get(r.colors, l, "") for l in labels]
    )

    # Define Sankey trace
    sankey_trace = sankey(
        domain = attr(
            x = [0, 1],
            y = [0, 1]
        ),
        node = node_data,
        link = link_data
    )

    # Create layout
    layout = Layout(
        title_text = "Sankey diagram",
    )

    # Create plot
    plot(sankey_trace, layout)
end

