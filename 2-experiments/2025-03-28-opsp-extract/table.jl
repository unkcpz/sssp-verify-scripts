# %%

using JSON3
using DataFrames

# %%
json_file = "/home/jyu/project/sssp-project/sssp-verify-scripts/2-experiments/2024-07-plot-element-summary/plot-pp-table-summary/table-report.json"
json_data = JSON3.read(open(json_file))

df = DataFrame()

# Collect all possible column names
column_names = Set{Symbol}()
for (_, value) in json_data
    union!(column_names, keys(value))
    if haskey(value, :nu_confs)
        union!(column_names, keys(value[:nu_confs]))
    else
        push!(column_names, :nu_confs)  # Ensure 'nu_confs' exists
    end
end
push!(column_names, :filename)  # Ensure 'filename' exists

# Ensure `df` is initialized with the correct columns
df =
    isempty(df) ? DataFrame(; [name => Union{Missing,Any}[] for name in column_names]...) :
    df

# Process each JSON entry
for (key, value) in json_data
    row = Dict{Symbol,Any}(name => missing for name in column_names)  # Default all fields to `missing`

    # Copy existing values from `value`
    for (k, v) in value
        if k == :nu_confs
            continue  # Handle nu_confs separately
        end
        row[k] = v
    end

    # Handle `nu_confs`
    if haskey(value, :nu_confs)
        for (k, v) in value[:nu_confs]
            row[k] = v
        end
    end

    row[:filename] = key  # Add filename column

    push!(df, row; promote = true)  # Ensure insertion works
end


rename!(df, :efficiency => :ecut_eff)
rename!(df, :precision => :ecut_prec)

# %%
function eos_score(nu::Union{Missing,Nothing,Real})
    if ismissing(nu) || nu === nothing
        return missing
    end
    nu = nu - 0.1
    return max(nu, 0.0)  # Ensures non-negative output
end

function ecut_score(ecut::Union{Missing,Nothing,Real}, w_ecut::Real)
    if ismissing(ecut) || ecut === nothing
        return missing
    end
    ecut = ecut - 30
    return max(ecut * w_ecut, 0.0)  # Ensures non-negative output
end


confs = [:SC, :BCC, :FCC, :DC, :XO, :XO2, :XO3, :X2O, :X2O3, :X2O5];

# Filter by element

function filter_df(df, element; criteria = :ecut_eff, w_ecut = 1 / 100, verbose = false)
    if element !== Nothing
        edf = df[df.element .== element, :]
    else
        edf = df[:, :]
    end

    # Compute eos_score for all confs while preserving missing values
    transform!(edf, confs .=> ByRow(x -> something(eos_score(x), missing)))

    # Compute eos_score function, handling missing values
    new_confs = Symbol.(string.(confs) .* "_function")
    transform!(
        edf,
        new_confs =>
            ByRow(
                (row...) -> begin
                    row_vals = collect(skipmissing(row))  # Remove missing values
                    if isempty(row_vals)
                        missing
                    else
                        (sum(row_vals) - maximum(row_vals)) / (length(row_vals) - 1)
                    end
                end,
            ) => :eos_score,
    )

    # Compute final score with cutoff consideration
    transform!(
        edf,
        [:eos_score, criteria] =>
            ByRow((x1, x2) -> something(x1, 0) + something(ecut_score(x2, w_ecut), 0)) =>
                :score,
    )

    out_cols = [:element, :abbr_name, :z_valence, :ecut_eff, :ecut_prec, :eos_score, :score]

    if verbose
        append!(out_cols, confs)
    end

    select!(edf, out_cols)

    # Sort the DataFrame while keeping missing values visible
    df_sorted = sort(edf, [:score], rev = false, by = x -> something(x, Inf))

    return df_sorted
end


# %%
# eleemnt to check:
# W, As, Au ...
# In order:
# H, He, Li, Be, ...
using PeriodicTable

df_final = filter_df(df, Nothing; criteria = :ecut_eff, w_ecut = 1 / 100, verbose = true);
rename!(df_final, :abbr_name => :name);
df_final[!, :z_valence] = parse.(Int, df_final[!, :z_valence]);  # Convert to integer
df_final.ZA = [
    haskey(PeriodicTable.elements, Symbol(el)) ?
    PeriodicTable.elements[Symbol(el)].number : missing for el in df_final.element
]
df_final = sort(df_final, [:ZA], rev = false);
# d_final = sort(df_final, [:z_valence], rev = false)
# Apply rounding to selected columns
for col in confs
    df_final[!, col] .= round.(df_final[!, col], digits = 2)
end
println(df_final)

# %% element inspect
filtered_df = df_final[df_final.element .== "C", :];
filtered_df
out_cols = [:element, :name, :z_valence, :eos_score, :score]
append!(out_cols, confs)
select!(filtered_df, out_cols)

# %% filter and write csv
#
using CSV, DataFrames

# CSV.write("filtered_pseudopotential_data_N.csv", filtered_df; transform=(col, val) -> something(val, missing))
CSV.write(stdout, filtered_df; transform = (col, val) -> something(val, missing))


# %%
# plot eos
using Plots
using Plots.PlotMeasures
gr()

function birch_murnaghan(V, E0, V0, B0, B1)
    # Return the energy for given volume (V - it can be a vector) according to
    # the Birch Murnaghan function with parameters E0,V0,B0,B01.
    r = (V0 / V)^(2.0 / 3.0)
    return E0 + 9.0 / 16.0 * B0 * V0 * ((r - 1.0)^3 * B1 + (r - 1.0)^2 * (6.0 - 4.0 * r))
end

palette10 = [
    "#1f77b4",  # blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # olive
    "#17becf",   # cyan
]

function plot_eos(eos_data, element; save = false)
    confs = [:SC, :BCC, :FCC, :DC, :XO, :XO2, :XO3, :X2O, :X2O3, :X2O5]
    data_ref = eos_data[:REF]
    pps = Dict(k => v for (k, v) in eos_data if k != :REF)
    p = plot(
        layout = (5, 2),
        size = (1000, 1200),
        plot_title = "EOS summary: $(element)",
        left_margin = 4.2mm,
    )

    # Loop over 10 confs
    # for conf in confs
    for (i_, conf) in enumerate(confs)
        # x, y labels
        plot!(xlabel = "Cell volume per atom [A 3]")
        plot!(ylabel = "Energy per atom [eV]")


        ref = false
        i = 1
        for (pname, data_pp) in pps
            x_pp = data_pp[conf][:volumes]
            xinter = maximum(x_pp) - minimum(x_pp)

            y_pp = data_pp[conf][:energies] .- data_pp[conf][:E0]
            if ref == false
                ref_V0 = data_ref[conf][:V0]
                ref_B0 = data_ref[conf][:B0]
                ref_B1 = data_ref[conf][:B1]
                x_smooth_ref = range(
                    minimum(x_pp) - xinter / 12,
                    stop = maximum(x_pp) + xinter / 12,
                    length = 100,
                )
                y_smooth_ref = @.birch_murnaghan(x_smooth_ref, 0.0, ref_V0, ref_B0, ref_B1)
                plot!(
                    x_smooth_ref,
                    y_smooth_ref,
                    subplot = i_,
                    linestyle = :dash,
                    color = :black,
                    label = "AE $conf",
                )
                ref = true
            end

            V0 = data_pp[conf][:V0]
            B0 = data_pp[conf][:B0]
            B1 = data_pp[conf][:B1]
            nu = data_pp[conf][:nu]
            nu = round(nu; digits = 2)
            x_smooth_pp = range(
                minimum(x_pp) - xinter / 12,
                stop = maximum(x_pp) + xinter / 12,
                length = 100,
            )
            y_smooth_pp = @.birch_murnaghan(x_smooth_pp, 0.0, V0, B0, B1)
            label = join(split(String(pname), '.')[4:end], '.')

            scatter!(x_pp, y_pp, subplot = i_, color = :red, markersize = 2, label = "")
            plot!(
                x_smooth_pp,
                y_smooth_pp,
                subplot = i_,
                label = "$(label) ($(nu))",
                color = palette10[i],
            )
            i = i + 1

        end


    end

    if save
        @info savefig(p, "/tmp/x.pdf")
    else
        p
    end
end

# %%
# EOS data extract
using JSON3
element = "N"
jfh = "/home/jyu/project/sssp-project/sssp-verify-scripts/2-experiments/2025-03-28-opsp-extract/$(element).EOS.json"
eos_data = JSON3.read(open(jfh))
p = plot_eos(eos_data, element)
display(p)

# %%
p = plot_eos(eos_data, "Au", "SG15", "GBRV-1.X")
display(p)

# %%
eos_data[:Au]

# %%
