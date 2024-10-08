# %%

using JSON3
using DataFrames

# %%
json_file = "/home/jyu/project/sssp-project/sssp-verify-scripts/2-experiments/2024-07-plot-element-summary/plot-pp-table-summary/table-report.json"
json_data = JSON3.read(open(json_file))

df = DataFrame()

for (key, value) in json_data
    # Create a copy of the data without the 'nu_confs' field
    if !haskey(value, :nu_confs)
        continue
    end

    nu_confs = value[:nu_confs]
    row = copy(value)
    delete!(row, :nu_confs)  # Remove 'nu_confs' field

    for (k, v) in nu_confs
        row[k] = v
    end

    for key in keys(row)
        if isnothing(row[key])
            row[key] = missing
        end
    end

    # Add the key (filename) as a separate column
    row[:filename] = key

    if nrow(df) == 0
        df = DataFrame(row)
    else
        push!(df, row; promote = true)
    end
end

rename!(df, :efficiency => :ecut_eff)
rename!(df, :precision => :ecut_prec)

println(df)


# %%
function eos_score(nu)
    nu = nu - 0.1
    if nu < 0
        return 0.0
    end

    nu
end

function ecut_score(ecut, w_ecut)
    if ismissing(ecut)
        return missing
    end

    ecut = ecut - 30
    if ecut < 0
        return 0.0
    end

    ecut * w_ecut
end

confs = [:SC, :BCC, :FCC, :DC, :XO, :XO2, :XO3, :X2O, :X2O3, :X2O5];

# Filter by element
function filter_df(df, element; criteria = :ecut_eff, w_ecut = 1 / 100, verbose = false)
    edf = df[df.element .== element, :]
    # compute eos score for all confs
    transform!(edf, confs .=> x -> @.eos_score(x))

    # compute eos score
    new_confs = Symbol.(string.(confs) .* "_function")
    transform!(
        edf,
        new_confs =>
            ((row...) -> (sum(row) - maximum(row)) / (length(row) - 1)) => :eos_score,
    )

    # compute final score with taking cutoff into account
    transform!(
        edf,
        [:eos_score, criteria] => ((x1, x2) -> x1 .+ @.ecut_score(x2, w_ecut)) => :score,
    )

    out_cols =
        [:element, :abbr_name, :z_valence, criteria, :eos_score, :score]

    if verbose
        append!(out_cols, confs)
    end

    select!(edf, out_cols)

    # Sort the non-missing DataFrame
    df_sorted = sort(edf, :score)
    df_sorted
end

# %%
# eleemnt to check: 
# W, As, Au ...
# In order:
# H, He, Li, Be, ...
df_final = filter_df(df, "As"; criteria = :ecut_eff, w_ecut = 1 / 100, verbose = false);
println(df_final)

# %%
# EOS data extract
jfh = "/home/jyu/project/sssp-project/sssp-verify-scripts/2-experiments/2024-07-plot-element-summary/plot-pp-table-summary/eos.json"
eos_data = JSON3.read(open(jfh))

# %%
p = plot_eos(eos_data, "Fe", "PSL-US-v1-low", "GBRV-1.X")
display(p)

# %%
# plot eos
using Plots
gr()

function birch_murnaghan(V, E0, V0, B0, B1)
    # Return the energy for given volume (V - it can be a vector) according to
    # the Birch Murnaghan function with parameters E0,V0,B0,B01.
    r = (V0 / V)^(2.0 / 3.0)
    return E0 + 9.0 / 16.0 * B0 * V0 * ((r - 1.0)^3 * B1 + (r - 1.0)^2 * (6.0 - 4.0 * r))
end

function plot_eos(eos_data, element, pp1, pp2)
    data_ref = eos_data[Symbol(element)][:REF]
    data_pp1 = eos_data[Symbol(element)][Symbol(pp1)]
    data_pp2 = eos_data[Symbol(element)][Symbol(pp2)]

    p = plot(layout = (5, 2), size = (1000, 1200))

    # Loop over 10 confs
    # for conf in confs
    for (i_, conf) in enumerate(confs)
        x_pp1 = data_pp1[conf][:volumes]
        V0 = data_pp1[conf][:V0]
        B0 = data_pp1[conf][:B0]
        B1 = data_pp1[conf][:B1]
        xinter = maximum(x_pp1) - minimum(x_pp1)
        x_smooth_pp1 = range(
            minimum(x_pp1) - xinter / 12,
            stop = maximum(x_pp1) + xinter / 12,
            length = 100,
        )
        y_smooth_pp1 = @.birch_murnaghan(x_smooth_pp1, 0.0, V0, B0, B1)

        y_pp1 = data_pp1[conf][:energies] .- data_pp1[conf][:E0]
        scatter!(x_pp1, y_pp1, subplot = i_, color = :red, markersize = 2, label = "")
        plot!(x_smooth_pp1, y_smooth_pp1, subplot = i_, label = pp1, color = :blue)

        x_pp2 = data_pp2[conf][:volumes]
        V0 = data_pp2[conf][:V0]
        B0 = data_pp2[conf][:B0]
        B1 = data_pp2[conf][:B1]
        xinter = maximum(x_pp2) - minimum(x_pp2)
        x_smooth_pp2 = range(
            minimum(x_pp2) - xinter / 12,
            stop = maximum(x_pp2) + xinter / 12,
            length = 100,
        )
        y_smooth_pp2 = @.birch_murnaghan(x_smooth_pp2, 0.0, V0, B0, B1)

        y_pp2 = data_pp2[conf][:energies] .- data_pp2[conf][:E0]
        scatter!(x_pp2, y_pp2, subplot = i_, color = :red, markersize = 2, label = "")
        plot!(x_smooth_pp2, y_smooth_pp2, subplot = i_, label = pp2, color = :red)

        ref_V0 = data_ref[conf][:V0]
        ref_B0 = data_ref[conf][:B0]
        ref_B1 = data_ref[conf][:B1]
        x_smooth_ref = range(
            minimum(x_pp2) - xinter / 12,
            stop = maximum(x_pp2) + xinter / 12,
            length = 100,
        )
        y_smooth_ref = @.birch_murnaghan(x_smooth_ref, 0.0, ref_V0, ref_B0, ref_B1)
        plot!(x_smooth_ref, y_smooth_ref, subplot = i_, linestyle = :dash, color = :black, label = "AE $conf")

        # x, y labels
        plot!(xlabel="Cell volume per atom [A 3]")
        plot!(ylabel="Energy per atom [eV]")
    end

    p
    # savefig(p, "/tmp/x.pdf")
end

# %%
p = plot_eos(eos_data, "Au", "SG15", "GBRV-1.X")
display(p)

# %%
eos_data[:Au]

# %%
