using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using Distributions, Random
using Plots, LaTeXStrings

# computes gap, q⁺, and q⁻ for 'optimal' and FIFO-includion block
function simulate_block(val_dist, tx_size_dist, N, k, T)
    gaps = zeros(T)
    q_up = zeros(T)
    q_low = zeros(T)
    for t in 1:T
        tx_vals = abs.(rand(val_dist, N))
        tx_sizes = rand(tx_size_dist, N)

        b = k*mean(tx_size_dist)

        tx_eff = tx_vals ./ tx_sizes
        sorted_inds = sortperm(tx_eff, rev=true)

        block_fifo_end = findlast(cumsum(tx_sizes) .<= b)
        block_fifo = collect(1:block_fifo_end)

        block_greedy_end = findlast(cumsum(tx_sizes[sorted_inds]) .<= b)
        block_greedy = sorted_inds[1:block_greedy_end]

        gaps[t] = sum(tx_vals[block_greedy]) / sum(tx_vals[block_fifo])
        q_up[t] = mean(tx_vals[block_greedy])
        q_low[t] = (sum(tx_vals) - sum(tx_vals[block_greedy])) /(N - k)
    end

    return gaps, q_up, q_low
end

Random.seed!(1234)
fig_path = joinpath(@__DIR__, "figs")

# Parameters
B⁺ = 3
B⁻ = 1
η = B⁺ / B⁻
N = 1_000
ks = 10:10:1000
T = 100

tx_size_dist = Uniform(B⁻, B⁺)
dists = [Exponential(2.5), LogNormal(1,1), Rayleigh(1)]
dist_labels = ["Exponential" "LogNormal" "Rayleigh"]

## Light distributions
x = 0:0.01:10
dist_plt = plot()
for (i, val_dist) in enumerate(dists)
    plot!(dist_plt, x, pdf.(val_dist, x), label="$(dist_labels[i])", lw=2)
end
display(dist_plt)
savefig(dist_plt, joinpath(fig_path, "dists.pdf"))


# Run trials
gaps = zeros(T, length(ks), length(dists))
q_ups = zeros(T, length(ks), length(dists))
q_lows = zeros(T, length(ks), length(dists))
for (i, val_dist) in enumerate(dists), (j, k) in enumerate(ks)
    ret = simulate_block(val_dist, tx_size_dist, N, k, T)
    gaps[:, j, i] .= ret[1]
    q_ups[:, j, i] .= ret[2]
    q_lows[:, j, i] .= ret[3]
end

μs = mean(gaps, dims=1)[1,:,:]
meds = median(gaps, dims=1)[1,:,:]
σs = std(gaps, dims=1)[1,:,:]
q_ups = mean(q_ups, dims=1)[1,:,:]
q_lows = mean(q_lows, dims=1)[1,:,:]

# Plot results
plts = []
colors = [:coral :indigo :firebrick :mediumblue]
for i in 1:length(dists)
    q_up = q_ups[:, i]
    q_low = q_lows[:, i]
    b1 = @. 1/η * q_up / ( (q_up - q_low) * ks / N + q_low )

    plt = plot(
        ks,
        μs[:,i],
        lw=2,
        ribbon=σs[:,i],  
        fillalpha=0.5,
        xaxis=:log,
        ylabel=L"$p^\star/p^\mathrm{fifo}$",
        xlabel="k = average block size (in # txs)",
        legend=:false,
        dpi=300,
        margin=3Plots.PlotMeasures.mm,
        color = colors[i]
    ) 
    plot!(plt, ks, b1, lw=2, label="Gap bound", ls=:dash, color=colors[i])
    
    push!(plts, plt)
end

[display(plt) for plt in plts]
for (i, plt) in enumerate(plts)
    savefig(plt, joinpath(fig_path, "gap-$(dist_labels[i]).pdf"))
end


## Heavy distributions
heavy_dist_plt = plot()
heavy_dists = [Levy(0, 1), Pareto(0.5)]
heavy_dist_labels = ["Levy" "Pareto"]
for (i, heavy_dists) in enumerate(heavy_dists)
    plot!(heavy_dist_plt, x, pdf.(heavy_dists, x), label="$(heavy_dist_labels[i])", lw=2)
end
display(heavy_dist_plt)
savefig(heavy_dist_plt, joinpath(fig_path, "heavy_dists.pdf"))

# Run trials
gaps_heavy = zeros(T, length(ks), length(heavy_dists))
q_ups_heavy = zeros(T, length(ks), length(heavy_dists))
q_lows_heavy = zeros(T, length(ks), length(heavy_dists))
for (i, val_dist) in enumerate(heavy_dists), (j, k) in enumerate(ks)
    ret = simulate_block(val_dist, tx_size_dist, N, k, T)
    gaps_heavy[:, j, i] .= ret[1]
    q_ups_heavy[:, j, i] .= ret[2]
    q_lows_heavy[:, j, i] .= ret[3]
end

μs_heavy = mean(gaps_heavy, dims=1)[1,:,:]
meds_heavy = median(gaps_heavy, dims=1)[1,:,:]
q_ups_heavy = mean(q_ups_heavy, dims=1)[1,:,:]
q_lows_heavy = mean(q_lows_heavy, dims=1)[1,:,:]

# Plot results
plts = []
colors = [:red :mediumblue]
for i in 1:length(heavy_dists)
    q_up = q_ups_heavy[:, i]
    q_low = q_lows_heavy[:, i]
    b1 = @. 1/η * q_up / ( (q_up - q_low) * ks / N + q_low )

    plt = plot(
        ks,
        μs_heavy[:,i],
        lw=3,
        fillalpha=0.5,
        xaxis=:log,
        yaxis=:log,
        # title="Gap between Greedy and FIFO ordering",
        ylabel=L"$p^\star/p^\mathrm{fifo}$",
        xlabel="k = average block size (in # txs)",
        legend=:false,
        dpi=300,
        # minorgrid=true,
        margin=3Plots.PlotMeasures.mm,
        color = colors[i]
    )
    plot!(plt, ks, meds_heavy[:,i], lw=3, label="Gap bound", ls=:dot, color=colors[i])
    plot!(plt, ks, b1, lw=2, label="Gap bound", ls=:dash, color=colors[i])
    
    push!(plts, plt)
end

[display(plt) for plt in plts]
for (i, plt) in enumerate(plts)
    savefig(plt, joinpath(fig_path, "gap-$(heavy_dist_labels[i]).pdf"))
end
