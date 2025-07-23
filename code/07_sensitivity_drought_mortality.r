
using Pkg

Pkg.activate("../")

using Plots, DataFrames, Distributions,
    SpecialFunctions, NLsolve, CSV, Random,
    ForwardDiff, ProgressBars,
    QuadGK, JuMP, Ipopt, StatsBase, Infiltrator,
    IterTools

Ninit = 1.0

## ecological parameters
μ::Float64 = 0.11   ## mortality rate
E::Float64 = 0.5   ## evapotranspiration rate
l::Float64 = 1.5   ## leaf area allometric constant
b::Float64 = 3.0   ## biomass allometric constant
F::Float64 = 10.0  ## fecundity per unit biomass
W₀::Float64 = 0.4  ## initial water content (default)
θ_fc::Float64 = 0.4 ## field capacity
P::Float64 = 10.0

## include function headers
include("simulator/utility_functions.jl")
include("simulator/simulation_functions.jl")
include("simulator/eq_functions.jl")
include("simulator/meta_functions.jl")

## set plotting variables
theme(:default)
my_cgrad = cgrad(:roma)

cg = cgrad(["#253494", "#41b6c4", "#a1dab4"], 3, rev = false)

sdnew = generate_spp_data(6, 0.7, 1, 1.0 / ((15 - 1) / 2), F, μ,
                          3.0, 0.4, 0.0, 0.25, 0.03, 1.5, 5.0, 2.0, 4.0, 0.9, 1.1, true)
psi = -10.0 .+ (4.5 ./ sdnew.δ) .^ 0.3
sdnew.Wᵢ .= w_psi.(psi)

sdnew


sim = sim_water_ppa(sdnew, 500, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    false, 1000, collect(range(0.001, 0.03, length = nrow(sdnew))),
                    100.0)

sim1000 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 500, collect(range(0.001, 0.03, length = nrow(sdnew))),
                    100.0)

sim500 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 250, collect(range(0.001, 0.03, length = nrow(sdnew))),
                    100.0)

sim100 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 100, collect(range(0.001, 0.03, length = nrow(sdnew))),
                    100.0)

sim50 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 50, collect(range(0.001, 0.03, length = nrow(sdnew))),
                    100.0)

sim25 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 25, collect(range(0.001, 0.03, length = nrow(sdnew))),
                    100.0)


p1 = plot_simulation_dynamics(sim, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-2.0, 6.0], colorbar = :false, legend = :false, xlab = "")
#annotate!(350,5.5, "No drought", 10)

p2 = plot_simulation_dynamics(sim1000, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-2.0, 6.0], colorbar = :false, xlab = "", ylab = "", legend = :false)
#annotate!(250,5.5, "100 yrs", 10)

p3 = plot_simulation_dynamics(sim500, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-2.0, 6.0], colorbar = :false, xlab = "", ylab = "", legend = :false)
#annotate!(250,5.5, "50 yrs", 10)

p4 = plot_simulation_dynamics(sim100, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-10.0, 6.0], colorbar = :false, legend = :false)
#annotate!(250,5.5, "20 yrs", 10)

p5 = plot_simulation_dynamics(sim50, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-10.0, 6.0], colorbar = :false, ylab = "", legend = :false)
#annotate!(250,5.5, "10 years", 10)

p6 = plot_simulation_dynamics(sim25, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-10.0, 6.0], colorbar = :false, ylab = "", legend = :false)
#annotate!(250,5.5, "5 years", 10)


plot(p1, p2, p3, p4, p5, p6, layout = (2,3))

savefig("../figures/figureS6.pdf")





##---------------------------------------------------------------
## Window size
##---------------------------------------------------------------



function calc_lower(R, ϵ, HSM, w₁ = 0.1, w₂ = 0.07, τ₁ = 0.05, τ₂ = 0.09)
    ((w₁ - w₂) / (R * (1 - (τ₁ / τ₂))))
end

function calc_upper(R, ϵ, HSM, w₁ = 0.1, w₂ = 0.07, τ₁ = 0.05, τ₂ = 0.09)
    ((log(w₂) - log(w₁ - HSM)) / (ϵ * (1 - τ₁ / 0.1)))
end


function calc_window(R, ϵ, HSM, w₁ = 0.1, w₂ = 0.07, τ₁ = 0.05, τ₂ = 0.09)
    ((log(w₂) - log(w₁ - HSM)) / (ϵ * (1 - τ₁ / 0.1))) - ((w₁ - w₂) / (R * (1 - (τ₁ / τ₂))))
end

## HSM
HSM =  [0.0:0.001:0.099;]
window = calc_window.(1.0, 10.0, HSM)
lower = calc_lower.(1.0, 10.0, HSM)
upper = calc_upper.(1.0, 10.0, HSM)

plot(HSM, window, frame = :box, linewidth = 0, color = :black, legend = :none)
Plots.abline!(0, 0, line = :dash, color = :gray, linewidth = 2)
hsm1 = plot!(HSM, window, frame = :box, linewidth = 10, color = :black, legend = :none,ylim = [-0.2, 0.9])

plot(HSM, lower, frame = :box, linewidth = 10, color = "#33a02c", legend = :none)
hsm2 = plot!(HSM, upper, linewidth = 10, color = "#1f78b4",ylim = [-0.2, 0.9])

plot(hsm2, hsm1, size = (1000, 500))
savefig("../figures/figureS9/figureS9a.pdf")


## evaporation
ϵ = [0.01:0.01:20.0;]
window = calc_window.(1.0, ϵ, 0.04)
lower = calc_lower.(1.0, ϵ, 0.04)
upper = calc_upper.(1.0, ϵ, 0.04)

plot(ϵ, window, frame = :box, linewidth = 0, color = :black, legend = :none)
Plots.abline!(0, 0, line = :dash, color = :gray, linewidth = 2)
evap1 = plot!(ϵ, window, frame = :box, linewidth = 10, color = :black, legend = :none,ylim = [-0.2, 0.9])

plot(ϵ, lower, frame = :box, linewidth = 10, color = "#33a02c", legend = :none)
evap2 = plot!(ϵ, upper, linewidth = 10, color = "#1f78b4",ylim = [-0.2, 0.9])

plot(evap2, evap1, size = (1000, 500))
savefig("../figures/figureS9/figureS9b.pdf")



## rainfall
R = [1.0:0.01:5.0;]
window = calc_window.(R, 10.0, 0.04)
lower = calc_lower.(R, 10.0, 0.04)
upper = calc_upper.(R, 10.0, 0.04)

plot(R, window, frame = :box, linewidth = 0, color = :black, legend = :none)
Plots.abline!(0, 0, line = :dash, color = :gray, linewidth = 2)
R1 = plot!(R, window, frame = :box, linewidth = 10, color = :black, legend = :none)

plot(R, lower, frame = :box, linewidth = 10, color = "#33a02c", legend = :none)
R2 = plot!(R, upper, linewidth = 10, color = "#1f78b4")

plot(R2, R1, size = (1000, 500))
savefig("../figures/figureS9/figureS9c.pdf")




# Define the vectors

HSM = collect(range(0.04 ,0.09, length = 40))
ϵ = collect(range(0.1, 5.0, length = 40))

# Generate all combinations using a list comprehension
combinations = [(h, e) for h in HSM for e in ϵ]

# Create DataFrame
df = DataFrame(HSM = first.(combinations), ϵ = last.(combinations))
df.window .= Vector{Float64}(undef, nrow(df))

for i in 1:nrow(df)
    df[i,:window] = calc_window(1.0, df[i,:ϵ], df[i,:HSM])
end

hsm_values = unique(df.HSM)
ϵ_values = unique(df.ϵ)

# Create a matrix with results
result_matrix = reshape(df.window, (length(ϵ_values), length(hsm_values)))

# Transpose to align with correct axes for visualization
result_matrix = permutedims(result_matrix)

result_matrix[result_matrix .< 0.0] .= NaN

# Plot heatmap
heatmap(hsm_values, ϵ_values, result_matrix', xlabel="HSM", ylabel="ϵ", title="Window size",
        colorbar_scale = :log10)
savefig("../figures/figureS9/S9d.pdf")


# Define the vectors

HSM = collect(range(0.04 ,0.09, length = 40))
R = collect(range(1.0, 5.0, length = 40))

# Generate all combinations using a list comprehension
combinations = [(h, e) for h in HSM for e in R]

# Create DataFrame
df = DataFrame(HSM = first.(combinations), R = last.(combinations))
df.window .= Vector{Float64}(undef, nrow(df))

for i in 1:nrow(df)
    df[i,:window] = calc_window(df[i,:R], 10.0, df[i,:HSM])
end

hsm_values = unique(df.HSM)
R_values = unique(df.R)

# Create a matrix with results
result_matrix = reshape(df.window, (length(R_values), length(hsm_values)))

# Transpose to align with correct axes for visualization
result_matrix = permutedims(result_matrix)

result_matrix[result_matrix .< 0.0] .= NaN

# Plot heatmap
heatmap(hsm_values, ϵ_values, result_matrix', xlabel="HSM", ylabel="Annual rainfall", title="Heatmap of calc_window results",
        colorbar_scale = :log10)
savefig("../figures/figureS9/S9e.pdf")




# Define the vectors

R = collect(range(1.0, 5.0, length = 40))
ϵ = collect(range(0.1, 5.0, length = 40))

# Generate all combinations using a list comprehension
combinations = [(h, e) for h in R for e in ϵ]

# Create DataFrame
df = DataFrame(R = first.(combinations), ϵ = last.(combinations))
df.window .= Vector{Float64}(undef, nrow(df))

for i in 1:nrow(df)
    df[i,:window] = calc_window(df[i,:R], df[i,:ϵ], 0.04)
end

R_values = unique(df.R)
ϵ_values = unique(df.ϵ)

# Create a matrix with results
result_matrix = reshape(df.window, (length(ϵ_values), length(R_values)))

# Transpose to align with correct axes for visualization
result_matrix = permutedims(result_matrix)

result_matrix[result_matrix .< 0.0] .= NaN

# Plot heatmap
heatmap(hsm_values, R_values, result_matrix', xlabel="Annual rainfall", ylabel="ϵ", title="Window size",
        colorbar_scale = :log10)
savefig("../figures/figureS9/S9f.pdf")





##---------------------------------------------------------------
## The dumb way
##---------------------------------------------------------------


function drought_mortality(w::Vector{Float64}, k::Float64 = 1000.0)
    dm = 1 ./ (1 .+ exp.(-k .* (w .- mean(w))))
end

dm = drought_mortality(sdnew.Wᵢ, 2000.0)
plot(sdnew.Wᵢ, dm, seriestype = :scatter)


sim = sim_water_ppa(sdnew, 1000, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    false, 1000, dm)

sim1000 = sim_water_ppa(sdnew, 1000, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 1000, dm)

sim500 = sim_water_ppa(sdnew, 1000, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 500, dm)

sim100 = sim_water_ppa(sdnew, 1000, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 100, dm)

sim50 = sim_water_ppa(sdnew, 1000, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 50, dm)

sim25 = sim_water_ppa(sdnew, 1000, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    true, 25, dm)


p1 = plot_simulation_dynamics(sim, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-1.0, 6.0], colorbar = :false, legend = :false, xlab = "")
annotate!(350,5.5, "No drought", 10)

p2 = plot_simulation_dynamics(sim1000, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-1.0, 6.0], colorbar = :false, xlab = "", ylab = "", legend = :false)
annotate!(250,5.5, "200 yrs", 10)

p3 = plot_simulation_dynamics(sim500, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-1.0, 6.0], colorbar = :false, xlab = "", ylab = "", legend = :false)
annotate!(250,5.5, "100 yrs", 10)

p4 = plot_simulation_dynamics(sim100, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-1.0, 6.0], colorbar = :false, legend = :false)
annotate!(250,5.5, "20 yrs", 10)

p5 = plot_simulation_dynamics(sim50, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-1.0, 6.0], colorbar = :false, ylab = "", legend = :false)
annotate!(250,5.5, "10 years", 10)

p6 = plot_simulation_dynamics(sim25, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-1.0, 6.0], colorbar = :false, ylab = "", legend = :false)
annotate!(250,5.5, "5 years", 10)


plot(p1, p2, p3, p4, p5, p6, layout = (2,3))
savefig("../figures/figureS6.pdf")


##---------------------------------------------------------------
## for many communities
##---------------------------------------------------------------
