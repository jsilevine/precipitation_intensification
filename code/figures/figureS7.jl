## Load required packages
using Pkg
Pkg.activate("../../")  # Activate the project environment

## Import libraries for plotting, data manipulation, and numerical methods
using Plots, DataFrames, Distributions, Random, SpecialFunctions, ProgressBars

cd("../")

## Initialize simulation parameters
Ninit = 1.0  # Initial population size

## Ecological parameters
μ::Float64 = 0.11   ## Mortality rate
E::Float64 = 0.5    ## Evapotranspiration rate
l::Float64 = 1.5    ## Leaf area allometric constant
b::Float64 = 3.0    ## Biomass allometric constant
F::Float64 = 10.0   ## Fecundity per unit biomass
W₀::Float64 = 0.4   ## Initial water content (default)
θ_fc::Float64 = 0.4 ## Field capacity
P::Float64 = 10.0   ## Precipitation rate

## Include function headers for simulations and utilities
include("simulator/utility_functions.jl")
include("simulator/simulation_functions.jl")
include("simulator/eq_functions.jl")
include("simulator/meta_functions.jl")

## Set plotting variables
theme(:default)  ## Set default plotting theme
my_cgrad = cgrad(:roma)  ## Define custom color gradient
cg = cgrad(["#253494", "#41b6c4", "#a1dab4"], 3, rev = false)  ## Another color gradient

## Generate species data for simulation
sdnew = generate_spp_data(6, 0.7, 1, 1.0 / ((15 - 1) / 2), F, μ,
                          3.0, 0.4, 0.0, 0.25, 0.03, 1.5, 5.0, 2.0, 4.0, 0.9, 1.1, true)

## Calculate initial water potential and update species data
psi = -10.0 .+ (4.5 ./ sdnew.δ) .^ 0.3
sdnew.Wᵢ .= w_psi.(psi)

## Run simulations for different drought durations
sim = sim_water_ppa(sdnew, 500, nrow(sdnew),
                    1.0, μ, F,
                    5.0, 0.25, 0.2, zeros(1,1),
                    true, 0.2, b, 0.1,
                    false, 0.5, 0.9,
                    false, 0.7,
                    false, 0.7,
                    false, 1000, collect(range(0.001, 0.03, length = nrow(sdnew))),
                    100.0);

sim1000 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                        1.0, μ, F,
                        5.0, 0.25, 0.2, zeros(1,1),
                        true, 0.2, b, 0.1,
                        false, 0.5, 0.9,
                        false, 0.7,
                        false, 0.7,
                        true, 500, collect(range(0.001, 0.03, length = nrow(sdnew))),
                        100.0);

sim500 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                       1.0, μ, F,
                       5.0, 0.25, 0.2, zeros(1,1),
                       true, 0.2, b, 0.1,
                       false, 0.5, 0.9,
                       false, 0.7,
                       false, 0.7,
                       true, 250, collect(range(0.001, 0.03, length = nrow(sdnew))),
                       100.0);

sim100 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                       1.0, μ, F,
                       5.0, 0.25, 0.2, zeros(1,1),
                       true, 0.2, b, 0.1,
                       false, 0.5, 0.9,
                       false, 0.7,
                       false, 0.7,
                       true, 100, collect(range(0.001, 0.03, length = nrow(sdnew))),
                       100.0);

sim50 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                      1.0, μ, F,
                      5.0, 0.25, 0.2, zeros(1,1),
                      true, 0.2, b, 0.1,
                      false, 0.5, 0.9,
                      false, 0.7,
                      false, 0.7,
                      true, 50, collect(range(0.001, 0.03, length = nrow(sdnew))),
                      100.0);

sim25 = sim_water_ppa(sdnew, 500, nrow(sdnew),
                      1.0, μ, F,
                      5.0, 0.25, 0.2, zeros(1,1),
                      true, 0.2, b, 0.1,
                      false, 0.5, 0.9,
                      false, 0.7,
                      false, 0.7,
                      true, 25, collect(range(0.001, 0.03, length = nrow(sdnew))),
                      100.0);

## Plot simulation dynamics for different drought durations
p1 = plot_simulation_dynamics(sim, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-2.0, 6.0], colorbar = :false, legend = :false, xlab = "")

p2 = plot_simulation_dynamics(sim1000, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-2.0, 6.0], colorbar = :false, xlab = "", ylab = "", legend = :false)

p3 = plot_simulation_dynamics(sim500, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-2.0, 6.0], colorbar = :false, xlab = "", ylab = "", legend = :false)

p4 = plot_simulation_dynamics(sim100, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-10.0, 6.0], colorbar = :false, legend = :false)

p5 = plot_simulation_dynamics(sim50, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-10.0, 6.0], colorbar = :false, ylab = "", legend = :false)

p6 = plot_simulation_dynamics(sim25, false, "", my_cgrad, 3.5, true)
plot!(ylim = [-10.0, 6.0], colorbar = :false, ylab = "", legend = :false)

## Combine all plots into a grid layout and save the figure
plot(p1, p2, p3, p4, p5, p6, layout = (2,3))
savefig("../figures/figureS7/figureS7.pdf")
