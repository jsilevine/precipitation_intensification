##---------------------------------------------------------------
## Figure S5: Feasibility window calculations for HSM, Evaporation, and Rainfall
##            Accounting for Drought Mortality -- heatmaps
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

## Load required packages
using Pkg
Pkg.activate("../../")  # Activate the project environment

## Import libraries for plotting, data manipulation, and numerical methods
using DataFrames, Plots

cd("../")

## Set plotting variables
theme(:default);  ## Set default plotting theme

##---------------------------------------------------------------
## 01. Define functions for calculating window size
##---------------------------------------------------------------

## Function to calculate the lower bound of the window
function calc_lower(R, ϵ, HSM, w₁ = 0.1, w₂ = 0.07, τ₁ = 0.05, τ₂ = 0.09)
    ((w₁ - w₂) / (R * (1 - (τ₁ / τ₂))))
end

## Function to calculate the upper bound of the window
function calc_upper(R, ϵ, HSM, w₁ = 0.1, w₂ = 0.07, τ₁ = 0.05, τ₂ = 0.09)
    ((log(w₂) - log(w₁ - HSM)) / (ϵ * (1 - τ₁ / 0.1)))
end

## Function to calculate the window size
function calc_window(R, ϵ, HSM, w₁ = 0.1, w₂ = 0.07, τ₁ = 0.05, τ₂ = 0.09)
    calc_upper(R, ϵ, HSM, w₁, w₂, τ₁, τ₂) - calc_lower(R, ϵ, HSM, w₁, w₂, τ₁, τ₂)
end

##---------------------------------------------------------------
## 02. Window Size Calculations for HSM
##---------------------------------------------------------------

## Define HSM values
HSM = [0.0:0.001:0.099;]

## Calculate window size, lower bound, and upper bound
window = calc_window.(1.0, 10.0, HSM)
lower = calc_lower.(1.0, 10.0, HSM)
upper = calc_upper.(1.0, 10.0, HSM)

## Plot window size
plot(HSM, window, frame = :box, linewidth = 0, color = :black, legend = :none)
Plots.abline!(0, 0, line = :dash, color = :gray, linewidth = 2)
hsm1 = plot!(HSM, window, frame = :box, linewidth = 10, color = :black, legend = :none, ylim = [-0.2, 0.9])

## Plot lower and upper bounds
plot(HSM, lower, frame = :box, linewidth = 10, color = "#33a02c", legend = :none)
hsm2 = plot!(HSM, upper, linewidth = 10, color = "#1f78b4", ylim = [-0.2, 0.9])

## Combine plots and save figure
plot(hsm2, hsm1, size = (1000, 500))
savefig("../figures/figureS5/figureS5a.pdf")

##---------------------------------------------------------------
## 03. Window Size Calculations for Evaporation
##---------------------------------------------------------------

## Define evaporation values
ϵ = [0.01:0.01:20.0;]

## Calculate window size, lower bound, and upper bound
window = calc_window.(1.0, ϵ, 0.04)
lower = calc_lower.(1.0, ϵ, 0.04)
upper = calc_upper.(1.0, ϵ, 0.04)

## Plot window size
plot(ϵ, window, frame = :box, linewidth = 0, color = :black, legend = :none)
Plots.abline!(0, 0, line = :dash, color = :gray, linewidth = 2)
evap1 = plot!(ϵ, window, frame = :box, linewidth = 10, color = :black, legend = :none, ylim = [-0.2, 0.9])

## Plot lower and upper bounds
plot(ϵ, lower, frame = :box, linewidth = 10, color = "#33a02c", legend = :none)
evap2 = plot!(ϵ, upper, linewidth = 10, color = "#1f78b4", ylim = [-0.2, 0.9])

## Combine plots and save figure
plot(evap2, evap1, size = (1000, 500))
savefig("../figures/figureS5/figureS5b.pdf")

##---------------------------------------------------------------
## 04. Window Size Calculations for Rainfall
##---------------------------------------------------------------

## Define rainfall values
R = [1.0:0.01:5.0;]

## Calculate window size, lower bound, and upper bound
window = calc_window.(R, 10.0, 0.04)
lower = calc_lower.(R, 10.0, 0.04)
upper = calc_upper.(R, 10.0, 0.04)

## Plot window size
plot(R, window, frame = :box, linewidth = 0, color = :black, legend = :none)
Plots.abline!(0, 0, line = :dash, color = :gray, linewidth = 2)
R1 = plot!(R, window, frame = :box, linewidth = 10, color = :black, legend = :none)

## Plot lower and upper bounds
plot(R, lower, frame = :box, linewidth = 10, color = "#33a02c", legend = :none)
R2 = plot!(R, upper, linewidth = 10, color = "#1f78b4")

## Combine plots and save figure
plot(R2, R1, size = (1000, 500))
savefig("../figures/figureS5/figureS5c.pdf")
