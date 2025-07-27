##---------------------------------------------------------------
## Figure S6: Feasibility window calculations for HSM, Evaporation, and Rainfall
##            Accounting for Drought Mortality
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

## Define the vectors for HSM and ϵ
HSM = collect(range(0.04, 0.09, length=40))  ## HSM values from 0.04 to 0.09
ϵ = collect(range(0.1, 5.0, length=40))      ## ϵ values from 0.1 to 5.0

## Generate all combinations of HSM and ϵ
combinations = [(h, e) for h in HSM for e in ϵ]

## Create a DataFrame to store combinations and results
df = DataFrame(HSM=first.(combinations), ϵ=last.(combinations))
df.window .= Vector{Float64}(undef, nrow(df))  ## Initialize the window column

## Calculate window values for each combination
for i in 1:nrow(df)
    df[i, :window] = calc_window(1.0, df[i, :ϵ], df[i, :HSM])  ## Example calculation
end

## Extract unique values for HSM and ϵ
hsm_values = unique(df.HSM)
ϵ_values = unique(df.ϵ)

## Reshape the results into a matrix for visualization
result_matrix = reshape(df.window, (length(ϵ_values), length(hsm_values)))
result_matrix = permutedims(result_matrix)  ## Transpose for correct axes

## Replace negative values with NaN
result_matrix[result_matrix .< 0.0] .= NaN

## Plot heatmap
heatmap(hsm_values, ϵ_values, result_matrix', xlabel="HSM", ylabel="ϵ", title="Window size",
        colorbar_scale=:log10)

## Save the figure
savefig("../figures/figureS6/S6a.pdf")

##---------------------------------------------------------------
## 03. Window Size Calculations for Evaporation
##---------------------------------------------------------------

## Define the vectors for HSM and R
HSM = collect(range(0.04, 0.09, length=40))  ## HSM values from 0.04 to 0.09
R = collect(range(1.0, 5.0, length=40))      ## R values from 1.0 to 5.0

## Generate all combinations of HSM and R
combinations = [(h, e) for h in HSM for e in R]

## Create a DataFrame to store combinations and results
df = DataFrame(HSM=first.(combinations), R=last.(combinations))
df.window .= Vector{Float64}(undef, nrow(df))  ## Initialize the window column

## Calculate window values for each combination
for i in 1:nrow(df)
    df[i, :window] = calc_window(df[i, :R], 10.0, df[i, :HSM])  ## Example calculation
end

## Extract unique values for HSM and R
hsm_values = unique(df.HSM)
R_values = unique(df.R)

## Reshape the results into a matrix for visualization
result_matrix = reshape(df.window, (length(R_values), length(hsm_values)))
result_matrix = permutedims(result_matrix)  ## Transpose for correct axes

## Replace negative values with NaN
result_matrix[result_matrix .< 0.0] .= NaN

## Plot heatmap
heatmap(hsm_values, R_values, result_matrix', xlabel="HSM", ylabel="Annual rainfall",
        title="Heatmap of calc_window results", colorbar_scale=:log10)

## Save the figure
savefig("../figures/figureS6/S6b.pdf")

##---------------------------------------------------------------
## 04. Window Size Calculations for Rainfall
##---------------------------------------------------------------

## Define the vectors for R and ϵ
R = collect(range(1.0, 5.0, length=40))      ## R values from 1.0 to 5.0
ϵ = collect(range(0.1, 5.0, length=40))      ## ϵ values from 0.1 to 5.0

## Generate all combinations of R and ϵ
combinations = [(h, e) for h in R for e in ϵ]

## Create a DataFrame to store combinations and results
df = DataFrame(R=first.(combinations), ϵ=last.(combinations))
df.window .= Vector{Float64}(undef, nrow(df))  ## Initialize the window column

## Calculate window values for each combination
for i in 1:nrow(df)
    df[i, :window] = calc_window(df[i, :R], df[i, :ϵ], 0.04)  ## Example calculation
end

## Extract unique values for R and ϵ
R_values = unique(df.R)
ϵ_values = unique(df.ϵ)

## Reshape the results into a matrix for visualization
result_matrix = reshape(df.window, (length(ϵ_values), length(R_values)))
result_matrix = permutedims(result_matrix)  ## Transpose for correct axes

## Replace negative values with NaN
result_matrix[result_matrix .< 0.0] .= NaN

## Plot heatmap
heatmap(R_values, ϵ_values, result_matrix', xlabel="Annual rainfall", ylabel="ϵ",
        title="Window size", colorbar_scale=:log10)

## Save the figure
savefig("../figures/figureS6/S6c.pdf")
