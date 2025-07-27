
using Pkg

Pkg.activate("../../")

using Plots, DataFrames, Distributions,
    SpecialFunctions, Random,ProgressBars

cd("../")

Ninit = 1.0

## ecological parameters
μ::Float64 = 0.11   ## mortality rate
E::Float64 = 0.5   ## evapotranspiration rate
l::Float64 = 1.5   ## leaf area allometric constant
b::Float64 = 3.0   ## biomass allometric constant
F::Float64 = 3.0  ## fecundity per unit biomass
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

function calc_τ(Aₘ::Float64, r::Float64, α::Float64, δ::Float64,
                F::Float64, μ::Float64, T::Float64, b::Float64 = 3.0)
    (T / Aₘ) * (α * δ * ((μ^b)/(F * gamma(b)))^(1/(b-1)) + r)
end

total_list = collect(range(0.2, stop = 0.6, length = 5));
P_list = collect(range(2.0, stop = 10.0, length = 50));
param = reshape(collect(Base.product(total_list, P_list)), (length(total_list) * length(P_list), 1));

W₀_list = Vector{Float64}(undef, length(param))
for i in 1:length(param)
    W₀_list[i] = param[i][1] / param[i][2]
end
W₀_list = W₀_list .+ 0.05

function w₁_max(w₀, τ₁, τ₂, w₂)
    w₀ - (τ₁ / τ₂) * (w₀ - w₂)
end

rain = DataFrame(w₀ = W₀_list,
                 total = Vector{Float64}(undef, length(W₀_list)),
                 P = Vector{Float64}(undef, length(W₀_list)))
for i in 1:nrow(rain)
    rain[i,:total] = param[i][1]
    rain[i,:P] = param[i][2]
end

τ₁ = 0.05
τ₂ = 0.075
w₂ = 0.08

rain.w₁_max = Vector{Float64}(undef, nrow(rain))
for i in 1:nrow(rain)
    rain[i,:w₁_max] = w₁_max(rain[i,:w₀], τ₁, τ₂, w₂)
end

p1 = plot(rain.P, rain.w₁_max, group = rain.total, frame = :box, grid = :false, widen = :false,
     linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true), line_z = rain.total, colorbar = :false,
          ylim = [minimum(rain.w₁_max), maximum(rain.w₁_max)], xtickfontcolor = :white)

rain.T .= 1 ./ rain.P

δ₁ = 4.0
ζ = 0.2

function δ_max(δ₁, ζ, ρ, E, T)
    (δ₁ * (1 + ζ)) / ((ρ / (E * T)) - ζ)
end

rain.δ_max = Vector{Float64}(undef, nrow(rain))
for i in 1:nrow(rain)
    rain[i,:δ_max] = δ_max(δ₁, ζ, rain[i,:w₀] - 0.05, E, rain[i,:T])
end

p2 = plot(rain.P, rain.δ_max, group = rain.total, frame = :box, grid = :false, widen = :false,
     linewidth = 5, seriescolor = my_cgrad, line_z = rain.total, colorbar = :false,
     plot_title = "Maximum δ")

psi = -10.0 .+ (4.5 ./ rain.δ_max) .^ 0.3
rain.wᵢ_min .= w_psi.(psi, 0.01, 0.05, -10.0, 5.5)

function psi_w(w, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    psimin * ((w - wmin) / (wmax - wmin)) ^ (-1 / lambda)
end


psi2 = psi_w.(rain.w₁_max, 0.01, 0.05, -10.0, 5.5)
δ_min = 4.5 .* (psi2 .+ 10.0) .^ (-1 / 0.3)
rain.δ_min .= δ_min

rain.g_max .= (0.25 .- 0.03) ./ (rain.δ_min .* 1.0)
rain.g_min .= (0.25 .- 0.03) ./ (rain.δ_max .* 1.0)

p3 = plot(rain.P, rain.wᵢ_min, group = rain.total, frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor = my_cgrad, line_z = rain.total, colorbar = :false,
          ylim = [minimum(rain.wᵢ_min), maximum(rain.w₁_max)])
cgrad(:roma, 5, categorical = true)[2]

t = 0.2
p4 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:wᵢ_min],
          ribbon = (rain[rain.total .== t,:wᵢ_min] .- rain[rain.total .== t,:wᵢ_min],
                    rain[rain.total .== t,:w₁_max] .- rain[rain.total .== t,:wᵢ_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true)[1], colorbar = :false, legend = :false,
          ylim = [minimum(rain.wᵢ_min), maximum(rain.w₁_max)], ytickfontcolor = :white, xtickfontcolor = :white)
annotate!(7.5, 0.16, text("ρ = 0.2", 10))

t = 0.3
p5 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:wᵢ_min],
          ribbon = (rain[rain.total .== t,:wᵢ_min] .- rain[rain.total .== t,:wᵢ_min],
                    rain[rain.total .== t,:w₁_max] .- rain[rain.total .== t,:wᵢ_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true)[2], colorbar = :false, legend = :false,
          ylim = [minimum(rain.wᵢ_min), maximum(rain.w₁_max)], ytickfontcolor = :white, xtickfontcolor = :white)
annotate!(7.5, 0.16, text("ρ = 0.3", 10))

t = 0.4
p6 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:wᵢ_min],
          ribbon = (rain[rain.total .== t,:wᵢ_min] .- rain[rain.total .== t,:wᵢ_min],
                    rain[rain.total .== t,:w₁_max] .- rain[rain.total .== t,:wᵢ_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true)[3], colorbar = :false, legend = :false,
          ylim = [minimum(rain.wᵢ_min), maximum(rain.w₁_max)])
annotate!(7.5, 0.16, text("ρ = 0.4", 10))

t = 0.5
p7 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:wᵢ_min],
          ribbon = (rain[rain.total .== t,:wᵢ_min] .- rain[rain.total .== t,:wᵢ_min],
                    rain[rain.total .== t,:w₁_max] .- rain[rain.total .== t,:wᵢ_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true)[4], colorbar = :false, legend = :false,
          ylim = [minimum(rain.wᵢ_min), maximum(rain.w₁_max)], ytickfontcolor = :white)
annotate!(7.5, 0.16, text("ρ = 0.5", 10))

t = 0.6
p8 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:wᵢ_min],
          ribbon = (rain[rain.total .== t,:wᵢ_min] .- rain[rain.total .== t,:wᵢ_min],
                    rain[rain.total .== t,:w₁_max] .- rain[rain.total .== t,:wᵢ_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor =cgrad(:roma, 5, categorical = true)[5], colorbar = :false, legend = :false,
          ylim = [minimum(rain.wᵢ_min), maximum(rain.w₁_max)], ytickfontcolor = :white)
annotate!(7.5, 0.16, text("ρ = 0.6", 10))

plot(p1, p4, p5, p6, p7, p8, layout = (2,3))

savefig("../figures/figureS3/figureS3.pdf")



p1 = plot(rain.P, rain.g_max, group = rain.total, frame = :box, grid = :false, widen = :false,
     linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true), line_z = rain.total, colorbar = :false,
          ylim = [minimum(rain.g_min), maximum(rain.g_max)], xtickfontcolor = :white)

t = 0.2
p4 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:g_min],
          ribbon = (rain[rain.total .== t,:g_min] .- rain[rain.total .== t,:g_min],
                    rain[rain.total .== t,:g_max] .- rain[rain.total .== t,:g_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true)[1], colorbar = :false, legend = :false,
          ylim = [minimum(rain.g_min), maximum(rain.g_max)], ytickfontcolor = :white, xtickfontcolor = :white)
annotate!(7.5, 0.6, text("ρ = 0.2", 10))

t = 0.3
p5 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:g_min],
          ribbon = (rain[rain.total .== t,:g_min] .- rain[rain.total .== t,:g_min],
                    rain[rain.total .== t,:g_max] .- rain[rain.total .== t,:g_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true)[2], colorbar = :false, legend = :false,
          ylim = [minimum(rain.g_min), maximum(rain.g_max)], ytickfontcolor = :white, xtickfontcolor = :white)
annotate!(7.5, 0.6, text("ρ = 0.3", 10))

t = 0.4
p6 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:g_min],
          ribbon = (rain[rain.total .== t,:g_min] .- rain[rain.total .== t,:g_min],
                    rain[rain.total .== t,:g_max] .- rain[rain.total .== t,:g_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true)[3], colorbar = :false, legend = :false,
          ylim = [minimum(rain.g_min), maximum(rain.g_max)])
annotate!(7.5, 0.6, text("ρ = 0.4", 10))

t = 0.5
p7 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:g_min],
          ribbon = (rain[rain.total .== t,:g_min] .- rain[rain.total .== t,:g_min],
                    rain[rain.total .== t,:g_max] .- rain[rain.total .== t,:g_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor = cgrad(:roma, 5, categorical = true)[4], colorbar = :false, legend = :false,
          ylim = [minimum(rain.g_min), maximum(rain.g_max)], ytickfontcolor = :white)
annotate!(7.5, 0.6, text("ρ = 0.5", 10))

t = 0.6
p8 = plot(rain[rain.total .== t,:P], rain[rain.total .== t,:g_min],
          ribbon = (rain[rain.total .== t,:g_min] .- rain[rain.total .== t,:g_min],
                    rain[rain.total .== t,:g_max] .- rain[rain.total .== t,:g_min]),
          frame = :box, grid = :false, widen = :false,
          linewidth = 5, seriescolor =cgrad(:roma, 5, categorical = true)[5], colorbar = :false, legend = :false,
          ylim = [minimum(rain.g_min), maximum(rain.g_max)], ytickfontcolor = :white)
annotate!(7.5, 0.6, text("ρ = 0.6", 10))

plot(p1, p4, p5, p6, p7, p8, layout = (2,3))
savefig("../figures/figureS4/figureS4.pdf")
