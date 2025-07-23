##---------------------------------------------------------------
## Code to run theoretical analysis, and produce relevant figures
##
## author: jacob levine; jacoblevine@princeton.edu
##---------------------------------------------------------------

using Pkg

Pkg.activate("../")

using Plots, DataFrames, Distributions,
    SpecialFunctions, NLsolve, CSV, Random,
    ForwardDiff, ProgressBars,
    QuadGK, JuMP, Ipopt, StatsBase, Infiltrator

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
include("../../water_comp_perennials/simulator/water_light/utility_functions.jl")
include("../../water_comp_perennials/simulator/water_light/simulation_functions.jl")
include("../../water_comp_perennials/simulator/water_light/eq_functions.jl")
include("../../water_comp_perennials/simulator/water_light/meta_functions.jl")

## set plotting variables
theme(:default)
my_cgrad = cgrad(:roma)

cg = cgrad(["#253494", "#41b6c4", "#a1dab4"], 3, rev = false)

meq = multi_eq(Nspp = 30, Niter = 6,
               θ_fc = 0.2, mintotal = 0.2, maxtotal = 0.25,
               lengthtotal = 3,
               minP = 2.0, maxP = 10.0, lengthP = 45,
               F = 3.0, μ = μ, n_hts = 1, uf = 0.1, δ_max = 5.0, psi_exp = 1.0)

smeq = summarize_multi_eq(meq)
sub = smeq[smeq.var .== "n", :]
p2a = plot(sub.P, sub.mean, group = sub.map, line_z = sub.map, ribbon = sub.sd, fill_z = sub.map, linewidth = 5,
          seriescolor = cg, fillalpha = 0.3, grid = :false, colorbar = false, frame = :box, xlims = [1.0, 15.0], ylims = [0.0, 30.0],
          legend = :outerright, legendtitle = "Monthly \n precip.", legendtitlefontsize = 7,
          xlab = "Storm frequency (storms per month)", ylab = "Species richness")

## initialize data
Niter = meq[2]; Nspp = meq[4]; tmp = Vector{Int64}; Npar = size(meq[1])[2];
summar = DataFrame(map = repeat(meq[3][:,3], outer = 4),
                    P = repeat(meq[3][:,2], outer = 4),
                    var = repeat(["n", "min", "max", "avg"], inner = Npar),
                    mean = Vector{Float64}(undef, Npar*4),
                    sd = Vector{Float64}(undef, Npar*4))

nfeas_temp = Vector{Float64}(undef, Niter); maxdelta_temp = Vector{Float64}(undef, Niter);
mindelta_temp = Vector{Float64}(undef, Niter); avgdelta_temp = Vector{Float64}(undef, Niter);

for i in 1:size(meq[1])[2]
    for j in 1:Niter
        tmp = findall(meq[1][[(Nspp*(j-1))+1:1:Nspp*j;],i] .> 0.0)
        maxdelta_temp[j] = maximum(meq[7][j][tmp,:Wᵢ])
        mindelta_temp[j] = minimum(meq[7][j][tmp,:Wᵢ])
        avgdelta_temp[j] = mean(meq[7][j][tmp,:Wᵢ])
        nfeas_temp[j] = length(tmp)
    end
    summar[i, :mean] = mean(nfeas_temp); summar[i, :sd] = std(nfeas_temp)
    summar[Npar+i, :mean] = mean(mindelta_temp); summar[Npar+i, :sd] = std(mindelta_temp)
    summar[Npar*2+i, :mean] = mean(maxdelta_temp); summar[Npar*2+i, :sd] = std(maxdelta_temp)
    summar[Npar*3+i, :mean] = mean(avgdelta_temp); summar[Npar*3+i, :sd] = std(avgdelta_temp)
end


sub = summar[summar.var .== "min", :]
sub1 = summar[summar.var .== "max", :]
sub2 = summar[summar.var .== "avg", :]

avg = sub2[sub2.map .== unique(sub2.map)[1], :mean]
minim = sub[sub.map .== unique(sub.map)[1], :mean]
maxim = sub1[sub1.map .== unique(sub1.map)[1], :mean]

p2b = plot(sub[sub.map .== unique(sub.map)[1], :P], avg, ribbon = (avg .- minim, maxim .- avg),
           xlim = [2.0, 10.0], ylim = [0.057, 0.086], legend = :false, frame = :box, grid = :false, widen = :false,
           linewidth = 5,fill = "#253494", color = "#253494",  fillalpha = 0.5)
#yflip!()

avg = sub2[sub2.map .== unique(sub2.map)[2], :mean]
minim = sub[sub.map .== unique(sub.map)[2], :mean]
maxim = sub1[sub1.map .== unique(sub1.map)[2], :mean]

p2c = plot(sub[sub.map .== unique(sub.map)[2], :P], avg, ribbon = (avg .- minim, maxim .- avg),
           xlim = [2.0, 10.0], ylim = [0.057, 0.086], legend = :false, frame = :box, grid = :false, widen = :false,
           linewidth = 5, fill = "#41b6c4", color = "#41b6c4", fillalpha = 0.5)
#yflip!()

avg = sub2[sub2.map .== unique(sub2.map)[3], :mean]
minim = sub[sub.map .== unique(sub.map)[3], :mean]
maxim = sub1[sub1.map .== unique(sub1.map)[3], :mean]

p2d = plot(sub[sub.map .== unique(sub.map)[3], :P], avg, ribbon = (avg .- minim, maxim .- avg),
           xlim = [2.0, 10.0], ylim = [0.057, 0.086], legend = :false, frame = :box, grid = :false, widen = :false,
           linewidth = 5, fill = "#a1dab4", color = "#a1dab4",fillalpha = 0.5)
#yflip!()

panel2 = plot(p2b, p2c, p2d, layout = (1,3), size = (1200, 400))

savefig("../figures/figure3a.pdf")

## initialize data
Niter = meq[2]; Nspp = meq[4]; tmp = Vector{Int64}; Npar = size(meq[1])[2];
summar = DataFrame(map = repeat(meq[3][:,3], outer = 4),
                    P = repeat(meq[3][:,2], outer = 4),
                    var = repeat(["n", "min", "max", "avg"], inner = Npar),
                    mean = Vector{Float64}(undef, Npar*4),
                    sd = Vector{Float64}(undef, Npar*4))

nfeas_temp = Vector{Float64}(undef, Niter); maxdelta_temp = Vector{Float64}(undef, Niter);
mindelta_temp = Vector{Float64}(undef, Niter); avgdelta_temp = Vector{Float64}(undef, Niter);

for i in 1:size(meq[1])[2]
    for j in 1:Niter
        tmp = findall(meq[1][[(Nspp*(j-1))+1:1:Nspp*j;],i] .> 0.0)
        tt = (meq[7][j][tmp,:Aₘ] .- meq[7][j][tmp,:r]) ./ (meq[7][j][tmp,:δ] .* meq[7][j][tmp,:α])
        maxdelta_temp[j] = maximum(tt)
        mindelta_temp[j] = minimum(tt)
        avgdelta_temp[j] = mean(tt)
        nfeas_temp[j] = length(tmp)
    end
    summar[i, :mean] = mean(nfeas_temp); summar[i, :sd] = std(nfeas_temp)
    summar[Npar+i, :mean] = mean(mindelta_temp); summar[Npar+i, :sd] = std(mindelta_temp)
    summar[Npar*2+i, :mean] = mean(maxdelta_temp); summar[Npar*2+i, :sd] = std(maxdelta_temp)
    summar[Npar*3+i, :mean] = mean(avgdelta_temp); summar[Npar*3+i, :sd] = std(avgdelta_temp)
end

sub = summar[summar.var .== "min", :]
sub1 = summar[summar.var .== "max", :]
sub2 = summar[summar.var .== "avg", :]

avg = sub2[sub2.map .== unique(sub2.map)[1], :mean]
minim = sub[sub.map .== unique(sub.map)[1], :mean]
maxim = sub1[sub1.map .== unique(sub1.map)[1], :mean]

p2b = plot(sub[sub.map .== unique(sub.map)[1], :P], avg, ribbon = (avg .- minim, maxim .- avg),
           xlim = [2.0, 10.0], ylim = [minimum(sub.mean), maximum(sub1.mean)], legend = :false, frame = :box, grid = :false, widen = :false,
           linewidth = 5,fill = "#253494", color = "#253494",  fillalpha = 0.5)

avg = sub2[sub2.map .== unique(sub2.map)[2], :mean]
minim = sub[sub.map .== unique(sub.map)[2], :mean]
maxim = sub1[sub1.map .== unique(sub1.map)[2], :mean]

p2c = plot(sub[sub.map .== unique(sub.map)[2], :P], avg, ribbon = (avg .- minim, maxim .- avg),
           xlim = [2.0, 10.0], ylim = [minimum(sub.mean), maximum(sub1.mean)], legend = :false, frame = :box, grid = :false, widen = :false,
           linewidth = 5, fill = "#41b6c4", color = "#41b6c4", fillalpha = 0.5)

avg = sub2[sub2.map .== unique(sub2.map)[3], :mean]
minim = sub[sub.map .== unique(sub.map)[3], :mean]
maxim = sub1[sub1.map .== unique(sub1.map)[3], :mean]

p2d = plot(sub[sub.map .== unique(sub.map)[3], :P], avg, ribbon = (avg .- minim, maxim .- avg),
           xlim = [2.0, 10.0], ylim = [minimum(sub.mean), maximum(sub1.mean)], legend = :false, frame = :box, grid = :false, widen = :false,
           linewidth = 5, fill = "#a1dab4", color = "#a1dab4",fillalpha = 0.5)

panel2 = plot(p2b, p2c, p2d, layout = (1,3), size = (1200, 400))

savefig("../figures/figure3b.pdf")


##---------------------------------------------------------------
## Trying something new
##---------------------------------------------------------------

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
          ylim = [minimum(rain.wᵢ_min), maximum(rain.w₁_max)], xtickfontcolor = :white)

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

savefig("../figures/figureS4.pdf")



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
savefig("../figures/figureS5.pdf")


##---------------------------------------------------------------
## Drought mortality
##---------------------------------------------------------------
