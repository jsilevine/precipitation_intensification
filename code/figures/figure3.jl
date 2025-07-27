##---------------------------------------------------------------
## figure3.jl: Theoretical analysis and generation of figure 3
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

using Pkg

Pkg.activate("../../")

using Plots, DataFrames, Distributions,
    SpecialFunctions, NLsolve, CSV, Random,
    ForwardDiff, ProgressBars,
    QuadGK, JuMP, Ipopt, StatsBase, Infiltrator

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

cg = cgrad(["#253494", "#41b6c4", "#a1dab4"], 3, rev = false)

meq = multi_eq(Nspp = 30, Niter = 6,
               θ_fc = 0.2, mintotal = 0.17, maxtotal = 0.22,
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

savefig("../figures/figure3/figure3a.pdf")

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

savefig("../figures/figure3/figure3b.pdf")
