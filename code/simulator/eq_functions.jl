##---------------------------------------------------------------
## WATER_ONLY_SIMULATOR -- eq_functions.jl
##
## By: Jacob Levine -- jacoblevine@princeton.edu
## November 2022
##
## This script contains functions related to calculated equilibrium values
## and determining coexistence outcomes. All functions defined here are primarily
## for finding analytical solutions, rather than performing simulations
##---------------------------------------------------------------

"""
    minfeas(μ::Float64 = 0.1, F::Float64 = 1.0, b::Float64 = 2.5, analytical = false)

Calculates the minimum feasible growth rate, `g`, given a mortality rate `μ`, fecundity rate `F`, and
allometric exponent `b`. If `analytical` is `true`, the denominator of the growth rate expression is
calculated using a gamma function (as in the continuous time model). If `analytical` is `false`, the
denominator is calculated using the polylogarithm function (as in the discrete time model).
"""
function minfeas(μ::Float64 = 0.1, F::Float64 = 1.0, b::Float64 = 2.5, analytical = false)

    if analytical
        s = gamma(b, 0.0) * μ^-(b)
    else
        s = sigma(μ, b)
    end

   (1/(F * s))^(1/(b-1))

end


"""
    calc_eq_leaf_area(eqN::Vector{Float64}, F::Float64, μ::Float64)

Calculates the equilibrium leaf area given a vector of equilibrium abundances `eqN.`
"""
function calc_eq_leaf_area(eqN::Vector{Float64}, F::Float64, μ::Float64)
    eqN ./ (μ .* F)
end


"""
    μ_feas(spp_data::DataFrame, F::Float64, μ::Float64, analytical::Bool = false)

Determine which species in `spp_data` are feasible in a non-competitive context. Operates by calling
the `minfeas()` function, which calculates the minimum feasible growth rate given the mortality rate, `μ`,
fecundity rate, `F`, and allometric exponent `b`. This function returns a list of indices indicating which
species are feasible. Indices refer to species ordering in `spp_data`.
"""
function μ_feas(spp_data::DataFrame, F::Float64, μ::Float64, analytical::Bool = false)
    mf = minfeas(μ, F, b, true)
    findall((spp_data.Aₘ ./ (spp_data.α .* spp_data.δ) ) .> mf)
end



"""
    calc_ρ_cc(T::Float64, spp_data::DataFrame, F::Float64, μ::Float64, b::Float64)

Calculate the value of `ρ` (the minimum storm size at which the canopy closes) as a function of
the time between storm events `T`, species characteristics `spp_data`, fecundity `F`, mortality rate `μ`,
and allometric exponent `b`.

"""
function calc_ρ_cc(T::Float64, spp_data::DataFrame, F::Float64, μ::Float64, b::Float64)
    E * spp_data.τ[1] + spp_data.Wᵢ[1] - spp_data.Wᵢ[nrow(spp_data)]
end



"""
    calc_T_cc(ρ::Float64, spp_data::DataFrame, F::Float64, μ::Float64, b::Float64)

Calculate the value of `T` (the time between storm events) as a function of
storm size `ρ`, species characteristics `spp_data`, fecundity `F`, mortality rate `μ`,
and allometric exponent `b`.

"""
function calc_T_cc(ρ::Float64, spp_data::DataFrame, F::Float64, μ::Float64, b::Float64)

    function st(T)
        [calc_τ(spp_data.C₁[1], spp_data.C₂[1], F, μ, T[1], b) - (spp_data.Wᵢ[nrow(spp_data)] + ρ - spp_data.Wᵢ[1]) / E]
    end

    nlsolve(st, [0.1])
end



"""
    calc_ρ_cx(spp_data::DataFrame, E::Float64, T::Float64)

Calculate the minimum value of `ρ` (storm size) at which a species is
excluded from the community due to capped transpiration.

"""
function calc_ρ_cx(wᵢ::Float64, δᵢ::Float64, w₁::Float64, δ₁::Float64,
                   Aₘ::Float64, r::Float64, E::Float64, T::Float64)
    E * (T * (δ₁ / δᵢ) * (1 - (r / Aₘ)) - r) + w₁ - wᵢ
end



"""
    calc_T_cx(spp_data::DataFrame, E::Float64, T::Float64)

Calculate the maximum value of `T` (storm interval length) at which a species
is excluded from the community due to capped transpiration.

"""
function calc_T_cx(wᵢ::Float64, δᵢ::Float64, w₁::Float64, δ₁::Float64,
                   Aₘ::Float64, r::Float64, E::Float64, ρ::Float64)
    ((Aₘ * δᵢ) / ((1 - r) * δ₁)) * ((wᵢ + ρ - w₁) / ϵ + r)
end


## assumes all species in spp_data are feasible.
function est_hstar(spp_data::DataFrame, ρ::Float64, F::Float64, μ::Float64,
                   T::Float64, uf::Float64, E::Float64, b::Float64 = 2.5, θ_fc = 0.4)

    t = (minimum([spp_data[nrow(spp_data), :Wᵢ] + ρ, θ_fc]) - spp_data[1,:Wᵢ]) / E
    g = calc_g(t, spp_data.Aₘ[1], spp_data.r[1], spp_data.α[1], spp_data.δ[1], T)
    spp_data.α[1] * sqrt((1 / (μ / (uf * g) - μ / (g))) * log((F * g^(b-1) * gamma(b)) / μ^b))

end

## assumes all species in spp_data are feasible.
function calc_hstar(spp_data::DataFrame, ρ::Float64, F::Float64, μ::Float64,
                   T::Float64, uf::Float64, E::Float64, b::Float64 = 2.5, θ_fc = 0.4)

    t = (minimum([spp_data[nrow(spp_data), :Wᵢ] + ρ, θ_fc]) - spp_data[1,:Wᵢ]) / E
    if t > T
        t = T
    end
    g = calc_g(t, spp_data.Aₘ[1], spp_data.r[1], spp_data.α[1], spp_data.δ[1], T)

    if g < 0
        return 0
    else

        function f(h)

            xs = (h[1] / spp_data.α[1]) .^ 2

            function int(τ)
                exp.(-μ .* τ) .* (xs .+ g .* τ) .^ (b-1)
            end

            [F * exp.(-μ .* (xs ./ (uf .* g))) .* quadgk(int, 0.0, 1e4)[1] .- 1.0]

        end

        nlsolve(f, [0.1]).zero[1]

    end
end

function est_ess_α(ρ::Float64, T::Float64, w₁::Float64, w₂::Float64, δ::Float64, Aₘ::Float64,
                   r::Float64, b::Float64, E::Float64, F::Float64, μ::Float64)

    C = F * ((1 / T) * (Aₘ*((w₂ + ρ - w₁ )/ E) + r))^(b-1) * gamma(b)
    (C^(1 / (b-1)) * exp(-1/2)) / (δ * μ^(b / (b-1)))

end


function calc_τ_cc(hs::Float64, T::Float64, Aₘ::Float64, r::Float64,
                   α::Float64, δ::Float64, μ::Float64, uf::Float64)

    function f(g)

        xs = (hs / α) .^ 2.0

        function int(τ)
            exp.(-μ .* τ) .* (xs .+ g[1] .* τ) .^ 2.0
        end

        [F * exp.(-μ .* (xs ./ (uf .* g[1]))) .* quadgk(int, 0.0, 1e4)[1] .- 1.0]

    end

    g = nlsolve(f, [0.1]).zero[1]
    (T / Aₘ) * (g * α * δ + r)

end

"""
    feas(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04,
    F::Float64 = 10.0, μ::Float64 = 0.1, θ_fc::Float64 = 0.4, uf::Float64 = 0.1,
    analytical::Bool = false, understory_transpiration::Bool = false)

Use convex hull algorithm to determine the feasibility of each species in `sd`. `sd` should be a `DataFrame`
as generated by `generate_spp_data()`. Returns a vector of indices giving the subset of species identifiers
in `sd` which are feasible given their break-even time and critical water content.

If there is a priority effect for early species, function returns two vectors: one indicating all species which
persist absent priority effects, and one which indicates the species which persist when priority effects occur.
"""
function feas(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04,
              F::Float64 = 10.0, μ::Float64 = 0.1, θ_fc::Float64 = 0.4, uf::Float64 = 0.1,
              analytical::Bool = false, understory_transpiration::Bool = false,
              E::Float64 = 0.5, b::Float64 = 3.0)

    W₀ = minimum([minimum(sd[:, :Wᵢ]) + ρ, θ_fc])

    ## first test for feasibility in water-only scenario
    infeas = μ_feas(sd, F, μ, true)
    first_feas = minimum(findall(sd.Wᵢ .< W₀), init = 9999)
    if first_feas == 9999 ## case where no species are feasible
        return []
    else

        ## create dataframe of species characteristics (including w₀) for convex hull algorithm
        points = DataFrame(τ = vcat(sd.τ[infeas], 0.0), Wᵢ = vcat(sd.Wᵢ[infeas], W₀),
                           spp = vcat(sd.spp[infeas], 0))
        ## remove all infeasible species from points dataframe
        points = points[first_feas:nrow(points), :]

        sort!(points, :spp) ## sort in reverse drought tolerance order


        if nrow(points) == 1 ## case where all species infeasible
            return []
        else

            coex = [] ## create empty object to populate with coexisting species indices
            e = 1
            ## loop over species starting with least drought tolerant
            while e < nrow(points)-1

                ## calculate slope between first point and all others
                points[:,:slope] .= (points.Wᵢ[e] .- points.Wᵢ) ./ (points.τ[e] .- points.τ)
                cand = points.slope[e+1:nrow(points)]

                ## find the species with the steepest slope (must be negative)
                k = findall(points.slope .== minimum(cand[cand .< 0.0]))[1]
                if k == nrow(points) ## if k is the most drought tolerant species
                    break
                else
                    push!(coex, points[k,:spp]) ## add species to list of coexisting species
                    e = k ## start next iteration from this species
                end
            end
            push!(coex, points[nrow(points),:spp]) ## add latest feasible species to list (automatically coexists)

            ## Now determine whether the species which coexist on water alone close the canopy
            if ρ > calc_ρ_cc(T, sd[coex, :], F, μ, b)

                 ## loop through species from 1 to Q, determining whether they are feasible with closed
                 ## canopy
                 sd_sub = sd[infeas, :]

                 ## assign variables requires outside local scope
                 sd_sub2 = DataFrame()
                 hs = 0.0
                 τ_temp = Vector{Float64}(undef, nrow(sd_sub))

                 first_found = false
                 first_id = 1
                 while !first_found

                     ## calculate hs, check if species Q is feasible, if not, repeat until it is
                     cor_hs = false
                     last = nrow(sd_sub)
                     while !cor_hs

                         ## second subset (narrow this )
                         sd_sub2 = sd_sub[first_id:last,:]

                         #@infiltrate nrow(sd_sub2) == 0

                         ## calculate H*
                         hs = calc_hstar(sd_sub2, ρ, F, μ, T, uf, E, b)
                         if hs > 1e-15 ## some imprecision in hs calculation
                             ## check if latest species is infeasible

                             if nrow(sd_sub2) == 1
                                 cor_hs = true
                             else

                                 hs_up = calc_hstar(sd_sub2[1:nrow(sd_sub2)-1,:], ρ, F, μ, T, uf, E, b)

                                 τ_temp = Vector{Float64}(undef, nrow(sd_sub2))
                                 ## calculate τ for all other species
                                 for i in 1:nrow(sd_sub2)
                                     τ_temp[i] = calc_τ_cc(hs_up, T, sd_sub2[i,:Aₘ], sd_sub2[i, :r],
                                                           sd_sub2[i, :α], sd_sub2[i, :δ], μ, uf)
                                 end

                                 if τ_temp[nrow(sd_sub2)] >= T
                                     last -= 1
                                 else
                                     cor_hs = true
                                 end

                             end
                         else
                             first_id += 1
                             break
                         end
                     end

                     if cor_hs

                         τ_temp = Vector{Float64}(undef, nrow(sd_sub2))
                         ## calculate τ for all other species
                         for i in 1:nrow(sd_sub2)
                             τ_temp[i] = calc_τ_cc(hs, T, sd_sub2[i,:Aₘ], sd_sub2[i, :r],
                                                   sd_sub2[i, :α], sd_sub2[i, :δ], μ, uf)
                         end

                         W₀ = minimum([minimum(sd_sub2[:, :Wᵢ]) + ρ, θ_fc])

                         ## create dataframe of species characteristics (including w₀) for convex hull algorithm
                         points = DataFrame(τ = vcat(τ_temp, 0.0), Wᵢ = vcat(sd_sub2.Wᵢ, W₀),
                                            spp = vcat(sd_sub2.spp, 0))
                         sort!(points, :spp) ## sort in reverse drought tolerance order

                         ## calculate slope between first point and all others
                         points[:,:slope] .= (points.Wᵢ[1] .- points.Wᵢ) ./ (points.τ[1] .- points.τ)
                         cand = points.slope[2:nrow(points)]
                         k = findall(points.slope .== minimum(cand[cand .< 0.0]))[1]
                         if k == 2
                             first_found = true
                         else
                             first_id += 1
                         end
                     end
                 end

                 ## create dataframe of species characteristics (including w₀) for convex hull algorithm
                 points = DataFrame(τ = vcat(τ_temp, 0.0), Wᵢ = vcat(sd_sub2.Wᵢ, W₀),
                                    spp = vcat(sd_sub2.spp, 0))
                 sort!(points, :spp) ## sort in reverse drought tolerance order

                 coex = [] ## create empty object to populate with coexisting species indices
                 e = 1
                 ## loop over species starting with least drought tolerant
                 while e < nrow(points)-1

                     ## calculate slope between first point and all others
                     points[:,:slope] .= (points.Wᵢ[e] .- points.Wᵢ) ./ (points.τ[e] .- points.τ)
                     cand = points.slope[e+1:nrow(points)]

                     ## find the species with the steepest slope (must be negative)
                     k = findall(points.slope .== minimum(cand[cand .< 0.0]))[1]
                     if k == nrow(points) ## if k is the most drought tolerant species
                         break
                     else
                         push!(coex, points[k,:spp]) ## add species to list of coexisting species
                         e = k ## start next iteration from this species
                     end
                 end
                 push!(coex, points[nrow(points),:spp]) ## add latest feasible species to list (automatically coexists)

                 if nrow(sd_sub2) > 1

                     ## finally, check if late species can invade earliest feasible species

                     #@infiltrate
                     ntau = τ_temp[sd_sub2.spp .∈ [coex]]
                     sd_sub2 = sd[coex,:]
                     sd_sub2.τ .= ntau

                     ## calculate H*
                     hs = calc_hstar(DataFrame(sd_sub2[1,:]), ρ, F, μ, T, uf, E, b)
                     τ_temp2 = Vector{Float64}(undef, nrow(sd_sub2))
                     ## calculate τ for all other species
                     for i in 1:nrow(sd_sub2)
                         τ_temp2[i] = calc_τ_cc(hs, T, sd_sub2[i,:Aₘ], sd_sub2[i, :r],
                                                sd_sub2[i, :α], sd_sub2[i, :δ], μ, uf)
                     end

                     if any(τ_temp2 .> T)

                         hs = calc_hstar(DataFrame(sd_sub2[τ_temp2 .< T,:]), ρ, F, μ, T, uf, E, b)
                         ## need to recalculate τ AGAIN if more than one species incl
                         τ_temp2 = Vector{Float64}(undef, nrow(sd_sub2))
                         ## calculate τ for all other species
                         for i in 1:nrow(sd_sub2)
                             τ_temp2[i] = calc_τ_cc(hs, T, sd_sub2[i,:Aₘ], sd_sub2[i, :r],
                                                    sd_sub2[i, :α], sd_sub2[i, :δ], μ, uf)
                         end

                         return Union{Vector{Int}, Vector{Float64}, Bool, Bool}[sd_sub2[τ_temp2 .< T,:spp],
                                                                                τ_temp2[τ_temp2 .< T], true, true]
                     else
                         return Union{Vector{Int}, Vector{Float64}, Bool, Bool}[Int.(coex), sd_sub2[:, :τ], true, false]
                     end

                else
                     return Union{Vector{Int}, Vector{Float64}, Bool, Bool}[Int.(coex), sd_sub2[:, :τ], true, false]
                end

                ## if canopy never closes, return list of species which coexist on water alone
            else
                return Union{Vector{Int}, Vector{Float64}, Bool, Bool}[sd[coex,:spp], sd[coex, :τ], false, false]
            end
        end
    end

    println(coex)

end;

"""
    calc_Λ(ρ::Float64, w::Vector{Float64}, τ::Float64, E::Float64, θ_fc::Float64 = 0.6)

Calculates the equilibrium canopy area for each species in a community given the storm size `ρ`,
a vector of species' critical water contents `w`, a vector of species' break-even times `τ`, the transpiration
constant `E`, and the field capactity `θ_fc`.

Note that this function does not do a feasibility screen. The supplied vectors of species characteristics
must be pre-screened.
"""
function calc_Λ(ρ::Float64, w::Vector{Float64}, τ::Vector{Float64}, E::Float64, θ_fc::Float64 = 0.6)
    Λ = Vector{Float64}(undef, length(τ))
    w = vcat(minimum([θ_fc, w[length(w)]+ρ]), w)
    τ = vcat(0.0, τ)

    ## loop through all species except most drought tolerant species
    for i in 2:length(w)-1
        Λ[i-1] = (1/E) * ((w[i-1] - w[i]) / (τ[i] - τ[i-1]) - (w[i] - w[i+1]) / (τ[i+1] - τ[i]))
    end

    ## calculate Λ for most drought tolerant species
    Λ[length(w)-1] = (1 / E) * ((w[length(w)-1] - w[length(w)]) / (τ[length(τ)] - τ[length(τ)-1]))

    return Λ

end


"""
    calc_eqN(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04, F::Float64 = F, E::Float64 = E,
             θ_fc::Float64 = 0.4, μ::Float64 = 0.1, analytical::Bool = false, understory_transpiration::Bool = false,
             uf::Float64 = 0.1)

Calculate equilibrium population density for species in `spp_data`. Determines whether species close the canopy, and then
calculates resulting equilibria accordingly.

"""
function calc_eqN(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04, F::Float64 = F, E::Float64 = E,
                  θ_fc::Float64 = 0.4, μ::Float64 = 0.1, analytical::Bool = false, understory_transpiration::Bool = false,
                  uf::Float64 = 0.1, b::Float64 = 3.0)

    output = Vector{Int64}(undef, nrow(sd))
    eqN = Vector{Float64}(undef, nrow(sd))
    sort!(sd, :Wᵢ, rev = true) ## ensure sd is ordered correctly

    ## check feasibility of each species in sd
    fs = feas(sd, T, ρ, F, μ, θ_fc, uf, analytical, false, E, b)

    if isempty(fs[1]) ## if no species are feasible set each pop. density to 0
        sd[:,:eqN] .= 0.0
    else

        ## calculate starting soil water content, capped at θ_fc (field capacity)
        W₀ = minimum([θ_fc, sd[fs[1][length(fs[1])], :Wᵢ] + ρ])

        ## set equilibrium pop. density of all species not in fs (feasible spp) to 0
        eqN[(!in).(sd.spp, Ref(fs[1]))] .= 0.0

        Λ = calc_Λ(ρ, sd.Wᵢ[fs[1]], fs[2], E, θ_fc)
        eqN[fs[1]] .= F .* Λ .* T .* (exp(μ * T) / (exp(μ * T) - 1))
        sd[:,:eqN] .= eqN

    end

    ## sort output and return
    output = sort(sd, :spp)
    return output

end;

"""
    check_eq_agreement(iter::Int64, Nspp::Int64, Nyr::Int64 = 2000,
                            W0::Float64 = 0.6)

Checks agreement between equilibrium population density calculations and simulations.
"""
function check_eq_agreement(iter::Int64, Nspp::Int64, Nyr::Int64 = 300, θ_fc::Float64 = 0.4,
                            P::Float64 = 10.0, mean_p::Float64 = 0.4,
                            μ::Float64 = 0.15, F::Float64 = 10.0, uf::Float64 = 0.1, n_hts::Int64 = 1,
                            b::Float64 = 3.0)
    eq_sim = Vector{Float64}(undef, iter*Nspp)
    eq_an = Vector{Float64}(undef, iter*Nspp)
    mt = mortality_table(Int(Nyr * P), μ, repeat([1.0 / P], inner = Int(Nyr * P)))
    @Threads.threads for i in 1:iter

        no_pe = false
        spd = DataFrame()
        while !no_pe
            spd = generate_spp_data(Nspp, 0.7, n_hts, 1.0 / P, F, μ, b, 0.4, 0.0,
                                         0.06, 0.03, 1.5, 5.0, 0.5, 1.5)
            fs = feas(spd, 1.0 / P, mean_p / P, F, μ, θ_fc, uf, false, false, E, b)
            if !fs[4]
                no_pe = true
            end
        end

        eq_sim[[((i-1)*Nspp)+1:1:i*Nspp;]] =
            Vector(sim_water_ppa(spd, Nyr, Nspp, 1.0, μ, F, P, mean_p, θ_fc,
                    mt, false, 0.4, 3.0, uf, false)[2][Int(round(Nyr*P)),2:Nspp+1])
        eq_an[[((i-1)*Nspp)+1:1:i*Nspp;]] =
           calc_eqN(spd, 1.0 / P, mean_p / P, F, E, θ_fc, μ, false, false, uf, b)[:,:eqN]
    end
    return DataFrame(eq_sim = eq_sim, eq_an = eq_an)
end;


"""
    plot_eq_agreement(data::DataFrame, save::Bool = false, filename::String = "")

Generates a plot of calculated vs. simulated equilibrium population density.
"""
function plot_eq_agreement(data::DataFrame, save::Bool = false, filename::String = "")
   p = plot(framestyle = :box, grid = false,
            legend = :none, ylab = "Equilibrium density (simulated)", xlab = "Equilibrium density (predicted)",
            ylim = [0, maximum(data.eq_sim)+0.1*maximum(data.eq_sim)], xlim = [0, maximum(data.eq_an)+0.1*maximum(data.eq_an)])
   p = Plots.abline!(1,0)
   p = plot!(data.eq_an, data.eq_sim, seriestype = "scatter")
   if save
       savefig(filename)
   end
   return p
end



"""
    calc_eq_biomass(spp_data::DataFrame, ρ::Float64, E::Float64, T::Float64, F::Float64,
                    μ::Float64, uf::Float64, θ_fc::Float64 = 0.6)

Calculates equilibrium biomass of species in `spp_data.`
"""
function calc_eq_biomass(spp_data, ρ, E, T, F, μ, uf, θ_fc::Float64 = 0.6, analytical::Bool = false, b::Float64 = 3.0)

    ## check feasibility of each species in spp_data
    fs = feas(spp_data, T, ρ, F, μ, θ_fc, uf, analytical, false, E, b)

    ## create empty vector of eq. biomass values, to be populated later
    eqB = zeros(nrow(spp_data))

    if isempty(fs) ## if no species are feasible set each pop. density to 0
        return eqB
    else

        ## remove infeasible species from dataframe
        spp_data = spp_data[fs[1], :]

        ## If canopy is open
        if ρ < calc_ρ_cc(T, spp_data, F, μ, b)

            Λ = calc_Λ(ρ, spp_data.Wᵢ, spp_data.τ, E, θ_fc)

            r = F .* Λ

            function int(t)
                exp(-μ * t) * t^b * (1 / (F * T^b * sigma(μ*T, b-1))) ^ (b / (b-1))
            end

            eqB[fs[1]] .= r .* quadgk(int, 0.0, 1e3)[1]
            return eqB

        else

            ## calculate equilibrium canopy area of each species
            Λ = calc_Λ(ρ, spp_data.Wᵢ, fs[2], E, θ_fc)
            ## calculate initial density of each cohort at equilibrium (aka equilibrium reproduction)
            r = F .* Λ

            hs = calc_hstar(spp_data, ρ, F, μ, T, uf, E, 3.0)
            xs = (hs ./ spp_data.α) .^ 2

            gs = (1 ./ (spp_data.α .* spp_data.δ)) .* (spp_data.Aₘ .* fs[2] ./ T .- spp_data.r)

            for i in 1:nrow(spp_data)

                ## first integral sums individual biomass values of all cohorts with x < xstar (understory cohorts),
                ## discounted for mortality
                function int1(t)
                    exp(-μ * t) * (uf * gs[i] * t) ^ b
                end

                ## second integral sums individual biomass values of all cohorts in overstory, discounted for
                ## mortality
                function int2(t)
                    exp(-μ * t) * (xs[i] + gs[i] * (t - (xs[i] / (uf * gs[i])))) ^ b
                end

                ## multiply sum of the individual biomass integrals by the starting population density of each cohort, r
                eqB[fs[1][i]] = (r * (quadgk(int1, 0.0, (xs[i] / (uf * gs[i])))[1] +
                    quadgk(int2, (xs[i] / (uf * gs[i])), 1e3)[1]))[1]

            end

            return eqB
        end
    end
end



"""
    check_eq(bd, tol = 1e-15)

Checks whether system has reached equilibrium given a DataFrame of density dynamics, as
generated by `sim_water_ppa()` and a tolerance `tol.`
"""
function check_eq(bd, tol = 1e-15)
    check = (Array(bd[nrow(bd),[2:1:ncol(bd);]]) .-
        Array(bd[nrow(bd)-1,[2:1:ncol(bd);]])) ./
        Array(bd[nrow(bd)-1,[2:1:ncol(bd);]])
    if all(check[.!isnan.(check)] .< tol)
        return true
    else return false
    end
end;

##---------------------------------------------------------------
## Scratch, keeping for reference
##---------------------------------------------------------------


"""
    feas(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04,
    F::Float64 = 10.0, μ::Float64 = 0.1, θ_fc::Float64 = 0.4, uf::Float64 = 0.1,
    analytical::Bool = false, understory_transpiration::Bool = false)

Use convex hull algorithm to determine the feasibility of each species in `sd`. `sd` should be a `DataFrame`
as generated by `generate_spp_data()`. Returns a vector of indices giving the subset of species identifiers
in `sd` which are feasible given their break-even time and critical water content.
"""
function feas_old(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04,
              F::Float64 = 10.0, μ::Float64 = 0.1, θ_fc::Float64 = 0.4, uf::Float64 = 0.1,
              analytical::Bool = false, understory_transpiration::Bool = false)

    W₀ = minimum([minimum(sd[:, :Wᵢ]) + ρ, θ_fc])

    ## first test for feasibility in water-only scenario
    infeas = μ_feas(sd, F, μ, analytical)
    first_feas = minimum(findall(sd.Wᵢ .< W₀), init = 9999)
    if first_feas == 9999 ## case where no species are feasible
        return []
    else

        ## create dataframe of species characteristics (including w₀) for convex hull algorithm
        points = DataFrame(τ = vcat(sd.τ[infeas], 0.0), Wᵢ = vcat(sd.Wᵢ[infeas], W₀),
                           spp = vcat(sd.spp[infeas], 0))
        sort!(points, :Wᵢ, rev = true) ## sort in reverse drought tolerance order

        ## remove all infeasible species from points dataframe
        points = points[first_feas:nrow(points), :]

        if nrow(points) == 1 ## case where all species infeasible
            return []
        else

            coex = [] ## create empty object to populate with coexisting species indices
            e = 1
            ## loop over species starting with least drought tolerant
            while e < nrow(points)-1

                ## calculate slope between first point and all others
                points.slope .= (points.Wᵢ[e] .- points.Wᵢ) ./ (points.τ[e] .- points.τ)
                cand = points.slope[e+1:nrow(points)]

                ## find the species with the steepest slope (must be negative)
                k = findall(points.slope .== minimum(cand[cand .< 0.0]))[1]
                if k == nrow(points) ## if k is the most drought tolerant species
                    break
                else
                    push!(coex, points[k,:spp]) ## add species to list of coexisting species
                    e = k ## start next iteration from this species
                end
            end
            push!(coex, points[nrow(points),:spp]) ## add latest feasible species to list (automatically coexists)

            ## Now determine whether the species which coexist on water alone close the canopy
            if ρ > calc_ρ_cc(T, sd[coex, :], F, μ, b)

                ## determine best height strategy and eliminate all others
                ht_list = unique(sd[coex, :ht])
                xss_list = Vector{Float64}(undef, length(ht_list))
                sd_sub = sd[coex, :]
                ## loop through height strategies and calculate equilibrium canopy closure size
                for i in 1:length(ht_list)
                    sd_sub_sub = sd_sub[sd_sub.ht .== ht_list[i], :]
                    xss_list[i] = calc_xss(sd_sub_sub[sd_sub_sub.Wᵢ .== maximum(sd_sub_sub.Wᵢ), :C₁][1],
                                           sd_sub_sub[sd_sub_sub.Wᵢ .== maximum(sd_sub_sub.Wᵢ), :C₂][1],
                                           E, T, ρ, F, μ, uf, understory_transpiration)

                end

                if all(xss_list .== 0.0) ## indicates no strategy closes canopy
                else

                    ## identify the best height strategy
                    ht_id = findall(xss_list .^ ht_list .== maximum(xss_list .^ ht_list))

                    ## thin list of coexisting species
                    coex = coex[findall(sd[coex, :ht] .== ht_list[ht_id])]

                end

                ## next check if capped transpiration limits growth period of late species
                if length(coex) == 1 ## if only one species remains in coex list, we are good
                    return coex[sd[coex, :Wᵢ] .< W₀] ## provided that species' Wᵢ is lower than W₀
                ## if expression true, latest species is capped
                elseif ρ > calc_ρ_cx(sd[coex, :], E, T)
                    ## get minimum critical water content from coexisting species
                    Wᵢ_min = sd[coex, :Wᵢ][length(coex)]

                    if length(coex) > 2

                        ## loop through species, eliminating those which are limited by capped transpiration
                        for i in [length(coex)-1:-1:2;]
                            if ρ > calc_ρ_cx(sd[coex[1:i], :], E, T)
                                Wᵢ_min = sd[coex, :Wᵢ][i]
                            else
                                break
                            end
                        end

                    end

                    return coex[sd[coex,:Wᵢ] .< W₀ .&& sd[coex, :Wᵢ] .> Wᵢ_min]
                ## otherwise return current list
                else
                    return coex[sd[coex,:Wᵢ] .< W₀]
                end

            ## if canopy never closes, return list of species which coexist on water alone
            else
                return coex[sd[coex,:Wᵢ] .< W₀]
            end

        end

    end

    println(coex)

end;


"""
    calc_eqN(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04, F::Float64 = F, E::Float64 = E,
             θ_fc::Float64 = 0.4, μ::Float64 = 0.1, analytical::Bool = false, understory_transpiration::Bool = false,
             uf::Float64 = 0.1)

Calculate equilibrium population density for species in `spp_data`. Determines whether species close the canopy, and then
calculates resulting equilibria accordingly.

"""
function calc_eqN_old(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04, F::Float64 = F, E::Float64 = E,
                  θ_fc::Float64 = 0.4, μ::Float64 = 0.1, analytical::Bool = false, understory_transpiration::Bool = false,
                  uf::Float64 = 0.1)

    output = Vector{Int64}(undef, nrow(sd))
    eqN = Vector{Float64}(undef, nrow(sd))
    sort!(sd, :Wᵢ, rev = true) ## ensure sd is ordered correctly

    ## check feasibility of each species in sd
    fs = feas(sd, T, ρ, F, μ, θ_fc, uf, analytical, understory_transpiration)

    if isempty(fs) ## if no species are feasible set each pop. density to 0
        sd.eqN .= 0.0
    else

        ## calculate starting soil water content, capped at θ_fc (field capacity)
        W₀ = minimum([θ_fc, sd[fs[length(fs)], :Wᵢ] + ρ])

        ## set equilibrium pop. density of all species not in fs (feasible spp) to 0
        eqN[(!in).(sd.spp, Ref(fs))] .= 0.0

        ## check if canopy closed
        if ρ > calc_ρ_cc(T, sd[fs,:], F, μ, b)

            ## calculate equilibrium population density for closed canopy system
            eqN[fs] .= calc_eqN_cc(sd[fs, :], ρ, E, T, F, μ, θ_fc)

            ## and set to output
            sd.eqN = eqN

        ## if canopy is open
        else

            ## if analytical is set to true, treat system as continuously reproducting
            if analytical

                ## calculate constant that will be reused in subsequent expressions
                cnst = F * T / E

                ## if only one species has a feasible equilibrium, calculate equilibrium reproduction
                if length(fs) == 1
                    eqN[fs[1]] = cnst * ((W₀ - sd[fs[1], :Wᵢ]) / (sd[fs[1], :τ]))
                else
                    ## loop through feasible species and calculate equilibrium reproduction, expression depends on relative phenology
                    for s in [1:1:length(fs);]
                        if s == 1
                            eqN[fs[s]] = cnst * ((W₀ - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ]) -
                                                 (sd[fs[s], :Wᵢ] - sd[fs[s+1], :Wᵢ]) / (sd[fs[s+1], :τ] - sd[fs[s], :τ]))
                        elseif s == length(fs)
                            eqN[fs[s]] = cnst * ((sd[fs[s-1], :Wᵢ] - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ] - sd[fs[s-1], :τ]))
                        else
                            eqN[fs[s]] = cnst * (((sd[fs[s-1], :Wᵢ] - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ] - sd[fs[s-1], :τ])) -
                                                 ((sd[fs[s], :Wᵢ] - sd[fs[s+1], :Wᵢ]) / (sd[fs[s+1], :τ] - sd[fs[s], :τ])))
                        end
                    end
                end

            ## if analytical is set to false, treat system as if reproduction occurs once per storm interval
            else

                ## calculate constant that will be reused in subsequent expressions
                cnst = F * T / E

                if length(fs) == 1
                    eqN[fs[1]] = cnst * ((W₀ - sd[fs[1], :Wᵢ]) / (sd[fs[1], :τ]))
                else
                    for s in [1:1:length(fs);]
                        if s == 1
                            eqN[fs[s]] = cnst * ((W₀ - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ]) -
                                (sd[fs[s], :Wᵢ] - sd[fs[s+1], :Wᵢ]) / (sd[fs[s+1], :τ] - sd[fs[s], :τ]))
                        elseif s == length(fs)
                            eqN[fs[s]] = cnst * ((sd[fs[s-1], :Wᵢ] - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ] - sd[fs[s-1], :τ]))
                        else
                            eqN[fs[s]] = cnst * (((sd[fs[s-1], :Wᵢ] - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ] - sd[fs[s-1], :τ])) -
                                    ((sd[fs[s], :Wᵢ] - sd[fs[s+1], :Wᵢ]) / (sd[fs[s+1], :τ] - sd[fs[s], :τ])))
                        end
                    end
                end
            end
            ## multiply equilibrium density by following constant to account for mortality
            sd.eqN = eqN .* (exp(μ * T) / (exp(μ * T) - 1))
        end
    end

    ## sort output and return
    output = sort(sd, :spp)
    return output

end;

"""
    calc_Λ(ρ::Float64, w::Vector{Float64}, τ::Float64, E::Float64, θ_fc::Float64 = 0.6)

Calculates the equilibrium canopy area for each species in a community given the storm size `ρ`,
a vector of species' critical water contents `w`, a vector of species' break-even times `τ`, the transpiration
constant `E`, and the field capactity `θ_fc`.

Note that this function does not do a feasibility screen. The supplied vectors of species characteristics
must be pre-screened.
"""
function calc_Λ_old(ρ::Float64, w::Vector{Float64}, τ::Vector{Float64}, E::Float64, θ_fc::Float64 = 0.6)
    Λ = Vector{Float64}(undef, length(τ))
    w = vcat(minimum([θ_fc, w[length(w)]+ρ]), w)
    τ = vcat(0.0, τ)

    ## loop through all species except most drought tolerant species
    for i in 2:length(w)-1
        Λ[i-1] = (1/E) * ((w[i-1] - w[i]) / (τ[i] - τ[i-1]) - (w[i] - w[i+1]) / (τ[i+1] - τ[i]))
    end

    ## calculate Λ for most drought tolerant species
    Λ[length(w)-1] = (1 / E) * ((w[length(w)-1] - w[length(w)]) / (τ[length(τ)] - τ[length(τ)-1]))

    return Λ

end



"""
    calc_eqN_cc(spp_data::DataFrame, ρ::Float64, E::Float64,
                T::Float64, F::Float64, μ::Float64, θ_fc::Float64 = 0.6)

Calculates equilibrium abundance of species in `spp_data` provided those species close the canopy.
This function does not itself check whether the species close the canopy, rather it is intended as a
helper function for the more general function `calc_eqN()`.

"""
function calc_eqN_cc_old(spp_data, ρ, E, T, F, μ, θ_fc::Float64 = 0.6)
    gstar = calc_gstar(spp_data, ρ, E, T, θ_fc)
    τ = Vector{Float64}(undef, nrow(spp_data))
    for i in 1:nrow(spp_data)
        τ[i] = calc_τ_cc(gstar, T, spp_data.C₁[i], spp_data.C₂[1])
    end
    Λ = calc_Λ(ρ, spp_data.Wᵢ, τ, E, θ_fc)
    return F .* Λ .* T .* (exp(μ * T) / (exp(μ * T) - 1))
end
