using CSV
using Test
using DataFrames
using Random
using Distributions
using LinearAlgebra

mutable struct metapop
    n_pops::Int64
    n_alleles::Int64
    migration_rate::Float64
    eff_pop_size::Int64
    ct_map::Array{Float64}
    diskernel::Array{Float64}
    metapop(n_pops::Int64, n_alleles::Int64, migration_rate::Float64, eff_pop_size::Int64) = new(n_pops, n_alleles, migration_rate, eff_pop_size, zeros(Float64, n_pops, n_alleles), zeros(Float64, n_pops, n_pops))
end

function get_random_pmf(n_categories::Int64)
    r_interval::Array{Float64} = rand(Uniform(), n_categories)
    s::Float64 = sum(r_interval)
    for j = 1:n_categories
        r_interval[j] = r_interval[j] / s
    end
    return(r_interval)
end

function init_random_mp(mp::metapop)
    n_p::Int64 = mp.n_pops
    n_al::Int64 = mp.n_alleles
    for i = 1:n_p
        r_interval::Array{Float64} = get_random_pmf(n_al)
        for j = 1:n_al
            mp.ct_map[i,j] = r_interval[j]
        end
    end
end

function init_uniform_mp(mp::metapop)
    n_p::Int64 = mp.n_pops
    n_al::Int64 = mp.n_alleles
    eff_pop_size::Int64 = mp.eff_pop_size


    p = fill(1.0/n_al, n_al)


    for i = 1:n_p
        x = rand(Multinomial(eff_pop_size, p))
        for j = 1:n_al
            mp.ct_map[i,j] = x[j]
        end
    end
end


# get the new proportion of alleles in pop j based on drift
function get_new_pj(x_j::Array{Float64}, n_al::Int64, eff_pop_size::Int64)
    normalize!(x_j, 1)
    xj_new = rand(Multinomial(eff_pop_size, x_j))
    return(xj_new)
end

function draw_from_diskern_row(diskernel_row::Array{Float64})
    u = rand(Uniform())

    ind::Int64 = 1
    s::Float64 = 0.0
    while (s < u)
        s += diskernel_row[ind]
        if s > u
            return ind
        end
        ind += 1
    end

    return ind
end

function get_allele_list(pop_alleles::Array{Float64}, eff_pop_size::Int64)
    allele_list::Array{Int64} = zeros(eff_pop_size)

    n_alleles::Int64 = length(pop_alleles)

    ct::Int64 = 1
    for i = 1:n_alleles
        new = pop_alleles[i]
        for j = 1:new
            allele_list[ct] = i
            ct += 1
        end
    end
    return(allele_list)
end

function migration(pop_from::Int64, pop_from_cts::Array{Float64}, n_indivs_leaving::Int64, mp::metapop)
    diskernel_row::Array{Float64} = mp.diskernel[pop_from,:]
    n_pops::Int64 = mp.n_pops
    n_alleles::Int64 = mp.n_alleles
    n_haplo::Int64 = 2
    eff_pop_size::Int64 = sum(pop_from_cts)

    if (eff_pop_size/2 < n_indivs_leaving)
        n_indivs_leaving = eff_pop_size/2
    end


    for ind = 1:n_indivs_leaving
        pop_to = draw_from_diskern_row(diskernel_row)
        @assert pop_to != pop_from
        for h =1:n_haplo
            eff_pop_size = sum(pop_from_cts)
            allele_list::Array{Int64} = get_allele_list(pop_from_cts, eff_pop_size)
            random_index = rand(DiscreteUniform(1, eff_pop_size))

            allele_num = allele_list[random_index]

            pop_from_cts[allele_num] -= 1
            mp.ct_map[pop_from, allele_num] -= 1
            mp.ct_map[pop_to, allele_num] += 1

            @assert mp.ct_map[pop_from, allele_num] >= 0
            @assert mp.ct_map[pop_to, allele_num] >= 0
        end
    end

    return(mp)
end

function run_gen(mp::metapop)
    n_p::Int64 = mp.n_pops
    n_al::Int64 = mp.n_alleles
    migration_rate::Float64 = mp.migration_rate

    # drift

    post_drift::Array{Float64} = zeros(n_p, n_al)

    # store this in a new array
    for p = 1:n_p
        eff_pop_size::Int64 = sum(mp.ct_map[p,:])
        new_pj = get_new_pj(mp.ct_map[p, :], n_al, eff_pop_size)
        post_drift[p,:] = new_pj
    end

    #post_drift_now = deepcopy(post_drift)
    mp.ct_map = post_drift
    # migration

    if (false)
        for p = 1:n_p
            post_drift_this_pop = (post_drift_now[p,:])
            eff_pop_size::Int64 = sum(post_drift_this_pop)
            exp_n_alleles_leaving::Float64 = (eff_pop_size/2.0) * migration_rate
            rem::Float64 = exp_n_alleles_leaving - floor(exp_n_alleles_leaving)
            extra_mig::Int64 = 0
            if (rand(Uniform()) < rem)
                extra_mig = 1
            end
            n_indivs_leaving::Int64 = (floor(exp_n_alleles_leaving) + extra_mig)

            # pick random in 1...eff_pop, and that in the ct is the allele
            # update the real df, pass it the post drift matrix to compute where migrants go,


            migration(p, post_drift[p,:], n_indivs_leaving, mp)
        end
    end
end

function update_df(df::DataFrame, gen::Int64, state::Array{Float64}, eff_pop_size::Int64)
    xdim = size(state)[1]
    ydim = size(state)[2]

    for p = 1:xdim
        for al = 1:ydim
            push!(df.gen, gen)
            push!(df.pop, p)
            push!(df.locus, al)
            freq::Float64 = state[p,al]/eff_pop_size
            push!(df.freq, freq)
        end
    end
end

function update_df(df::DataFrame, gen::Int64, jost_d::Float64, gst::Float64)
    push!(df.gen, gen)
    push!(df.jostd, jost_d)
    push!(df.gst, gst)
end

function init_uniform_dispersal_kernel(mp::metapop)
    n_pops::Int64 = mp.n_pops
    uniform_prob::Float64 = 1.0 / (n_pops-1)

    for i = 1:n_pops
        for j = 1:n_pops
            if (i != j)
                mp.diskernel[i,j] = uniform_prob
            else
                mp.diskernel[i,j] = 0.0
            end
        end
    end
end

function calc_ht(state, pops, n_pops, n_alleles)
    outer_sum::Float64 = 0.0
    for al = 1:n_alleles
        inner_s::Float64 = 0.0
        for p = 1:n_pops
            eff_pop_size::Float64 = sum(state[p,:])

            inner_s += (state[p,al] / eff_pop_size)
        end
        inner_s = inner_s / n_pops

        outer_sum += (inner_s^2)
    end
    return(1.0 - outer_sum)
end

function calc_hs(state, pops, n_pops, n_alleles)
    s::Float64 = 0.0
    for p = 1:n_pops
        eff_pop_size::Float64 = sum(state[p,:])
        inner_s::Float64 = 0.0
        for al = 1:n_alleles
            inner_s += (state[p,al] / eff_pop_size)^2
        end


        s += (1.0 - inner_s)
    end
    s = s / n_pops
    return s
end


function calc_jost_d(state::Array{Float64}, pops::Array{Int64},  n_pops, n_alleles)
    H_T::Float64 = calc_ht(state, pops, n_pops, n_alleles)
    H_S::Float64 = calc_hs(state, pops, n_pops, n_alleles)

    jostD::Float64 = ((H_T - H_S)*n_pops) / ((1.0 - H_S)*(n_pops-1))
    return(jostD)
end

function calc_gst(state::Array{Float64}, pops::Array{Int64},  n_pops, n_alleles)
    H_T::Float64 = calc_ht(state, pops, n_pops, n_alleles)
    H_S::Float64 = calc_hs(state, pops, n_pops, n_alleles)
    gst::Float64 = (H_T - H_S) / (H_T)
    return(gst)
end


function run_dke(n_gen::Int64, n_pops::Int64, n_alleles::Int64, migration_rate::Float64, eff_pop_size::Int64, log_freq::Int64)
    mp::metapop = metapop(n_pops, n_alleles, migration_rate, eff_pop_size)

    #init_random_mp(mp)
    init_uniform_mp(mp)

    init_uniform_dispersal_kernel(mp)

    #df = DataFrame()
    #df.gen = []
    #df.pop = []
    #df.locus = []
    #df.freq = []

    df = DataFrame()
    df.jostd = []
    df.gst = []
    df.gen = []

    all_pops = collect(1:n_pops)
    for g = 0:n_gen
        if g % log_freq == 0
            state::Array{Float64} = mp.ct_map
            #update_df(df, g, state, eff_pop_size)
            jostd::Float64 = calc_jost_d(state, all_pops,  n_pops, n_alleles)
            gst::Float64 = calc_gst(state, all_pops, n_pops, n_alleles)
            update_df(df,g,jostd,gst)
        end
        run_gen(mp)
    end

    return df
end

#plot(df())
