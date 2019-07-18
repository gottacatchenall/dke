include("./src/dke.jl")

function run()

    n_als = (5, 10, 30)
    n_pops = (10, 30, 100)
    m = (0.001, 0.01, 0.1)
    eff_pop_size_list = (30,100,300,1000)
    n_gen = 300
    n_rep = 50

    metadata = DataFrame()
    metadata.runid = []
    metadata.m = []
    metadata.eff_pop = []
    metadata.n_pops = []
    metadata.n_alleles = []

    df = DataFrame()
    df.gst = []
    df.jostd = []
    df.gen = []
    df.id = []
    CSV.write("output.csv", df)

    idct::Int64 = 0
    for base_mig in m
        for n_al in n_als
            for n_pop in n_pops
                for eff_pop in eff_pop_size_list
                    for rep = 1:n_rep
                        df = run_dke(n_gen, n_pop, n_al, base_mig, eff_pop, 10)
                        df.ids = repeat([idct], nrow(df))

                        push!(metadata.runid, idct)
                        push!(metadata.m, base_mig)
                        push!(metadata.eff_pop, eff_pop)
                        push!(metadata.n_pops, n_pop)
                        push!(metadata.n_alleles, n_al)

                        idct += 1
                        CSV.write("output.csv", df, append=true)
                    end
                end
            end
        end
    end
    CSV.write("metadata.csv", metadata)
end


@time run()
#@time df = run_dke(n_gen, n_pops, n_al, m, eff_pop_size, 10)
#CSV.write("output.csv", df)
