using Printf
using TOML
using PyPlot

""" Performs a pruning analysis

data is plotted for different seeds and user enters start and end 
from which to take data. The result is written to a TOML file 
"""
function prune_analysis(data, file::AbstractString)
    # Plot timeseries
    fig, ax = subplots()
    for (key, val) in data
        ax.plot(1:length(val), val, label=key)
    end
    legend()
    show()

    # User input to decide pruning
    start_str = ""
    end_str = ""
    for (key, val) in data
        println("@ Seed " * key)
        
        print("  Start (default Inf): ")
        rd = readline()
        start = Inf
        if rd != ""
            start = rd
        end 
        print("  End   (default Inf): ")
        rd = readline()
        endd = Inf
        if rd != ""
            endd = rd
        end

        start_str *= @sprintf("%s = \"%s\"\n", key, string(start))
        end_str *= @sprintf("%s = \"%s\"\n", key, string(endd))
    end

    directory = dirname(file)
    mkpath(directory)
    open(file, "w") do fl
        write(fl, "[start]\n")
        write(fl, start_str * "\n")
        write(fl, "[end]\n")
        write(fl, end_str)
    end
    println("Pruning written to:")
    println(file)
end


function prune(data, prune_file::AbstractString)
    prune_dict = TOML.parsefile(prune_file)

    data_pruned = Dict()
    for (key, val) in data
        start_str = prune_dict["start"][string(key)]
        end_str = prune_dict["end"][string(key)]

        start = 0
        if start_str == "Inf"
            continue
        else
            start = parse(Int64, start_str)
        end

        dim = length(size(val))
        if dim == 2
            if end_str == "Inf"
                data_pruned[key] = val[start:end,:]
            else
                endd = parse(Int64, end_str)
                data_pruned[key] = val[start:endd,:]
            end
        elseif dim == 3
            if end_str == "Inf"
                data_pruned[key] = val[start:end,:,:]
            else
                endd = parse(Int64, end_str)
                data_pruned[key] = val[start:endd,:,:]
            end
        else
            error("unexpected data dimension")
        end            
    end
    return sort(collect(data_pruned), by = x->x[1])
end
