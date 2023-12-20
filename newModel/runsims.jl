using JLD2, ProgressBars
include("simulations.jl")

N = 15
preycount = 250:250
herdmean = 191:5:246
npreds = [75]
comm = [false]
f = "simDict1.jld2"
D = Dict()
if isfile(f)
    D = load_object(f)
    println("Loaded past dict with ", length(D) ," keys")
end

for pc ∈ preycount, hm ∈ herdmean, np ∈ npreds, c in comm
    if haskey(D, (N, pc, np, hm, c))
        continue
    end
    println("running:\t",(N, pc, np, hm, c))
    sims = []
    
    for _ in ProgressBar(1:N)
        push!(sims, sim(pc/hm, np;herdmean = hm, comm=c))
    end
    D[(N, pc, np, hm, c)] = sims
    save_object(f, D)
end