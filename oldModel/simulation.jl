using Distributions
using Plots
#using RecipesBase

Δt = 1/(24*60)

struct Herd
    count::Integer
    loc::Vector{Float64}
    v::Float64
    B::Float64
    cooldown::Float64
    t::Float64
end
Base.broadcastable(x::Herd) = Ref(x)

population(herds::Vector{Herd}) = sum(getfield.(herds,:count))

struct Predator
    loc::Vector{Float64}
    v::Float64
    B::Float64
    range::Float64
    cooldown::Float64
    t::Float64
    active::Bool
end
Base.broadcastable(x::Predator) = Ref(x)

veclen(x) = sqrt(sum(x .^ 2))

function normalize(x::Vector{Float64};to=1) 
    return to*x/veclen(x)
end

function moveherd(herd::Herd, preds::Vector{Predator}; herdvariance = 5, world = [10,10])
    newB = herd.B + randn()*sqrt(herdvariance*Δt)
    newCount = herd.count
    newLoc = herd.loc
    newt = herd.t - Δt
    #activepreds = [p for p in preds if p.active]
    #preddist = veclen.([herd.loc - p.loc for p in activepreds])
    preddist = veclen.([herd.loc - p.loc for p in preds])
    if length(preddist) != 0 && minimum(preddist) < .01 && newt < herd.cooldown
        # hunted
        newt = herd.cooldown + Δt*10
        newCount -= 1
    end
    newLoc = herd.loc + [cos(herd.B),sin(herd.B)]*herd.v*Δt
    if newt > 0
        newLoc = newLoc = herd.loc + 5*[cos(herd.B),sin(herd.B)]*herd.v*Δt
    end
    if abs(newLoc[1]) > world[1]/2 || abs(newLoc[1]) < 1
        newB = pi - mod(newB, 2*pi)
    elseif abs(newLoc[2]) > world[2]/2 || abs(newLoc[2]) < 1
        newB = 2*pi - mod(newB, 2*pi)
    end
    if any(abs.(newLoc) .> world/2 .+ .1)
        newLoc = rand(2) .* world - world/2
        newt = rand()*1/2
    end
    return Herd(newCount, newLoc, herd.v, newB, herd.cooldown, newt)
end

function movepredator(pred::Predator, herds::Vector{Herd}, preds::Vector{Predator}; predvariance = 1, comm = true, world = [10,10])
    newB = pred.B + randn()*sqrt(predvariance*Δt)
    # if veclen(pred.loc) < .01
    #     newB = rand()*2*pi
    # end
    newLoc = pred.loc + [cos(pred.B),sin(pred.B)]*pred.v*Δt
    if abs(newLoc[1]) > world[1]/2
        newB = pi - mod(newB, 2*pi)
    elseif abs(newLoc[2]) > world[2]/2
        newB = 2*pi - mod(newB, 2*pi)
    end
    newt = pred.t - Δt
    newactivity = false
    if newt > 0
        # going h/ waiting at (0,0)
        newB = atan(reverse(-newLoc)...)
    else
        # hunting
        herddist = veclen.([pred.loc - pos for pos in getfield.(herds,:loc)])
        activepreds = [p for p in preds if p.active && p.loc != pred.loc]
        preddist = veclen.([pred.loc - p.loc for p in activepreds])
        #println("c")
        if minimum(herddist) < pred.range
            # priority 1: go to herd in range
            closest = herds[argmin(herddist)]
            dir = closest.loc - pred.loc
            newactivity = true
            if veclen(dir) < .015
                # successful hunt
                newt = pred.cooldown
                newactivity = false
            end
            newB = atan(reverse(dir)...)
        elseif comm && length(preddist) != 0 && minimum(preddist) < pred.range
            # priority 2: go to active hunter
            dir = activepreds[argmin(preddist)].loc - pred.loc
            newB = atan(reverse(dir)...)
            #newactivity = true
        end
    end
    if any(abs.(newLoc) .> world/2 .+ .1)
        newLoc = [0,0]
        newt = rand()*1/2
    end
    return Predator(newLoc, pred.v, newB, pred.range, pred.cooldown, newt, newactivity)
end

function removeempty(herds::Vector{Herd})
    return [herd for herd in herds if herd.count>0]
end
    


function makeherd(herdmean, herdgen, world, speed, cooldown)
    n = herdmean
    if lowercase(herdgen) == "poisson"
        n = rand(X, 1)
    end
    loc = rand(2) .* world - world/2
    B = rand()*2*pi
    return Herd(n, loc, speed, B, cooldown, -1)
end

function makepredator(speed, range, cooldown)
    B = rand()*2*pi
    return Predator([0,0], speed, B, range, cooldown, rand()*1/2, false)
end

function initSim(nherds, npreds; herdmean=40, world = [10,10], herdspeed = 3, herdcool = 1/24, predspeed = 5, predrange = 1, predcool=16/24, herdgen="exact")
    herds = [makeherd(herdmean, herdgen, world, herdspeed, herdcool) for _ in 1:nherds]
    preds = [makepredator(predspeed, predrange, predcool) for _ in 1:npreds]
    return herds, preds
end

function sim(nherds, npreds; herdmean=40, world = [10,10],
    herdspeed = 3, herdcool = 1/24, predspeed = 5, predrange = 1, predcool=16/24,
    herdgen="exact", comm = true, thresh=1/4)
    # initialize herds
    herds, preds = initSim(nherds, npreds; herdmean, world, herdspeed, herdcool, predspeed, predrange, predcool, herdgen)
    threshhold = population(herds) * thresh
    steps = 0
    while population(herds) > threshhold
        herds = removeempty(moveherd.(herds, Ref(preds)))
        preds = movepredator.(preds, Ref(herds), Ref(preds); comm = comm)
        steps += 1
    end
    return steps * Δt
end

function plotsim(herds::Vector{Herd}, preds::Vector{Predator}; world::Vector=[10,10])
    t = eachrow(reduce(hcat,getfield.(herds,:loc)))
    scatter(t[1], t[2],
    xlims=(-world[1]/2,world[1]/2),
    ylims=(-world[2]/2,world[2]/2),
    label=nothing)
    t = eachrow(reduce(hcat,getfield.(preds,:loc)))
    scatter!(t[1], t[2], label=nothing)
end