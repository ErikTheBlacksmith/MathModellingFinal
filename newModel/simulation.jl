using Distributions
using Plots

Δt = 1/(24*60)


mutable struct Herd
    count::Integer
    loc::Vector{Float64}
    v::Float64
    B::Float64
    cooldown::Float64
    t::Float64
end

population(herds::Vector{Herd}) = sum(getfield.(herds,:count))

mutable struct Predator
    loc::Vector{Float64}
    v::Float64
    B::Float64
    range::Float64
    cooldown::Float64
    t::Float64
    pastloc::Union{Vector{Float64}, Nothing}
    state::Int8 # -1: going home. 0: home. 1: hunting. 2: chasing.
    following::Union{Herd, Nothing}
end

population(herds::Vector{Herd}) = sum(getfield.(herds,:count))

dist(x) = sqrt(sum(x .^2))
dist(x,y) = sqrt(sum((x-y) .^2))

function steptime!(herds::Vector{Herd}, preds::Vector{Predator};
    herdvariance = 5,
    predvariance = 1, predcomm = true,
    world = [10,10])

    # move herds
    stepherd!.(herds, Ref(herds), Ref(preds);
    herdvariance = herdvariance, world = world)

    # step predators
    steppred!.(preds, Ref(herds), Ref(preds);
    predvariance = predvariance, comm = predcomm,
    world = world)

end

function stepherd!(herd::Herd, herds::Vector{Herd}, preds::Vector{Predator};
    herdvariance = 5,
    world = [10,10])
    herd.B += randn()*sqrt(herdvariance*Δt)

    herd.t -= Δt
    if herd.t > 0
        herd.loc += 5*[cos(herd.B),sin(herd.B)]*herd.v*Δt
    else
        herd.loc += [cos(herd.B),sin(herd.B)]*herd.v*Δt
    end

    if abs(herd.loc[1]) > world[1]/2 || abs(herd.loc[1]) < 1
        herd.B = pi - mod(herd.B, 2*pi)
    elseif abs(herd.loc[2]) > world[2]/2 || abs(herd.loc[2]) < 1
        herd.B = 2*pi - mod(herd.B, 2*pi)
    end

    if any(abs.(herd.loc) .> world/2 .+ .1)
        herd.loc = rand(2) .* world - world/2
    end
end

function steppred!(pred, herds::Vector{Herd}, preds::Vector{Predator};
    predvariance = 1, comm = true,
    world = [10,10])

    if pred.state == -1
        # going home
        pred.following = nothing
        pred.B = atan(reverse(-pred.loc)...)
        pred.loc += [cos(pred.B),sin(pred.B)]*pred.v*Δt
        if dist(pred.loc) < .01
            pred.state = 0
        end
    elseif pred.state == 0
        # at home
        pred.following = nothing
        if !isnothing(pred.pastloc)
            pred.B = atan(reverse(pred.pastloc)...)
        end
        pred.t -= Δt
        if pred.t <= 0
            pred.state = 1
        end
    else
        # hunting
        herddist = dist.([pred.loc - pos for pos in getfield.(herds,:loc)])
        activepreds = [p for p in preds if p.state ==2 && p.loc != pred.loc]
        preddist = dist.([pred.loc - p.loc for p in activepreds])
        pred.state = 1
        if minimum(herddist) < pred.range
            # priority 1: go to herd in range
            closest = herds[argmin(herddist)]
            pred.following = closest
            dir = closest.loc - pred.loc
            pred.B = atan(reverse(dir)...)

            pred.state = 2
            if dist(dir) < .01
                # successful hunt
                pred.pastloc = pred.loc
                nearactivepreds = activepreds[preddist .< .1]
                push!(nearactivepreds, pred)

                closestherd = herds[argmin(herddist)]

                for np in nearactivepreds
                    closestherd.count -= 1
                    np.state = -1
                    np.t = np.cooldown
                end
                closestherd.t = closestherd.cooldown

            end
        elseif !isnothing(pred.following) && pred.following.count > 0
            pred.state = 2
            dir = pred.following.loc - pred.loc
            pred.B = atan(reverse(dir)...)
        elseif comm && length(preddist) > 0 && minimum(preddist) < pred.range
            # priority 2: go to active hunter
            #println("q")
            pred.following = preds[argmin(preddist)].following
        elseif pred.state == 2
            # actually not hunting lol 
            pred.state =1
            pred.following = nothing
        end

        pred.B += randn()*sqrt(predvariance*Δt)
        pred.loc += [cos(pred.B),sin(pred.B)]*pred.v*Δt
        if abs(pred.loc[1]) > world[1]/2
            pred.B = pi - mod(pred.B, 2*pi)
        elseif abs(pred.loc[2]) > world[2]/2
            pred.B = 2*pi - mod(pred.B, 2*pi)
        end
    end

    if any(abs.(pred.loc) .> world/2 .+ .1)
        pred.loc = [0,0]
    end
end

function makeherd(herdmean, herdgen, world, speed, cooldown)
    n = trunc(Int, herdmean)
    if rand() < n - herdmean
        n += 1
    end
    if lowercase(herdgen) == "poisson"
        n = rand(X, 1)
    end
    loc = rand(2) .* world - world/2
    B = rand()*2*pi
    return Herd(n, loc, speed, B, cooldown, -1)
end

function makepredator(speed, range, cooldown)
    B = rand()*2*pi
    return Predator([0,0], speed, B, range, cooldown, rand()*1/2, nothing, 0, nothing)
end

function removeempty(herds::Vector{Herd})
    return [herd for herd in herds if herd.count>0]
end

function initSim(nherds, npreds; herdmean=40, world = [10,10], herdspeed = 3, herdcool = 1/24, predspeed = 5, predrange = 1, predcool=16/24, herdgen="exact")
    herds = [makeherd(herdmean, herdgen, world, herdspeed, herdcool) for _ in 1:nherds]
    if mod(nherds,1) != 0
        push!(herds, makeherd(herdmean*mod(nherds,1), herdgen, world, herdspeed, herdcool))
    end
    preds = [makepredator(predspeed, predrange, predcool) for _ in 1:npreds]
    return herds, preds
end

function plotsim(herds::Vector{Herd}, preds::Vector{Predator}; world::Vector=[10,10], diff = false)
    t = eachrow(reduce(hcat,getfield.(herds,:loc)))
    scatter(t[1], t[2],
    xlims=(-world[1]/2,world[1]/2),
    ylims=(-world[2]/2,world[2]/2),
    label=nothing)

    if diff == false    
        t = eachrow(reduce(hcat,getfield.(preds,:loc)))
        scatter!(t[1], t[2], label=nothing)
    else
        predsm1 = [p for p in preds if p.state == -1]
        preds0 = [p for p in preds if p.state == 0]
        preds1 = [p for p in preds if p.state == 1]
        preds2 = [p for p in preds if p.state == 2 && isnothing(p.following)]
        preds3 = [p for p in preds if p.state == 2 && !isnothing(p.following)]


        tm1 = [[],[]]
        t0 = [[],[]]
        t1 = [[],[]]
        t2 = [[],[]]
        t3 = [[],[]]
        if length(predsm1) != 0
            tm1 = eachrow(reduce(hcat,getfield.(predsm1,:loc)))
        end

        if length(preds0) != 0
            t0 = eachrow(reduce(hcat,getfield.(preds0,:loc)))
        end

        if length(preds1) != 0
            t1 = eachrow(reduce(hcat,getfield.(preds1,:loc)))
        end

        if length(preds2) != 0
            t2 = eachrow(reduce(hcat,getfield.(preds2,:loc)))
        end

        if length(preds3) != 0
            t3 = eachrow(reduce(hcat,getfield.(preds3,:loc)))
        end

        scatter!(tm1[1], tm1[2], label=nothing)
        scatter!(t0[1], t0[2], label=nothing)
        scatter!(t1[1], t1[2], label=nothing)
        scatter!(t2[1], t2[2], label=nothing)
        scatter!(t3[1], t3[2], color=:red,label=nothing)
    end

end

function sim(nherds, npreds; herdmean=40, world = [10,10],
    herdspeed = 3, herdcool = 2/24, predspeed = 5, predrange = 1, predcool=16/24,
    herdgen="exact", comm = true, thresh=1/4)
    # initialize herds
    herds, preds = initSim(nherds, npreds; herdmean, world, herdspeed, herdcool, predspeed, predrange, predcool, herdgen)
    threshhold = population(herds) * thresh
    steps = 0
    while population(herds) > threshhold
        steptime!(herds, preds; predcomm= comm)
        herds = removeempty(herds)
        steps += 1
    end
    return steps * Δt
end