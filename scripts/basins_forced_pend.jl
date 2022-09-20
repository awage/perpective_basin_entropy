using DrWatson
@quickactivate "PerspectiveFigures" # exports DynamicalSystems, GLMakie and other goodies in `src`

using OrdinaryDiffEq
using DynamicalSystems
using CairoMakie
# using Plots


# Equations of motion:
function forced_pendulum(u, p, t)
    @inbounds begin
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] - sin(u[1])+ F*cos(omega*t)
    return SVector{2}(du1, du2)
    end
end

# We have to define a callback to wrap the phase in [-π,π]
function affect!(integrator)
    uu = integrator.u
    if integrator.u[1] < 0
        set_state!(integrator, SVector(uu[1] + 2π, uu[2]))
        u_modified!(integrator, true)
    else
        set_state!(integrator, SVector(uu[1] - 2π, uu[2]))
        u_modified!(integrator, true)
    end
end


function compute_basins_pend(d, F, ω, res)
    df = ContinuousDynamicalSystem(forced_pendulum,rand(2), [d, F, ω])
    condition(u,t,integrator) = (integrator.u[1] < -π  || integrator.u[1] > π)
    cb = DiscreteCallback(condition,affect!)
    diffeq = (reltol = 1e-9,  alg = Vern9(), callback = cb)
    xg = range(-pi,pi,length = res)
    yg = range(-4.,4.,length = res)
    smap = stroboscopicmap(df, 2*pi/ω; diffeq)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    bsn, att = basins_of_attraction(mapper)
    return bsn, att, (xg,yg)
end

function makesim(di::Dict)
    @unpack d, F, ω, res = di
    bsn, att, grid = compute_basins_pend(d, F, ω, res)
    return @strdict(bsn, att, grid, d, F, ω, res)
end

d = 0.2; F = 1.3636363636363635; ω = 0.5 # Parameters for Riddled Basins
res = 2000
params = @strdict d F ω res

data, file = produce_or_load(
    datadir("basins"), params, makesim;
    prefix = "pendulum", storepatch = false, suffix = "jld2", force = false
)


@unpack bsn, grid = data

# xg, yg = grid
xg = range(-pi,pi,length = res)
yg = range(-4.,4.,length = res)

fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], ylabel = "y", xlabel = "x", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)
heatmap!(ax, xg, yg, bsn', rasterize = 1)
save("basins_pendulum.svg",fig)



d = 0.2; F = 1.66; ω = 1. # Parameters for Wada Basins
res = 200
params = @strdict d F ω res

data, file = produce_or_load(
    datadir("basins"), params, makesim;
    prefix = "pendulum", storepatch = false, suffix = "jld2", force = false
)


@unpack bsn, grid = data

# xg, yg = grid
xg = range(-pi,pi,length = res)
yg = range(-4.,4.,length = res)

fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], ylabel = "y", xlabel = "x", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)
heatmap!(ax, xg, yg, bsn', rasterize = 1)
save("basins_pendulum_wada.svg",fig)
