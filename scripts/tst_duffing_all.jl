using DrWatson
@quickactivate "PerspectiveFigures" # exports DynamicalSystems, GLMakie and other goodies in `src`
using DynamicalSystems
using OrdinaryDiffEq:Vern9
using CairoMakie

@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end



function compute_basins_duffing(d, F, ω, res)
    ds = ContinuousDynamicalSystem(duffing, rand(2), [d, F, ω])
    xg = yg = range(-2.2,2.2,length = res)
    diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e6)
    smap = stroboscopicmap(ds, 2*pi/ω; diffeq)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    bsn, att = basins_of_attraction(mapper)
    return bsn, att, (xg,yg)
end

function makesim(di::Dict)
    @unpack d, F, ω, res = di
    bsn, att, grid = compute_basins_duffing(d, F, ω, res)
    return @strdict(bsn, att, grid, d, F, ω, res)
end

#F=0.3818791946308725; ω= 0.1966442953020134
# F=0.2771812080536913; ω=0.1;  # smooth boundary
ω=0.1617;F = 0.395; d = 0.15
# ω=0.3;F = 0.1

res = 200
params = @strdict d F ω res

data, file = produce_or_load(
    datadir("basins"), params, makesim;
    prefix = "duffing", storepatch = false, suffix = "jld2", force = false
)


@unpack bsn, grid = data

xg, yg = grid

fig = Figure(resolution = (1024, 768))
ax = Axis(fig[1,1], ylabel = "y", xlabel = "x", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabelsize = 20)
heatmap!(ax, xg, yg, bsn', rasterize = 4)
# save("diag_bif_henon.svg",fig)
save("basins_duffing.svg",fig)


