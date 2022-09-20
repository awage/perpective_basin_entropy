using DrWatson
@quickactivate "PerspectiveFigures" # exports DynamicalSystems, GLMakie and other goodies in `src`
using OrdinaryDiffEq:Vern9
using DynamicalSystems
using CairoMakie
# Equations of motion: E. Ott, et al. I Physica D 76 (1994) 384-410
function forced_particle!(du, u, p, t)
    γ=0.05  ; x̄ = 1.9  ; f₀=2.3  ; ω =3.5
    x₀=1. ; y₀=0.;
    x, y, dx, dy = u
    du[1] = dx
    du[2] = dy
    du[3] = -γ*dx -(-4*x*(1-x^2) + y^2) +  f₀*sin(ω*t)*x₀
    du[4] = -γ*dy -(2*y*(x+x̄)) +  f₀*sin(ω*t)*y₀
end


function _get_basins_ott(d)
    @unpack res = d
    xg = range(-2,2,length=res)
    yg = range(0.,2,length=res)
    df = ContinuousDynamicalSystem(forced_particle!,rand(4),(0.0,20.0))
    diffeq = (reltol = 1e-9,  alg = Vern9())
    ω =3.5
    smap = stroboscopicmap(df, 2π/ω; diffeq)
    psys = projected_integrator(smap, [1,2], [0., 0,])
    mapper = AttractorsViaRecurrences(psys, (xg, yg); horizon_limit = 10)
    basins, att = basins_of_attraction(mapper)
    return @strdict(basins, xg, yg)
end

res = 1500

data, file = produce_or_load(
    datadir("basins"), # path
    @dict(res), # container
    _get_basins_ott, # function
    prefix = "basin_ott", # prefix for savename
    force = false
)

@unpack basins, xg, yg = data
xg = range(-2,2,length = res)
yg = range(0.,2,length = res)


function print_fig(w,h,cmap)
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"y", xlabel = L"x", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    heatmap!(ax, xg, yg, basins, rasterize = 1, colormap = cmap)
    save("basins_riddle_ott.png",fig)
end