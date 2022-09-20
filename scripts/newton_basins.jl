using DrWatson
@quickactivate "PerspectiveFigures" # exports DynamicalSystems, GLMakie and other goodies in `src`

using DynamicalSystems
using CairoMakie

function newton_map(dz,z, p, n)
    f(x) = x^p[1]-1
    df(x)= p[1]*x^(p[1]-1)
    z1 = z[1] + im*z[2]
    dz1 = f(z1)/df(z1)
    z1 = z1 - dz1
    dz[1]=real(z1)
    dz[2]=imag(z1)
    return
end

# dummy function to keep the initializator happy
function newton_map_J(J,z0, p, n)
   return
end

function compute_basins_newton(N, res)
    ds = DiscreteDynamicalSystem(newton_map,[0.1, 0.2], [N] , newton_map_J)
    xg=range(-3.,3.,length = res)
    yg=range(-3.,3.,length = res)
    mapper = AttractorsViaRecurrences(ds, (xg, yg))
    bsn, att = basins_of_attraction(mapper)
    return bsn, att, (xg,yg)
end

function makesim(di::Dict)
    @unpack N, res = di
    bsn, att, grid = compute_basins_newton(N, res)
    return @strdict(bsn, att, grid, N, res)
end

res = 2000
N = 3
params = @strdict N res

data, file = produce_or_load(
    datadir("basins"), params, makesim;
    prefix = "newton", storepatch = false, suffix = "jld2", force = false
)

@unpack bsn, grid = data
function print_fig(w,h,cmap)
    xg, yg = grid
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"y", xlabel = L"x", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    save("basins_newton.png",fig)
end