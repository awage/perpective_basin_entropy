using DrWatson
@quickactivate "PerspectiveFigures" # exports DynamicalSystems, GLMakie and other goodies in `src`
using Colors
using ColorSchemes

cmap = ColorScheme([RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )
include("newton_basins.jl")
print_fig(600, 600, cmap) 
@show get_Sb(3, 2000)

include("ott_basins.jl")
cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
print_fig(600, 500, cmap, 1200) 
@show get_Sb(1200)

include("duffing_basins.jl")
cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
d = 0.1; F=0.1; ω=0.1;  # smooth boundary
print_fig(600, 500, cmap, d, F, ω, 600) 
@show get_Sb(d, F, ω, 1500)
d = 0.4; F=0.1; ω=0.1;  # smooth boundary
print_fig(600, 500, cmap, d, F, ω, 600) 
@show get_Sb(d, F, ω, 1000)

include("pendulum_basins.jl")
d = 0.2; F = 1.3636363636363635; ω = 0.5 # Parameters for Riddled Basins
cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
print_fig(600, 500, cmap, d, F, ω, 2000)
@show get_Sb(d, F, ω, 2000)
d = 0.2; F = 1.66; ω = 1. # Parameters for Wada Basins
cmap = ColorScheme([RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )
print_fig(600, 500, cmap, d, F, ω, 500)
@show get_Sb(d, F, ω, 1500)

include("cuencas_hh.jl")
print_fig(600, 600, 0.25, 800) 
print_fig(600, 600, 0.2, 800) 
