

include("newton_basins.jl")
cmap = :YlOrRd_3
print_fig(600, 500, cmap) 


include("ott_basins.jl")
cmap = :YlOrRd_3
print_fig(600, 500, cmap) 

include("duffing_basins.jl")
cmap = :YlOrRd_3
print_fig(600, 500, cmap) 


include("pendulum_basins.jl")
cmap = :YlOrRd_3
print_fig(600, 500, cmap) 
