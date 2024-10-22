# A simple demo that can be run without jupyter notebooks.
# First, activate the project (preferably after loading Revise), then include this file:
# > using Revise
# ] activate .
# > include("example/simple_demo.jl")
using SOLPS2imas: SOLPS2imas
using IMASggd: IMASggd
using Plots
using SOLPS2ctrl

sample_path = "$(@__DIR__)/../sample/ITER_Lore_2296_00000/"
SOLPS2imas_samples = splitdir(pathof(SOLPS2imas))[1] * "/../samples"

dd = SOLPS2ctrl.preparation(
    "Baseline2008-li0.70.x4.mod2.eqdsk",
    [sample_path, SOLPS2imas_samples]...;
    eqdsk_set_time=0.0,
)

grid_ggd = dd.edge_profiles.grid_ggd[1] # First grid_ggd time slice. It is allowed to vary in time
space = grid_ggd.space[1] # First space in this grid_ggd

# Choose backend
gr()           # Fast and can save pdf
# plotlyjs()   # Use for interactive plot, can only save png
n_e = IMASggd.get_prop_with_grid_subset_index(
    dd.edge_profiles.ggd[1].electrons.density,
    5,
)
plot(dd.edge_profiles.grid_ggd, n_e; colorbar_title="Electrons density / m^(-3)")
plot!(
    space,
    IMASggd.get_grid_subset(grid_ggd, 16);
    linecolor=:black,
    linewidth=2,
    linestyle=:solid,
    label="Separatix",
    legend=true,
)
