"""
Utilities for extrapolating profiles
"""

using Interpolations: Interpolations
using IMASggd:
    IMASggd, get_grid_subset, add_subset_element!, get_subset_boundary,
    project_prop_on_subset!, get_subset_centers, get_TPS_mats, subset_do,
    all__grid_subset_prop
using PolygonOps: PolygonOps
using JSON: JSON

export extrapolate_core,
    fill_in_extrapolated_core_profile!, mesh_psi_spacing, cached_mesh_extension!

"""
    cumul_integrate(
        x::AbstractVector,
        y::AbstractVector,
    )::AbstractVector

Computes cumulative integral of y(x) using trapezoidal rule.
Source code from obsolete NumericalIntegration.jl package.
https://github.com/dextorious/NumericalIntegration.jl/blob/540b6c4bee089dfef7b9ae46e8ff188382e7c42e/src/NumericalIntegration.jl#L290
Modified the code to remove use of @inbounds, @fastmath macros that Julia documentation
recommends to use with caution.
"""
function cumul_integrate(
    x::AbstractVector,
    y::AbstractVector,
)::AbstractVector
    init = (x[2] - x[1]) * (y[1] + y[2])
    n = length(x)
    retarr = Vector{typeof(init)}(undef, n)
    retarr[1] = init
    for i ∈ 2:n
        retarr[i] = retarr[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1])
    end

    return 0.5 * retarr
end

"""
    extrapolate_core(
        edge_rho::Vector{Float64},
        edge_quantity::Vector{Float64},
        rho_output::Vector{Float64},
    )::Vector{Float64}

Function for assuming a core profile when given edge profile data.

Concept:

 1. value and derivative should be continuous when joining real data
 2. derivative at magnetic axis is known to be 0 when making profile vs. rho,
    by the definition of rho
 3. derivative probably does something fancier between the pedestal and axis than
    just linear interpolation, so add an extra point in there
 4. there's a joint between the steep pedestal and the shallow core that needs an
    extra knot to manage it properly
 5. after making up a very simple gradient profile out of a few line segments,
    integrate it to get the profile of the quantity in question
"""
function extrapolate_core(
    edge_rho::Vector{Float64},
    edge_quantity::Vector{Float64},
    rho_output::Vector{Float64},
)::Vector{Float64}
    grad = IMAS.gradient(edge_rho, edge_quantity)
    gf = grad[1]
    rf = edge_rho[1]
    gmid = -abs(gf) / 4.0
    rmid = rf / 2.0
    rped_enforce = rf - 0.08
    gped_enforce = gf + (gmid - gf) / (rmid - rf) * rped_enforce
    gped_max = maximum(grad) / 2.0
    gped_enforce = minimum([abs(gped_enforce), abs(gped_max)]) * sign(gf)

    gg = [0, gmid, gped_enforce, gf]
    rr = [0, rmid, rped_enforce, rf]
    # https://www.matecdev.com/posts/julia-interpolation.html
    itp = Interpolations.linear_interpolation(
        rr,
        gg;
        extrapolation_bc=Interpolations.Line(),
    )
    #itp = Interpolations.extrapolate(itp, Interpolations.Line())
    g = itp(rho_output)
    q_extend = cumul_integrate(rho_output, g)
    q_offset =
        edge_quantity[1] - Interpolations.linear_interpolation(rho_output, q_extend)(rf)
    q_extend .+= q_offset

    output_profile = Array{Float64}(undef, length(rho_output))
    output_profile[rho_output .< rf] = q_extend[rho_output .< rf]
    output_profile[rho_output .>= rf] =
        Interpolations.linear_interpolation(
            edge_rho,
            edge_quantity,
        ).(
            rho_output[rho_output .>= rf],
        )
    return output_profile
end

"""
    fill_in_extrapolated_core_profile!(
        dd::IMAS.dd,
        quantity_name::String;
        method::String="simple",
        eq_time_idx::Int=1,
        eq_profiles_2d_idx::Int=1,
        value_field::Symbol=:values,
        grid_ggd_idx::Int=1,
        cell_subset_idx::Int=5,
    )

This function accepts a DD that should be populated with `equilibrium` and
`edge_profiles` as well as a request for a quantity to extrapolate into the core. It
then maps `edge_profiles` data to rho, calls the function that performs the extrapolation (which is
not a simple linear extrapolation but has some trickery to attempt to make a somewhat
convincing profile shape), and writes the result to core_profiles. This involves a bunch
of interpolations and stuff.

Input arguments:

  - `dd`: an IMAS data dictionary
  - `quantity_name`: the name of a quantity in `edge_profiles.profiles_2d` and
    `core_profiles.profiles_1d`, such as "electrons.density"
  - `method`: Extrapolation method.
  - `eq_time_id`x: index of the equilibrium time slice to use. For a typical SOLPS run,
    the SOLPS mesh will be based on the equilibrium reconstruction at a single
    time, so the DD associated with the SOLPS run only needs one equilibrium
    time slice to be loaded. However, one could combine the complete
    equilibrium time series with the SOLPS run and then have to specify which
    slice of the equilibrium corresponds to the SOLPS mesh.
  - `eq_profiles_2d_idx`: index of the `profiles_2D` in equilibrium `time_slice`.
  - `value_field`: Symbolic name of the values field in quantity.
  - `grid_ggd_idx`: index of the `grid_ggd` to use. For a typical SOLPS run, the SOLPS
    grid is fixed, so this index defaults to 1. But in future, if a time varying grid is
    used, then this index will need to be specified.
  - `cell_subset_idx`: index of the subset of cells to use for the extrapolation. The
    default is 5, which is the subset of all cells. If `edge_profiles` data is instead
    present for a different subset, for instance, -5, which are b2.5 cells only, then
    this index should be set to -5.
"""
function fill_in_extrapolated_core_profile!(
    dd::IMAS.dd,
    quantity_name::String;
    method::String="simple",
    eq_time_idx::Int=1,
    eq_profiles_2d_idx::Int=1,
    value_field::Symbol=:values,
    grid_ggd_idx::Int=1,
    cell_subset_idx::Int=5,
)
    check_rho_1d(dd; time_slice=eq_time_idx, throw_on_fail=true)
    grid_ggd = dd.edge_profiles.grid_ggd[grid_ggd_idx]
    cell_subset = get_grid_subset(grid_ggd, cell_subset_idx)
    midplane_subset = get_grid_subset(grid_ggd, 11)

    if length(midplane_subset.element) < 1
        throw(
            ArgumentError(
                string(
                    "Midplane subset length was ", length(midplane_subset.element),
                    ". Unacceptable data in data dictionary.",
                ),
            ),
        )
    end
    if length(cell_subset.element) < 1
        throw(
            ArgumentError(
                string(
                    "Cell subset length was ", length(cell_subset.element),
                    ". Unacceptable data in data dictionary.",
                ),
            ),
        )
    end

    nt = length(dd.edge_profiles.ggd)
    if length(dd.core_profiles.profiles_1d) < nt
        resize!(dd.core_profiles.profiles_1d, nt)
    end
    TPS_mats = get_TPS_mats(grid_ggd, cell_subset_idx)

    mcc_qty = project_prop_on_subset!(
        dd.edge_profiles.ggd,
        quantity_name,
        cell_subset,
        midplane_subset;
        value_field,
        TPS_mats)

    for it ∈ 1:nt
        midplane_cell_centers, quantity = mcc_qty[it]
        # Now quantity is at the outboard midplane

        # Get the rho values to go with the midplane quantity values
        r = [midplane_cell_centers[i][1] for i ∈ eachindex(midplane_cell_centers)]
        z = [midplane_cell_centers[i][2] for i ∈ eachindex(midplane_cell_centers)]
        if length(r) != length(quantity)
            throw(
                DimensionMismatch(
                    string(
                        "Number of cell center coordinates (", length(r),
                        ") does not match number of cell values (", length(quantity),
                        ") in the midplane subset.",
                    ),
                ),
            )
        end
        eq_time_slice = dd.equilibrium.time_slice[eq_time_idx]
        eq_prof_2d = eq_time_slice.profiles_2d[eq_profiles_2d_idx]
        r_eq = eq_prof_2d.grid.dim1
        z_eq = eq_prof_2d.grid.dim2
        rho1_eq = eq_time_slice.profiles_1d.rho_tor_norm
        psia = eq_time_slice.global_quantities.psi_axis
        psib = eq_time_slice.global_quantities.psi_boundary
        psi1_eq = (eq_time_slice.profiles_1d.psi .- psia) ./ (psib - psia)
        psi2_eq = (eq_prof_2d.psi .- psia) ./ (psib - psia)
        # println(size(psi2_eq), ", ", size(r_eq), ", ", size(z_eq))
        rzpi = Interpolations.linear_interpolation((r_eq, z_eq), psi2_eq)
        in_bounds =
            (r .< maximum(r_eq)) .& (r .> minimum(r_eq)) .& (z .> minimum(z_eq)) .&
            (z .< maximum(z_eq))
        psi_for_quantity = 10.0 .+ zeros(length(r))
        psi_for_quantity[in_bounds] = rzpi.(r[in_bounds], z[in_bounds])
        rho_for_quantity = copy(psi_for_quantity)
        in_bounds = psi_for_quantity .<= 1.0
        dpsi = diff(psi1_eq)
        drho = diff(rho1_eq)
        prepend!(dpsi, [0.0])
        prepend!(drho, [0.0])
        rho_for_quantity[in_bounds] =
            Interpolations.linear_interpolation(
                psi1_eq,
                rho1_eq,
            ).(
                psi_for_quantity[in_bounds],
            )

        # Make sure the output 1D rho grid exists; create it if needed
        if IMAS.ismissing(dd.core_profiles.profiles_1d[it].grid, :rho_tor_norm)
            # If you don't like this default, then you should write grid.rho_tor_norm
            # before calling this function.
            dd.core_profiles.profiles_1d[it].grid.rho_tor_norm =
                collect(LinRange(0, 1, 201))
        end
        rho_core = dd.core_profiles.profiles_1d[it].grid.rho_tor_norm

        # Finally, we're ready to call the extrapolation function and write the result
        if method == "simple"
            quantity_core = extrapolate_core(rho_for_quantity, quantity, rho_core)
        else
            throw(ArgumentError(string(
                "Unrecognized extraplation method: ", method,
            )))
        end
        parent = dd.core_profiles.profiles_1d[it]
        tags = split(quantity_name, ".")
        for tag ∈ tags[1:(end-1)]
            parent = getproperty(parent, Symbol(tag))
        end
        setproperty!(parent, Symbol(tags[end]), quantity_core)
    end
end

"""
    extrapolate_edge_exp(
        quantity_edge::Vector{Float64},
        dqdpsi::Vector{Float64},
        psin_out::Vector{Float64},
    )::Matrix{Float64}

Exponential decay version of edge profile extrapolation. Should work well for
many quantities, including Te and ne. If the exponential profile has no background
offset, then its amplitude and scale length can be completely defined by matching
the quantity value and first derivative at the edge of the mesh.

Input Arguments:

  - `quantity_edge`: values of some physics quantity in cells along the outer edge of
    the mesh.

  - `dqdpsi`: Gradient of the quantity vs. psi, aligned perpendicular to the row of
    cells being used.
  - `psin_out`: Normalized psi values along a vector orthogonal to the row of cells
    along the edge. These `psi_N` values should be outside of the mesh (because the
    quantity is already known in the mesh).

The output will be a matrix.
"""
function extrapolate_edge_exp(
    quantity_edge::Vector{Float64},
    dqdpsi::Vector{Float64},
    psin_out::Vector{Float64},
)::Matrix{Float64}
    x = psin_out - 1.0
    lambda = -quantity_edge / dqdpsi
    q0 = quantity_edge ./ exp(-x' ./ lambda)
    return q0 * exp(-x ./ lambda)
end

"""
    prep_flux_map(dd::IMAS.dd; eq_time_idx::Int=1, eq_profiles_2d_idx::Int=1)

Reads equilibrium data and extracts/derives some useful quantities.
This is very basic, but it was being repeated and that's a no-no.
Returns:

  - R values of the equilibrium grid
  - Z values of the eq grid
  - normalized poloidal flux on the equilibrium grid
  - a linear interpolation of norm pol flux vs. R and Z, ready to be evaluated
"""
function prep_flux_map(dd::IMAS.dd; eq_time_idx::Int=1, eq_profiles_2d_idx::Int=1)
    eqt = dd.equilibrium.time_slice[eq_time_idx]
    p2 = eqt.profiles_2d[eq_profiles_2d_idx]
    r_eq = p2.grid.dim1
    z_eq = p2.grid.dim2
    psi = p2.psi
    psia = eqt.global_quantities.psi_axis
    psib = eqt.global_quantities.psi_boundary
    psin_eq = (psi .- psia) ./ (psib - psia)
    rzpi = Interpolations.linear_interpolation((r_eq, z_eq), psin_eq)
    return r_eq, z_eq, psin_eq, rzpi
end

"""
    mesh_psi_spacing(
        dd::IMAS.dd;
        eq_time_idx::Int=1,
        eq_profiles_2d_idx::Int=1,
        grid_ggd_idx::Int=1,
        space_idx::Int=1,
        avoid_guard_cell::Bool=true,
        spacing_rule="mean",
    )

Inspects the mesh to see how far apart faces are in psi_N.
Requires that GGD and equilibrium are populated.

Input Arguments:

  - `dd`: a data dictionary instance with required data loaded into it
  - `eq_time_idx`: index of the equilibrium time slice to use. For a typical SOLPS run,
    the SOLPS mesh will be based on the equilibrium reconstruction at a single
    time, so the DD associated with the SOLPS run only needs one equilibrium
    time slice to be loaded. However, one could combine the complete
    equilibrium time series with the SOLPS run and then have to specify which
    slice of the equilibrium corresponds to the SOLPS mesh.
  - `eq_profiles_2d_id`x: index of the `profiles_2D` in equilibrium `time_slice`.
  - `grid_ggd_idx`: index of the `grid_ggd` to use. For a typical SOLPS run, the SOLPS
    grid is fixed, so this index defaults to 1. But in future, if a time varying grid
    is used, then this index will need to be specified.
  - `space_idx`: index of the space to use. For a typical SOLPS run, there will be only
    one space so this index will mostly remain at 1.
  - `avoid_guard_cell`: assume that the last cell is a guard cell so take `end-2` and
    `end-1` instead of `end` and `end-1`
  - `spacing_rule`: "edge" or "mean" to make spacing of new cells (in `psi_N`) be the
    same as the spacing at the edge of the mesh, or the same as the average spacing
"""
function mesh_psi_spacing(
    dd::IMAS.dd;
    eq_time_idx::Int=1,
    eq_profiles_2d_idx::Int=1,
    grid_ggd_idx::Int=1,
    space_idx::Int=1,
    avoid_guard_cell::Bool=true,
    spacing_rule="mean",
)
    # Inspect input
    if length(dd.equilibrium.time_slice) < eq_time_idx
        throw(
            ArgumentError(
                string(
                    "DD equilibrium does not have enough time slices: ",
                    length(dd.equilibrium.time_slice),
                    " slices were present but slice index ",
                    eq_time_idx, " was requested.",
                ),
            ),
        )
    end
    bad_ggd = (length(dd.edge_profiles.grid_ggd) < grid_ggd_idx)
    if bad_ggd
        throw(ArgumentError(string(
            "Invalid GGD data.",
        )))
    end

    # Get flux map
    r_eq, z_eq, psin_eq, rzpi = prep_flux_map(dd; eq_time_idx, eq_profiles_2d_idx)

    # Get a row of cells. Since the mesh should be aligned to the flux surfaces,
    # it shouldn't matter which row is used, although the divertor rows might be
    # weird. So use the outboard midplane. That's always a solid choice.
    grid_ggd = dd.edge_profiles.grid_ggd[grid_ggd_idx]
    space = grid_ggd.space[space_idx]
    midplane_subset = get_grid_subset(grid_ggd, 11)
    midplane_cell_centers = get_subset_centers(space, midplane_subset)
    r_mesh = [midplane_cell_centers[i][1] for i ∈ eachindex(midplane_cell_centers)]
    z_mesh = [midplane_cell_centers[i][2] for i ∈ eachindex(midplane_cell_centers)]
    psin_mesh = rzpi.(r_mesh, z_mesh)
    # This should come out sorted, but with GGD, who knows.
    ii = sortperm(psin_mesh)
    psin = psin_mesh[ii]
    if spacing_rule == "edge"
        if avoid_guard_cell
            dpsin = psin[end-1] - psin[end-2]
        else
            dpsin = psin[end] - psin[end-1]
        end
    else
        if avoid_guard_cell
            dpsin = diff(psin[2:(end-1)])
        else
            dpsin = diff(psin)
        end
        dpsin = sum(dpsin) / length(dpsin)
    end
    return dpsin
end

"""
    pick_extension_psi_range(
        dd::IMAS.dd;
        eq_time_idx::Int=1,
        eq_profiles_2d_idx::Int=1,
        grid_ggd_idx::Int=1,
        space_idx::Int=1,
    )::Vector{Float64}

Defines the `psi_N` levels for an extended mesh. The range of `psi_N` levels starts
at the outer edge of the existing `edge_profiles` mesh at the midplane and goes
out to the most distant (in flux space) point on the limiting surface.
Returns a vector of `psi_N` levels.
"""
function pick_extension_psi_range(
    dd::IMAS.dd;
    eq_time_idx::Int=1,
    eq_profiles_2d_idx::Int=1,
    grid_ggd_idx::Int=1,
    space_idx::Int=1,
)::Vector{Float64}
    r_eq, z_eq, psin_eq, rzpi = prep_flux_map(dd; eq_time_idx, eq_profiles_2d_idx)

    # Use wall to find maximum extent of contouring project
    limiter = dd.wall.description_2d[1].limiter
    wall_r = limiter.unit[1].outline.r
    wall_z = limiter.unit[1].outline.z
    wall_psin = rzpi.(wall_r, wall_z)

    # Use ggd mesh to find inner limit of contouring project
    grid_ggd = dd.edge_profiles.grid_ggd[grid_ggd_idx]
    space = grid_ggd.space[space_idx]
    midplane_subset = get_grid_subset(grid_ggd, 11)
    midplane_cell_centers = get_subset_centers(space, midplane_subset)
    psin_midplane = rzpi.(midplane_cell_centers[end][1], midplane_cell_centers[end][2])

    # Choose contour levels
    # The psi spacing function doesn't do much that's unique anymore.
    # Should it be absorbed into this function intead?
    dpsin = mesh_psi_spacing(
        dd;
        eq_time_idx=eq_time_idx,
        eq_profiles_2d_idx=eq_profiles_2d_idx,
        grid_ggd_idx=grid_ggd_idx,
        space_idx=space_idx,
        avoid_guard_cell=true,
        spacing_rule="mean",
    )
    lvlstart = maximum(psin_midplane * sign(dpsin)) / sign(dpsin) + dpsin
    lvlend = maximum(wall_psin * sign(dpsin)) / sign(dpsin) + dpsin
    nlvl = Int(ceil((lvlend - lvlstart) / dpsin))
    psin_levels = collect(LinRange(lvlstart, lvlend, nlvl))
    # eqt = dd.equilibrium.time_slice[eq_time_idx]
    # if hasproperty(eqt.boundary_secondary_separatrix, :psi)
    #     secondary_psi = eqt.boundary_secondary_separatrix.psi
    #     secondary_psin = (secondary_psi - psia) / (psib - psia)
    # else
    #     secondary_psin = lvlend + dpsin
    # end
    return psin_levels
end

"""
    pick_mesh_ext_starting_points(
        grid_ggd::IMAS.edge_profiles__grid_ggd,
        space::IMAS.edge_profiles__grid_ggd___space,
    )::Tuple{Vector{Float64}, Vector{Float64}}

Picks starting points for the radial lines of the mesh extension. The strategy
is to start from the outer edge of the existing mesh and follow the steepest
gradient (of psi_N) to extend these gridlines outward.
Returns a tuple with vectors of R and Z starting points.
"""
function pick_mesh_ext_starting_points(
    grid_ggd::IMAS.edge_profiles__grid_ggd,
    space::IMAS.edge_profiles__grid_ggd___space,
)::Tuple{Vector{Float64}, Vector{Float64}}
    # Choose starting points for the orthogonal (to the contour) gridlines
    # Use the existing cells of the standard mesh
    all_cell_subset = get_grid_subset(grid_ggd, 5)
    all_cell_bdry = get_subset_boundary(space, all_cell_subset)
    core_bdry = get_grid_subset(grid_ggd, 15)
    outer_target = get_grid_subset(grid_ggd, 13)
    inner_target = get_grid_subset(grid_ggd, 14)
    wall_subset = subset_do(setdiff, all_cell_bdry, core_bdry)
    wall_subset = subset_do(setdiff, wall_subset, outer_target)
    wall_subset = subset_do(setdiff, wall_subset, inner_target)

    border_edges = wall_subset.element
    npol = length(border_edges)
    r = zeros(npol)
    z = zeros(npol)
    corner_idx = 1
    for i ∈ 1:npol
        edge_idx = border_edges[i].object[1].index
        node_idx = space.objects_per_dimension[2].object[edge_idx].nodes[corner_idx]
        geo = space.objects_per_dimension[1].object[node_idx].geometry
        r[i, 1] = geo[1]
        z[i, 1] = geo[2]
    end
    return r, z
end

"""
    mesh_ext_follow_grad(
        r_eq::Vector{Float64},
        z_eq::Vector{Float64},
        psin_eq::Matrix,
        rstart::Vector{Float64},
        zstart::Vector{Float64},
        nlvl::Int,
        dpsin::Float64,
        rzpi=nothing,
    )::Tuple{Matrix{Float64}, Matrix{Float64}}

Follows the steepest gradient from a set of starting points, dropping nodes at
approximately regular intervals in psi_N. Due to the numerical techniques used, the
node spacing may be imperfect (especially prone to error in regions where curvature
of `psi_N` is large compared to its gradient).

Input Arguments:

  - `r_eq`: Equilibrium reconstruction's grid, R coordinates
  - `z_eq`: Equilibrium reconstruction's grid, Z coordinates
  - `psin_eq`: Normalized poloidal flux in the equilibrium reconstruction as a function
    of R and Z
  - `rstar`t: R coordinates of starting points for the gradient following.
  - `zstart`: Z coordinates of starting points
  - `nlvl`: number of nodes to drop while following the gradient
  - `dpsin`: node spacing in delta psi_N
  - `rzpi`: linear interpolation of psin_eq() as a function of r_eq and z_eq. This was
    probably already computed and I think time would be saved by reusing it. If you
    don't already have it, you can pass in nothing and let this function calculate it.

Returns two matrices with R and Z coordinates of the mesh extension
"""
function mesh_ext_follow_grad(
    r_eq::Vector{Float64},
    z_eq::Vector{Float64},
    psin_eq::Matrix,
    rstart::Vector{Float64},
    zstart::Vector{Float64},
    nlvl::Int,
    dpsin::Float64,
    rzpi=nothing,
)::Tuple{Matrix{Float64}, Matrix{Float64}}
    npol = length(rstart)
    mesh_r = zeros((npol, nlvl))
    mesh_z = zeros((npol, nlvl))
    for i ∈ 1:npol
        mesh_r[i, 1] = rstart[i]
        mesh_z[i, 1] = zstart[i]
    end

    if rzpi === nothing
        rzpi = Interpolations.linear_interpolation((r_eq, z_eq), psin_eq)
    end

    # Step along the paths of steepest descent to populate the mesh.
    dpsindr, dpsindz = IMAS.gradient(r_eq, z_eq, psin_eq)
    dpdr = Interpolations.linear_interpolation((r_eq, z_eq), dpsindr)
    dpdz = Interpolations.linear_interpolation((r_eq, z_eq), dpsindz)
    rlim = (minimum(r_eq), maximum(r_eq))
    zlim = (minimum(z_eq), maximum(z_eq))
    pfr = rzpi.(mesh_r[:, 1], mesh_z[:, 1]) .< 1
    for i ∈ 1:npol
        if (mesh_r[i, 1] > rlim[1]) & (mesh_r[i, 1] < rlim[2]) &
           (mesh_z[i, 1] > zlim[1]) & (mesh_z[i, 1] < zlim[2])
            # direction = rzpi(mesh_r[i, 1], mesh_z[i, 1]) >= 1 ? 1 : -1
            direction = pfr[i] ? -1 : 1
            for j ∈ 2:nlvl
                # This is low resolution linear shooting that could go wrong
                mesh_r[i, j] = mesh_r[i, j-1]
                mesh_z[i, j] = mesh_z[i, j-1]
                upscale = 5  # Should improve gradient following in regions of low gradient (near X-point)
                for k ∈ 1:upscale
                    if (mesh_r[i, j] > rlim[1]) & (mesh_r[i, j] < rlim[2]) &
                       (mesh_z[i, j] > zlim[1]) & (mesh_z[i, j] < zlim[2])
                        dpr = dpdr(mesh_r[i, j], mesh_z[i, j])
                        dpz = dpdz(mesh_r[i, j], mesh_z[i, j])
                        d = dpsin * direction / upscale / (dpr .^ 2 + dpz .^ 2)
                        mesh_r[i, j] += d * dpr
                        mesh_z[i, j] += d * dpz
                    end
                end
            end
        end
        # println("i=",i,"; r=",mesh_r[i, 1],":", mesh_r[i, end],",z=",mesh_z[i, 1], ":", mesh_z[i, end])
    end
    return mesh_r, mesh_z
end

"""
    modify_mesh_ext_near_x!(
        eqt::IMAS.equilibrium__time_slice,
        mesh_r::Matrix{Float64},
        mesh_z::Matrix{Float64},
    )

Modifies an extended mesh near a secondary X-point to compensate for the
tendency of the mesh to go nuts near the X-point.

Input Arguments:

  - `eqt`: equilibrium.time_slice information
  - `mesh_r`: matrix of R values for the extended mesh
  - `mesh_z`: matrix of Z values for the extended mesh
"""
function modify_mesh_ext_near_x!(
    eqt::IMAS.equilibrium__time_slice,
    mesh_r::Matrix{Float64},
    mesh_z::Matrix{Float64},
)
    # There's a special path; the one that probably should've gone through the
    # secondary X-point (if there is one). Any numerical error will make this
    # path miss the X-point and go off to somewhere crazy instead, so instead of
    # that, let's draw a straight line from the edge of the SOLPS mesh to the
    # X-point.
    npol, nlvl = size(mesh_r)
    bssx = eqt.boundary_secondary_separatrix.x_point
    nsx = length(bssx)
    println("there are ", nsx, " secondary x points")
    if nsx > 0
        # There is hopefully only one X point on the secondary separatrix, but
        # just in case some joker made a secondary snowflake or something crazy,
        # we should pick which one we like best.
        if nsx == 1
            xidx = 1  # Easy peasy
        else
            score = zeros(nsx)
            rsx = [bssx[i].r for i ∈ 1:nsx]
            zsx = [bssx[i].z for i ∈ 1:nsx]
            # Being close to the vertically flipped position of the primary X-pt
            # seems nice; let's reward that
            bsx = eqt.boundary_separatrix.x_point
            npx = length(bsx)
            rpx = [bsx[i].r for i ∈ 1:npx]
            zpx = [bsx[i].z for i ∈ 1:npx]
            # Well, if we have a primary snowflake, there could be multiple
            # primary X-points, so now we have to pick our favorite primary.
            # Being close to the boundary is a good sign
            min_ds2 = ones(npx) * 1e6
            for j ∈ 1:nsx
                dr = [mesh_r[i, 1] - rpx[j] for i ∈ 1:npol]
                dz = [mesh_z[i, 1] - zpx[j] for i ∈ 1:npol]
                ds2 = dr .^ 2 .+ dz .^ 2
                min_ds2[j] = minimum(ds2)
            end
            primary_idx = argmin(min_ds2)
            rpx = rpx[primary_idx]
            zpx = zpx[primary_idx]
            # There can be only one

            # Distance between secondary X-points and vert flip of favorite primary
            ds2xx = (rsx .- rpx) .^ 2 .+ (zsx .- (-zpx)) .^ 2
            score += ds2xx * 5

            # Being close to the boundary is good
            min_ds2xb = ones(nsx) * 1e6
            for j ∈ 1:nsx
                dr = [mesh_r[i, 1] - rsx[j] for i ∈ 1:npol]
                dz = [mesh_z[i, 1] - zsx[j] for i ∈ 1:npol]
                ds2 = dr .^ 2 .+ dz .^ 2
                min_ds2xb[j] = minimum(ds2)
            end
            score += min_ds2xb
            # And the winner is the one with the lowest score
            xidx = argmin(score)
        end
        # Pick which edges are closest to the X-point
        rx = bssx[xidx].r
        zx = bssx[xidx].z
        ds2 = (mesh_r[:, 1] .- rx) .^ 2 .+ (mesh_z[:, 1] .- zx) .^ 2
        closest = argmin(ds2)
        # Work out the mesh spacing
        drm = mesh_r[closest, 2] - mesh_r[closest, 1]
        dzm = mesh_z[closest, 2] - mesh_z[closest, 1]
        dm = sqrt(drm^2 + dzm^2)
        # Pick the angle
        dr = rx - mesh_r[closest, 1]
        dz = zx - mesh_z[closest, 1]
        angle = atan(dz, dr)
        # Replace the points
        for i ∈ 2:nlvl
            mesh_r[closest, i] = mesh_r[closest, 1] + dm * (i - 1) * cos(angle)
            mesh_z[closest, i] = mesh_z[closest, 1] + dm * (i - 1) * sin(angle)
        end
    end
end

"""
    record_regular_mesh!(
        grid_ggd::IMAS.edge_profiles__grid_ggd,
        space::IMAS.edge_profiles__grid_ggd___space,
        mesh_r::Matrix{Float64},
        mesh_z::Matrix{Float64},
        cut::Int,
    )

Records arrays of mesh data from regular 2D arrays into the DD

  - `grid_ggd`: `grid_ggd` within edge_profiles
  - `space`: space in edge_profiles
  - `mesh_r`: Matrix of R values along a mesh. Should be 2D. The two dimensions are in
    the radial and poloidal directions.
  - `mesh_z`: Z values to go with mesh_r.
  - `cut`: Poloidal index of a cut between two groups of disconnected cells. Poloidal
    connections (faces, cells) will not be added between this poloidal index
    and the next index.
"""
function record_regular_mesh!(
    grid_ggd::IMAS.edge_profiles__grid_ggd,
    space::IMAS.edge_profiles__grid_ggd___space,
    mesh_r::Matrix{Float64},
    mesh_z::Matrix{Float64},
    cut::Int,
)
    npol, nlvl = size(mesh_r)
    o0 = space.objects_per_dimension[1]  # Nodes
    o1 = space.objects_per_dimension[2]  # Edges
    o2 = space.objects_per_dimension[3]  # Cells (2D projection)
    node_dim = 1
    edge_dim = 2
    cell_dim = 3

    n_old_node = length(o0.object)
    n_new_node = nlvl * npol
    n_nodes = n_old_node + n_new_node

    n_old_edge = length(o1.object)
    n_new_edge_i = nlvl * (npol - 2)
    n_new_edge_j = (nlvl - 1) * npol
    n_new_edge = n_new_edge_i + n_new_edge_j
    n_edges = n_old_edge + n_new_edge

    n_old_cell = length(o2.object)
    n_new_cell = (nlvl - 1) * (npol - 2)  # -2 to account for disconnect @ pfr
    n_cells = n_old_cell + n_new_cell

    # Define starting points and increments
    node_start = n_old_node + 1
    edge_start1 = n_old_edge + 1
    edge_start2 = edge_start1 + n_new_edge_i
    cell_start = n_old_cell + 1

    n_per_i = nlvl
    n_per_j = 1
    e1_per_i = nlvl  # Edges along the npol direction / i direction
    e2_per_i = nlvl - 1  # Edges along the nlvl direction / j direction
    c_per_i = nlvl - 1  # # of cells < # of nodes

    # Get new subsets ready
    n_existing_subsets = length(grid_ggd.grid_subset)
    new_subset_indices = [-1, -2, -3, -4, -5, -201, -202, -203, -204, -205]
    n_new_subsets = length(new_subset_indices)
    resize!(grid_ggd.grid_subset, n_existing_subsets + n_new_subsets)
    for i ∈ 1:n_new_subsets
        grid_ggd.grid_subset[n_existing_subsets+i].identifier.index =
            new_subset_indices[i]
    end
    ext_nodes_sub = get_grid_subset(grid_ggd, -201)
    ext_edges_sub = get_grid_subset(grid_ggd, -202)
    ext_xedges_sub = get_grid_subset(grid_ggd, -203)
    ext_yedges_sub = get_grid_subset(grid_ggd, -204)
    ext_cells_sub = get_grid_subset(grid_ggd, -205)

    # Preserve record of standard (non extended) mesh
    for i ∈ 1:5
        orig_sub = get_grid_subset(grid_ggd, i)
        grid_ggd.grid_subset[n_existing_subsets+i] = deepcopy(orig_sub)
        grid_ggd.grid_subset[n_existing_subsets+i].identifier.index = -i
    end
    all_nodes_sub = get_grid_subset(grid_ggd, 1)
    all_edges_sub = get_grid_subset(grid_ggd, 2)
    all_xedges_sub = get_grid_subset(grid_ggd, 3)
    all_yedges_sub = get_grid_subset(grid_ggd, 4)
    all_cells_sub = get_grid_subset(grid_ggd, 5)

    nodes = resize!(o0.object, n_nodes)
    edges = resize!(o1.object, n_edges)
    cells = resize!(o2.object, n_cells)
    space_idx = space.identifier.index
    for i ∈ 1:npol
        pastpc = i > cut
        for j ∈ 1:nlvl
            # Modified counters
            ii = i - 1  # Offset due to indexing from 1
            jj = j - 1  # Offset due to indexing from 1
            iii = ii - 1 - pastpc  # # of connections < than # of nodes, missing row @ PFR transition
            jjj = jj - 1  # # of connections < # of nodes

            # Nodes
            node_idx = node_start + ii * n_per_i + jj
            nodes[node_idx].geometry = [mesh_r[i, j], mesh_z[i, j]]
            add_subset_element!(ext_nodes_sub, space_idx, node_dim, node_idx)
            add_subset_element!(all_nodes_sub, space_idx, node_dim, node_idx)

            # Edges
            if (i > 1) & (i != (cut + 1))  # i-1 to i  in the npol direction
                edge_idx1 = edge_start1 + iii * e1_per_i + jj
                edges[edge_idx1].nodes = [node_idx, node_idx - n_per_i]
                resize!(edges[edge_idx1].boundary, 2)
                edges[edge_idx1].boundary[1].index = node_idx
                edges[edge_idx1].boundary[2].index = node_idx - n_per_i
                add_subset_element!(ext_edges_sub, space_idx, edge_dim, edge_idx1)
                add_subset_element!(ext_xedges_sub, space_idx, edge_dim, edge_idx1)
                add_subset_element!(all_edges_sub, space_idx, edge_dim, edge_idx1)
                add_subset_element!(all_xedges_sub, space_idx, edge_dim, edge_idx1)
            end
            if (j > 1)  # j-1 to j in the nlvl direction
                edge_idx2 = edge_start2 + ii * e2_per_i + jjj
                edges[edge_idx2].nodes = [node_idx, node_idx - n_per_j]
                resize!(edges[edge_idx2].boundary, 2)
                edges[edge_idx2].boundary[1].index = node_idx
                edges[edge_idx2].boundary[2].index = node_idx - n_per_j
                add_subset_element!(ext_edges_sub, space_idx, edge_dim, edge_idx2)
                add_subset_element!(ext_yedges_sub, space_idx, edge_dim, edge_idx2)
                add_subset_element!(all_edges_sub, space_idx, edge_dim, edge_idx2)
                add_subset_element!(all_yedges_sub, space_idx, edge_dim, edge_idx2)
            end

            # Cells
            if (i > 1) & (i != (cut + 1)) & (j > 1)
                cell_idx = cell_start + iii * c_per_i + jjj
                cells[cell_idx].nodes = [
                    node_idx,
                    node_idx - n_per_i,
                    node_idx - n_per_i - n_per_j,
                    node_idx - n_per_j,
                ]
                resize!(cells[cell_idx].boundary, 4)
                cells[cell_idx].boundary[1].index = edge_idx1
                cells[cell_idx].boundary[2].index = edge_idx2
                cells[cell_idx].boundary[3].index = edge_idx1 - 1
                cells[cell_idx].boundary[4].index = edge_idx2 - e2_per_i

                add_subset_element!(ext_cells_sub, space_idx, cell_dim, cell_idx)
                add_subset_element!(all_cells_sub, space_idx, cell_dim, cell_idx)
            end
        end
    end
end

"""
    convert_filename(filename::String)::String

Converts a filename into a string that doesn't have illegal characters.
The main application is removing the path separator from source files with full
paths so the full paths* can be part of the new filename. This way, the input
files used to form some data can be part of the cache name, allowing quick lookup:
the cache filename is defined by the input files, and if it doesn't exist, it
needs to be generated.
"""
function convert_filename(filename::String)::String
    filename_mod = replace(filename, "/" => "__")  # Illegal on *nix bc it's the path separator
    filename_mod = replace(filename_mod, ":" => "--")  # Illegal on mac and windows
    filename_mod = replace(filename_mod, "\\" => "__")  # Illegal on windows
    filename_mod = replace(filename_mod, "\"" => "''")  # Illegal on windows
    filename_mod = replace(filename_mod, "<" => "_lt_")  # Illegal on windows
    filename_mod = replace(filename_mod, ">" => "_gt_")  # Illegal on windows
    filename_mod = replace(filename_mod, "|" => "_pipe_")  # Illegal on windows
    filename_mod = replace(filename_mod, "?" => "_q_")  # Illegal on windows
    filename_mod = replace(filename_mod, "*" => "_a_")  # Illegal on windows
    return filename_mod
end

"""
    cached_mesh_extension!(
        dd::IMAS.dd,
        eqdsk_file::String,
        b2fgmtry::String;
        eq_time_idx::Int=1,
        eq_profiles_2d_idx::Int=1,
        grid_ggd_idx::Int=1,
        space_idx::Int=1,
        clear_cache::Bool=false,
    )::String

Adds an extended mesh to a data dictionary, possibly from a cached result.

Input Arguments:

  - `dd`: The data dictionary. It will be modified in place.
  - `eqdsk_file`: the name of the EQDSK file that was used to get equilibrium data in
    the dd.
  - `b2fgmtry`: the name of the SOLPS geometry file that was used to get GGD info in
    `edge_profiles` in the dd.
  - `eq_time_idx`: Index of the time slice in equilibrium
  - `eq_profiles_2d_idx`: Index of the 2D profile set in equilibrium
    (there is usually only one)
  - `grid_ggd_idx`: Index of the `grid_ggd` set in edge_profiles
  - `space_idx`: Index of the space
  - `clear_cache`: delete any existing cache file (for use in testing)
"""
function cached_mesh_extension!(
    dd::IMAS.dd,
    eqdsk_file::String,
    b2fgmtry::String;
    eq_time_idx::Int=1,
    eq_profiles_2d_idx::Int=1,
    grid_ggd_idx::Int=1,
    space_idx::Int=1,
    clear_cache::Bool=false,
)::String
    path = "$(@__DIR__)/../data/"
    cached_ext_name = path * string(hash(eqdsk_file * b2fgmtry)) * ".mesh_ext.json"
    if clear_cache
        rm(cached_ext_name; force=true)
        return cached_ext_name
    end
    if isfile(cached_ext_name)
        md = JSON.parsefile(cached_ext_name)
        pfr_transition = md["pfr_transition"]
        mesh_r = convert(Matrix{Float64}, mapreduce(permutedims, vcat, md["r"])')
        mesh_z = convert(Matrix{Float64}, mapreduce(permutedims, vcat, md["z"])')
        npol = md["npol"]
        nlvl = md["nlvl"]
        record_regular_mesh!(
            dd.edge_profiles.grid_ggd[grid_ggd_idx],
            dd.edge_profiles.grid_ggd[grid_ggd_idx].space[space_idx],
            mesh_r,
            mesh_z,
            pfr_transition,
        )
    else
        mesh_r, mesh_z, pfr_transition = mesh_extension_sol!(
            dd;
            eq_time_idx=eq_time_idx,
            eq_profiles_2d_idx=eq_profiles_2d_idx,
            grid_ggd_idx=grid_ggd_idx,
            space_idx=space_idx,
        )
        npol, nlvl = size(mesh_r)
        data = Dict()
        data["pfr_transition"] = pfr_transition
        data["npol"] = npol
        data["nlvl"] = nlvl
        data["r"] = mesh_r
        data["z"] = mesh_z
        # YAML.write_file(cached_ext_name, data)
        open(cached_ext_name, "w") do f
            return JSON.print(f, data)
        end
    end
    return cached_ext_name
end

"""
    mesh_extension_sol!(
        dd::IMAS.dd;
        eq_time_idx::Int=1,
        eq_profiles_2d_idx::Int=1,
        grid_ggd_idx::Int=1,
        space_idx::Int=1,
    )

Extends the mesh out into the SOL
"""
function mesh_extension_sol!(
    dd::IMAS.dd;
    eq_time_idx::Int=1,
    eq_profiles_2d_idx::Int=1,
    grid_ggd_idx::Int=1,
    space_idx::Int=1,
)
    grid_ggd = dd.edge_profiles.grid_ggd[grid_ggd_idx]
    space = grid_ggd.space[space_idx]
    eqt = dd.equilibrium.time_slice[eq_time_idx]

    r_eq, z_eq, psin_eq, rzpi = prep_flux_map(dd; eq_time_idx, eq_profiles_2d_idx)
    psin_levels = pick_extension_psi_range(
        dd;
        eq_time_idx,
        eq_profiles_2d_idx,
        grid_ggd_idx,
        space_idx,
    )
    nlvl = length(psin_levels)
    dpsin = psin_levels[2] - psin_levels[1]
    grad_start_r, grad_start_z = pick_mesh_ext_starting_points(grid_ggd, space)
    npol = length(grad_start_r)
    mesh_r, mesh_z = mesh_ext_follow_grad(
        r_eq,
        z_eq,
        psin_eq,
        grad_start_r,
        grad_start_z,
        nlvl,
        dpsin,
        rzpi,
    )
    modify_mesh_ext_near_x!(eqt, mesh_r, mesh_z)

    # Now we have all the nodes needed for the new mesh, but none are connected
    # yet. They are, however, organized nicely in order, so it shouldn't be too
    # hard. Any X-points in the domain should be aggressively ignored, so all the
    # connections should be nice and simple.
    # The PFR flag needs to be respected; pfr nodes don't connect to non-pfr nodes
    pfr = rzpi.(mesh_r[:, 1], mesh_z[:, 1]) .< 1
    pfr_transition = argmax(abs.(diff(pfr)))
    record_regular_mesh!(grid_ggd, space, mesh_r, mesh_z, pfr_transition)
    return mesh_r, mesh_z, pfr_transition
end

"""
    fill_in_extrapolated_edge_profile!(
        dd::IMAS.dd,
        quantity_name::String;
        method::String="simple",
        eq_time_idx::Int=1,
    )

JUST A PLACEHOLDER FOR NOW. DOESN'T ACTUALLY WORK YET.
"""
function fill_in_extrapolated_edge_profile!(
    dd::IMAS.dd,
    quantity_name::String;
    method::String="simple",
    eq_time_idx::Int=1,
)
    # Inspect input
    if length(dd.equilibrium.time_slice) < eq_time_idx
        throw(
            ArgumentError(
                string(
                    "DD equilibrium does not have enough time slices: ",
                    length(dd.equilibrium.time_slice),
                    " slices were present but slice index ",
                    eq_time_idx, " was requested.",
                ),
            ),
        )
    end
    bad_ggd = (length(dd.edge_profiles.grid_ggd) < 1)
    if bad_ggd
        throw(ArgumentError(string(
            "Invalid GGD data.",
        )))
    end

    if method == "simple"
        println("THERE IS NOT ACTUALLY AN EDGE EXTRAPOLATION METHOD YET.")
        println("THIS IS A PLACEHOLDER SO THE REST OF THE WORKFLOW CAN BE SET UP.")
        # In this case, there's a good chance that we'll need a few different classes
        # of extrapolation techniques for different quantities as the may behave
        # differently.
    else
        throw(ArgumentError(string(
            "Unrecognized extraplation method: ", method,
        )))
    end
end
