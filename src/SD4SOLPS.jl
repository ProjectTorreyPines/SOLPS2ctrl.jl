module SD4SOLPS

using OMAS: OMAS
using SOLPS2IMAS: SOLPS2IMAS
using EFIT: EFIT
using Interpolations: Interpolations
#import GGDUtils

export find_files_in_allowed_folders
export geqdsk_to_imas

include("$(@__DIR__)/supersize_profile.jl")
include("$(@__DIR__)/repair_eq.jl")
include("$(@__DIR__)/actuator_model.jl")

greet() = print("Hello World!")

"""
    find_files_in_allowed_folders()

Searches a list of allowed folders for a set of filenames that will provide information
about the SOLPS case. Returns a list of filenames with complete paths.

Example:
SD4SOLPS.find_files_in_allowed_folders(
"<your samples folder>/D3D_Ma_184833_03600",
eqdsk_file="g184833.03600",
)
"""
function find_files_in_allowed_folders(
    input_dirs...;
    eqdsk_file,
    recursive=true,
    allow_reduced_versions=false,
)
    files = ["b2fgmtry", "b2time.nc", "b2mn.dat", "gridspacedesc.yml", eqdsk_file]
    reduced_files =
        ["b2fgmtry_red", "b2time_red.nc", "b2mn.dat", "gridspacedesc.yml", eqdsk_file]
    output_files = fill("", length(files))
    if recursive
        dirs = []
        for dir ∈ input_dirs
            dirs =
                append!(dirs, [subdir[1] for subdir ∈ [item for item ∈ walkdir(dir)]])
        end
    else
        dirs = input_dirs
    end
    for i ∈ 1:length(files)
        for dir ∈ dirs
            file = dir * "/" * files[i]
            reduced_file = dir * "/" * reduced_files[i]
            if isfile(file)
                output_files[i] = file
                break
            elseif isfile(reduced_file)
                output_files[i] = reduced_file
                break
            end
        end
    end
    return output_files
end

"""
    geqdsk_to_imas()

Transfers the equilibrium reconstruction in an EFIT-style gEQDSK file into
the IMAS DD structure.
"""
function geqdsk_to_imas(eqdsk_file, dd; time_index=1)
    # https://github.com/JuliaFusion/EFIT.jl/blob/master/src/io.jl
    g = EFIT.readg(eqdsk_file)
    # Copying ideas from OMFIT: omfit/omfit_classes/omfit_eqdsk.py / to_omas()
    eq = dd.equilibrium
    resize!(eq.time_slice, 1)
    eqt = eq.time_slice[time_index]

    # 0D
    gq = eqt.global_quantities
    gq.magnetic_axis.r = g.rmaxis
    gq.magnetic_axis.z = g.zmaxis
    gq.psi_axis = g.simag
    gq.psi_boundary = g.sibry
    gq.ip = g.current
    eq.vacuum_toroidal_field.r0 = g.rcentr
    b0 = Array{Float64}(undef, time_index)
    b0[time_index] = g.bcentr
    eq.vacuum_toroidal_field.b0 = b0

    # 1D
    p1 = eqt.profiles_1d
    nprof = length(g.pres)
    psi = collect(LinRange(g.simag, g.sibry, nprof))
    p1.psi = psi
    p1.f = g.fpol
    p1.pressure = g.pres
    p1.f_df_dpsi = g.ffprim
    p1.dpressure_dpsi = g.pprime
    p1.q = g.qpsi
    if hasproperty(g, :rhovn)
        # rhovn is not in the original EFIT.jl but is added on a branch
        p1.rho_tor_norm = g.rhovn
    end

    # 2D
    resize!(eqt.profiles_2d, 1)
    p2 = eqt.profiles_2d[1]
    p2.grid.dim1 = collect(g.r)
    p2.grid.dim2 = collect(g.z)
    p2.psi = g.psirz  # Not sure if transpose is correct
    # missing j_tor = pcurrt

    # Derived
    psin1d = (psi .- g.simag) ./ (g.sibry - g.simag)
    gq.magnetic_axis.b_field_tor = g.bcentr * g.rcentr / g.rmaxis
    gq.q_axis = g.qpsi[1]
    gq.q_95 = Interpolations.LinearInterpolation(psin1d, g.qpsi)(0.95)
    qmin_idx = argmin(abs.(g.qpsi))
    gq.q_min.value = g.qpsi[qmin_idx]
    if hasproperty(g, :rhovn)
        gq.q_min.rho_tor_norm = g.rhovn[qmin_idx]
    end

    # Boundary / LCFS
    eqt.boundary.outline.r = g.rbbbs
    eqt.boundary.outline.z = g.zbbbs

    # Wall
    resize!(dd.wall.description_2d, 1)
    limiter = dd.wall.description_2d[1].limiter
    limiter.type.name = "first wall"
    limiter.type.index = 0
    limiter.type.description = "first wall"
    resize!(limiter.unit, 1)
    limiter.unit[1].outline.r = g.rlim
    return limiter.unit[1].outline.z = g.zlim
end

"""
    core_profile_2d(dd, prof_time_idx, eq_time_idx, quantity, r, z)

Reads a 1D core profile and flux map and returns a quantity at requested R,Z points
dd: a data dictionary instance with equilibrium and core profile data loaded
quantity: a string specifying the quantity to fetch
r: R locations of desired output points / m
z: Z locations of desired output points / m
"""
function core_profile_2d(dd, prof_time_idx, eq_time_idx, quantity, r, z)
    if !check_rho_1d(dd; time_slice=eq_time_idx)
        throw(ArgumentError("Equilibrium rho profile in input DD was missing."))
    end
    prof = dd.core_profiles.profiles_1d[prof_time_idx]
    rho_prof = prof.grid.rho_tor_norm
    quantity_fields = split(quantity, ".")
    p = prof
    for field ∈ quantity_fields
        p = getproperty(p, Symbol(field))
    end
    eqt = dd.equilibrium.time_slice[eq_time_idx]
    p1 = eqt.profiles_1d
    p2 = eqt.profiles_2d[1]
    gq = eqt.global_quantities
    psi_a = gq.psi_axis
    psi_b = gq.psi_boundary
    rhon_eq = p1.rho_tor_norm
    psi_eq = p1.psi
    psin_eq = (psi_eq .- psi_a) ./ (psi_b - psi_a)
    psirz = p2.psi
    psinrz = (psirz .- psi_a) ./ (psi_b - psi_a)
    r_eq = p2.grid.dim1
    z_eq = p2.grid.dim2
    extension = [1.0001, 1.1, 5]
    # rho_N isn't defined on open flux surfaces, so it is extended by copying psi_N
    psin_eq_ext = copy(psin_eq)
    append!(psin_eq_ext, extension)
    rhon_eq_ext = copy(rhon_eq)
    append!(rhon_eq_ext, extension)
    neg_extension = [-5, -0.0001]  # I guess this would be at a PF coil or something?
    prepend!(psin_eq_ext, neg_extension)
    prepend!(rhon_eq_ext, neg_extension)

    rho_prof_ext = vcat(rho_prof, extension)
    p_ext = vcat(p, zeros(size(extension)))
    rho_prof_ext = prepend!(rho_prof_ext, neg_extension)
    p_ext = prepend!(p_ext, zeros(size(neg_extension)))

    # Data are set up and ready to process

    # EDGE PROFILES (the input data):
    # rho_prof_ext: rho_N values associated with the profile of some quantity
    # p_ext : profile of some quantity vs. rho_prof

    # EQUILBRIUM (the connection from rho to R,Z via psi):
    # psin_eq_ext : 1D array of psi_N values that corresponds to rhon_eq_ext
    # rhon_eq_ext : 1D array of rho_N values that corresponds to psin_eq_ext
    #       --> connects rho to everything else via psi
    # psinrz : 2D array of psi_N values that corresponds to r_eq and z_eq
    # r_eq and z_eq : 1D coordinate arrays that correspond to psinrz

    # OUTPUT INSTRUCTIONS:
    # r and z : coordinates of output points where values of p are desired

    psi_at_requested_points =
        Interpolations.LinearInterpolation((r_eq, z_eq), psinrz).(r, z)
    rhonpsi = Interpolations.LinearInterpolation(psin_eq_ext, rhon_eq_ext)
    rho_at_requested_points = rhonpsi.(psi_at_requested_points)
    itp = Interpolations.LinearInterpolation(rho_prof, p)
    p_at_requested_points = itp.(rho_at_requested_points)
    return p_at_requested_points
end

"""
    preparation()

Gathers SOLPS and EFIT files and loads them into IMAS structure. Extrapolates
profiles as needed to get a complete picture.
"""
function preparation(
    eqdsk_file,
    dirs...;
    core_method::String="simple",
    edge_method::String="simple",
    filename::String="sd_input_data",
    output_format::String="json",
)
    b2fgmtry, b2time, b2mn, gridspec, eqdsk =
        find_files_in_allowed_folders(dirs...; eqdsk_file=eqdsk_file)
    println("Found source files:")
    println("    b2fgmtry = ", b2fgmtry)
    println("    b2time = ", b2time)
    println("    b2mn.dat = ", b2mn)
    println("    gridspec = ", gridspec)
    println("    eqdsk = ", eqdsk)

    dd = SOLPS2IMAS.solps2imas(b2fgmtry, b2time, gridspec, b2mn)
    geqdsk_to_imas(eqdsk, dd)
    println("Loaded input data into IMAS DD")

    fill_in_extrapolated_core_profile!(dd, "electrons.density"; method=core_method)
    fill_in_extrapolated_core_profile!(dd, "electrons.temperature"; method=core_method)
    # ... more profiles here as they become available in b2time
    println("Extrapolated core profiles")

    fill_in_extrapolated_edge_profile!(dd, "electrons.density"; method=core_method)
    # ... more profiles here
    println("Extrapolated edge profiles (but not really (placeholder only))")

    if output_format == "json"
        OMAS.imas2json(dd, filename * ".json")
    else
        throw(ArgumentError(string("Unrecognized output format: ", output_format)))
    end
    return dd
end

end # module SD4SOLPS
