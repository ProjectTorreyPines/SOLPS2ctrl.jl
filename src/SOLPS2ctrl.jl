module SOLPS2ctrl

using IMAS: IMAS
using SOLPS2imas: SOLPS2imas
using EFIT: EFIT
using Interpolations: Interpolations
using Dates
using Printf

export find_files_in_allowed_folders, geqdsk_to_imas!, preparation, write_pcs_config

include("$(@__DIR__)/supersize_profile.jl")
include("$(@__DIR__)/repair_eq.jl")
include("$(@__DIR__)/unit_utils.jl")
include("$(@__DIR__)/control_utils.jl")
include("$(@__DIR__)/plant.jl")
include("$(@__DIR__)/system_id.jl")
include("$(@__DIR__)/actuators.jl")
include("$(@__DIR__)/controllers.jl")

"""
    find_files_in_allowed_folders(
        input_dirs::String...;
        eqdsk_file::String,
        recursive::Bool=true,
    )

Searches a list of allowed folders for a set of filenames that will provide information
about the SOLPS case. Returns a list of filenames with complete paths.

Example:

```julia
SOLPS2ctrl.find_files_in_allowed_folders(
    \"<your samples folder>/D3D_Ma_184833_03600\";
    eqdsk_file=\"g184833.03600\",
)
```
"""
function find_files_in_allowed_folders(
    input_dirs::String...;
    eqdsk_file::String,
    recursive::Bool=true,
)::Vector{String}
    files = ["b2fgmtry", "b2time.nc", "b2mn.dat", eqdsk_file]
    reduced_files =
        ["b2fgmtry_red", "b2time_red.nc", "b2mn.dat", eqdsk_file]
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
    for i ∈ eachindex(files)
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
    geqdsk_to_imas!(
        eqdsk_file::String,
        dd::IMAS.dd;
        set_time::Union{Nothing, Float64}=nothing,
        time_index::Int=1,
    )

Transfers the equilibrium reconstruction from an EFIT-style gEQDSK file into
the IMAS DD structure.
"""
function geqdsk_to_imas!(
    eqdsk_file::String,
    dd::IMAS.dd;
    set_time::Union{Nothing, Float64}=nothing,
    time_index::Int=1,
    allow_boundary_flux_correction::Bool=false,
)
    # https://github.com/JuliaFusion/EFIT.jl/blob/master/src/io.jl
    g = EFIT.readg(eqdsk_file; set_time=set_time)
    gfilename = split(eqdsk_file, "/")[end]
    # Copying ideas from OMFIT: omfit/omfit_classes/omfit_eqdsk.py / to_omas()
    eq = dd.equilibrium
    if IMAS.ismissing(eq, :time)
        eq.time = Array{Float64}(undef, time_index)
    end
    eq.time[time_index] = g.time
    if length(eq.time_slice) < time_index
        resize!(eq.time_slice, time_index)
    end
    eqt = eq.time_slice[time_index]
    eqt.time = g.time

    source_for_summary = "gEQDSK file $gfilename loaded during SOLPS2ctrl workflow."

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

    if IMAS.ismissing(dd.summary, :time)
        dd.summary.time = Array{Float64}(undef, time_index)
    end
    dd.summary.time[time_index] = g.time
    ip = Array{Float64}(undef, time_index)
    ip[time_index] = g.current
    dd.summary.global_quantities.ip.value = ip
    dd.summary.global_quantities.r0.value = g.rcentr
    dd.summary.global_quantities.b0.value = b0
    summarize = ["ip", "r0", "b0"]

    # 2D
    resize!(eqt.profiles_2d, 1)
    p2 = eqt.profiles_2d[1]
    p2.grid.dim1 = collect(g.r)
    p2.grid.dim2 = collect(g.z)
    p2.psi = g.psirz  # Not sure if transpose is correct (I have been getting away with this for some time and suspect it's okay)
    p2.grid_type.index = 1  # 1 = rectangular, such as dim1 = R, dim2 = Z
    p2.grid_type.name = "R-Z grid for flux map"
    p2.grid_type.description = (
        "A rectangular grid of points in R,Z on which poloidal " *
        "magnetic flux psi is defined. The grid's dim1 is R, dim2 is Z."
    )
    # missing j_tor = pcurrt

    if allow_boundary_flux_correction
        # Check / correct simag. Intended for the case where simag is quoted imprecisely and the contour doesn't close
        level = gq.psi_boundary + 0
        paths, level = IMAS.flux_surface(eqt, level, :closed)
        count = 0  # prevents inf loop if something goes wrong
        while (length(paths) >= 1) && (count < 50)  # push boundary flux out until there is no closed boundary
            level += (gq.psi_boundary - gq.psi_axis) * 0.001
            paths, level = IMAS.flux_surface(eqt, level, :closed)
            count += 1
        end
        println(" --- ")
        while (length(paths) < 1) && (count < 200)  # pull boundary flux back in until there is a closed boundary
            level -= (gq.psi_boundary - gq.psi_axis) * 0.0001
            paths, level = IMAS.flux_surface(eqt, level, :closed)
            count += 1
        end
        println(
            "Correcting psi_boundary from $(gq.psi_boundary) to $level so contouring will find a closed flux surface.",
        )
        gq.psi_boundary = level
    end

    # 1D
    p1 = eqt.profiles_1d
    p1.psi = collect(g.psi)
    p1.f = g.fpol
    p1.pressure = g.pres
    p1.f_df_dpsi = g.ffprim
    p1.dpressure_dpsi = g.pprime
    p1.q = g.qpsi
    if hasproperty(g, :rhovn)
        p1.rho_tor_norm = g.rhovn
    end

    # Derived
    psin1d = (g.psi .- gq.psi_axis) ./ (gq.psi_boundary - gq.psi_axis)
    gq.magnetic_axis.b_field_tor = g.bcentr * g.rcentr / g.rmaxis
    gq.q_axis = g.qpsi[1]
    gq.q_95 = Interpolations.linear_interpolation(psin1d, g.qpsi)(0.95)
    qmin_idx = argmin(abs.(g.qpsi))
    gq.q_min.value = g.qpsi[qmin_idx]
    if hasproperty(g, :rhovn)
        gq.q_min.rho_tor_norm = g.rhovn[qmin_idx]
    end

    dd.summary.global_quantities.q_95.value = Array{Float64}(undef, time_index)
    dd.summary.global_quantities.q_95.value[time_index] = gq.q_95
    summarize = [summarize; "q_95"]

    # X-points
    xrs, xzs, xpsins, xseps = EFIT.x_points(g; within_limiter_only=false)
    if length(xrs) > 0
        bx = eqt.boundary.x_point
        resize!(bx, length(xrs))
        for i ∈ eachindex(xrs)
            bx[i].r = xrs[i]
            bx[i].z = xzs[i]
        end
        nprim = sum(xseps .== 1)
        if nprim > 0
            bsx = eqt.boundary_separatrix.x_point
            resize!(bsx, nprim)
            xrprim = xrs[xseps .== 1]
            xzprim = xzs[xseps .== 1]
            for i ∈ nprim
                bsx[i].r = xrprim[i]
                bsx[i].z = xzprim[i]
            end
        end
        nsec = sum(xseps .== 2)
        if nsec > 0
            bssx = eqt.boundary_secondary_separatrix.x_point
            resize!(bssx, nsec)
            xrsec = xrs[xseps .== 2]
            xzsec = xzs[xseps .== 2]
            for i ∈ nsec
                bssx[i].r = xrsec[i]
                bssx[i].z = xzsec[i]
            end
        end
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
    limiter.unit[1].outline.z = g.zlim

    for s ∈ summarize
        getproperty(dd.summary.global_quantities, Symbol(s)).source = source_for_summary
    end

    return
end

"""
    preparation(
        eqdsk_file::String,
        dirs::String...;
        core_method::String="simple",
        filename::String="sd_input_data",
        output_format::String="json",
        eqdsk_set_time::Union{Nothing, Float64}=nothing,
        eq_time_index::Int=1,
        allow_boundary_flux_correction::Bool=false,
    )::IMAS.dd

Gathers SOLPS and EFIT files and loads them into IMAS structure. Extrapolates
profiles as needed to get a complete picture.
"""
function preparation(
    eqdsk_file::String,
    dirs::Vector{String};
    core_method::String="simple",
    filename::String="sd_input_data",
    output_format::String="json",
    eqdsk_set_time::Union{Nothing, Float64}=nothing,
    eq_time_index::Int=1,
    allow_boundary_flux_correction::Bool=false,
)::IMAS.dd
    b2fgmtry, b2time, b2mn, eqdsk =
        find_files_in_allowed_folders(dirs...; eqdsk_file=eqdsk_file)
    println("Found source files:")
    println("    b2fgmtry = ", b2fgmtry)
    println("    b2time = ", b2time)
    println("    b2mn.dat = ", b2mn)
    println("    eqdsk = ", eqdsk)

    dd = IMAS.dd()
    geqdsk_to_imas!(
        eqdsk,
        dd;
        set_time=eqdsk_set_time,
        time_index=eq_time_index,
        allow_boundary_flux_correction=allow_boundary_flux_correction,
    )
    # Fill out more equilibrium data
    add_rho_to_equilibrium!(dd)  # Doesn't do anything if rho is valid
    dd.global_time = dd.equilibrium.time_slice[1].time
    for eqt ∈ dd.equilibrium.time_slice
        IMAS.flux_surfaces(eqt, dd.wall)
    end
    # Add SOLPS data
    dd = SOLPS2imas.solps2imas(b2fgmtry, b2time; b2mn=b2mn, ids=dd)
    println("Loaded input data into IMAS DD")

    # Core profiles
    # Set timing
    nt = length(dd.edge_profiles.ggd)
    if length(dd.core_profiles.profiles_1d) < nt
        resize!(dd.core_profiles.profiles_1d, nt)
    end
    if ismissing(dd.core_profiles, :time)
        dd.core_profiles.time = Array{Float64}(undef, nt)
    elseif length(dd.core_profiles.time) < nt
        resize!(dd.core_profiles.time, nt)
    end
    for it ∈ 1:nt
        dd.core_profiles.time[it] =
            dd.core_profiles.profiles_1d[it].time = dd.edge_profiles.ggd[it].time
    end
    # Extrapolate profiles
    core_profiles = ["electrons.density", "electrons.temperature"]
    extrapolated_core_profiles = []
    for core_profile ∈ core_profiles
        tags = split(core_profile, ".")
        parent = dd.edge_profiles.ggd[1]
        for tag ∈ tags[1:(end-1)]
            parent = getproperty(parent, Symbol(tag))
        end
        qty = getproperty(parent, Symbol(tags[end]), core_profile)
        if length(qty) > 0
            fill_in_extrapolated_core_profile!(dd, core_profile; method=core_method)
            append!(extrapolated_core_profiles, [core_profile])
            println("    Extrapolated $core_profile into the core.")
        else
            println(
                "  > Warning: quantity $core_profile was not usable and was not extrapolated into the core.",
            )
        end
    end
    # ... more profiles here as they become available in b2time
    println(
        "Extrapolated $(length(extrapolated_core_profiles))/$(length(core_profiles)) core profiles.",
    )

    cached_mesh_extension!(dd, eqdsk_file, b2fgmtry)
    fill_in_extrapolated_edge_profile!(dd, "electrons.density"; method=core_method)
    # ... more profiles here
    println("Extrapolated edge profiles (but not really (placeholder only))")

    print("Exporting to file: ")
    if output_format == "json"
        println(filename * ".json")
        IMAS.imas2json(dd, filename * ".json"; strict=true, freeze=false)
    else
        throw(ArgumentError(string("Unrecognized output format: ", output_format)))
    end
    return dd
end

function zero_pad_str(array::Vector, len::Int64, format=Printf.Format("%+13.6e"))
    padded = zeros(len)
    m = min(length(array), len)
    for i ∈ 1:m
        padded[i] = array[i]
    end
    str_rep = join([Printf.format(format, pad) for pad ∈ padded], " ")
    return str_rep
end

function zero_pad_str(
    array::Matrix,
    rows::Int64,
    cols::Int64,
    format=Printf.Format("%+13.6e"),
)
    padded = zeros(rows, cols)
    m = min(size(array)[1], rows)
    n = min(size(array)[2], cols)
    for i ∈ 1:m
        for j ∈ 1:n
            padded[i, j] = array[i, j]
        end
    end
    row_strings = ""
    for i ∈ 1:rows
        str_r = join([Printf.format(format, pad) for pad ∈ padded[i, :]], " ")
        row_strings = join([row_strings, str_r], "\n")
    end
    return row_strings
end

function username()
    possible_environment_vars_for_username = ["USER", "USERNAME", "LNAME", "LOGNAME"]
    for varname ∈ possible_environment_vars_for_username
        haskey(ENV, varname) && return ENV[varname]
    end
    return nothing
end

user = username()

function write_pcs_config(
    # The model
    A::Matrix,
    B::Matrix,
    C::Matrix,
    D::Matrix,
    Y2x::Matrix,
    U2x::Matrix,
    n_inputs::Int64,
    n_outputs::Int64,
    n_states::Int64,
    n_history::Int64,
    model_type_code::Int64,  # 0 = density, 1 = heat flux of div Te (direct engineering limit exhaust control), 2 = is detachment or dissipation like Afrac or frad (indirect exhaust control), 3 = other
    gas_atomic_number::Int64,
    model_time_step_sec::Float64,
    # Adaptation instructions
    adapt_input_Gp::Vector{Float64},  # units are inverse of input units
    adapt_input_tauI::Vector{Float64},  # seconds
    adapt_output_Gp::Vector{Float64},  # units are inverse of output units
    adapt_output_tauI::Vector{Float64},  # seconds
    adapt_factor_min::Float64,  # unitless
    adapt_factor_max::Float64,
    # Input description
    input_offset::Vector{Float64},
    input_factor::Vector{Float64},  # normalized input = (input - offset) * factor, one entry per input
    input_description::String,  # Describe things like "input 1 is deuterium, input 2 is neon, ..."; this should allow operator to connect PCS input signals to the model
    # Output description
    output_offset::Vector{Float64},
    output_factor::Vector{Float64},
    output_description::String,  # Describe things like "output 1 is density, output 2 is peak heatflux on outer divertor, ..." so that operator can connect PCS measuerments to the model
    # Controller
    Gp::Vector{Float64},  # One per actively controlled input, units are gas flow units divided by output units
    tauI::Vector{Float64}, # One per input, units are seconds
    # Preparation records
    model_description::String="No description given",
    model_counter::Int64=0,  # increment this if you make multiple similar models on the same day and need to differentiate them
    prepared_by::String=user,  # Your name
    secondary_gas_atomic_number::Int64=0,
)
    # PCS settings
    DENSITYCS_N_PREDICTOR_INPUTS = 5
    DENSITYCS_N_PREDICTOR_OUTPUTS = 4
    DENSITYCS_N_PREDICTOR_STATES = 6
    DENSITYCS_N_PREDICTOR_HISTORY = 200
    DENSITYCS_N_JMATRIX_COLS = 800
    DENSITYCS_N_KMATRIX_COLS = 1000

    date = Dates.format(Date(today()), "yyyymmdd")

    serial =
        date * string(n_states) * string(n_inputs) * string(n_outputs) *
        string(n_history; pad=4) * string(model_type_code) *
        string(gas_atomic_number; pad=2) * string(secondary_gas_atomic_number; pad=2) *
        string(model_counter; pad=2)
    filename = "pvlc_config_" * serial * ".txt"
    ff = Printf.Format("%+13.6e")

    open(filename, "w") do io
        println(io, "PVLC config")
        println(io, serial)
        println(io, "Model description:   " * model_description)
        println(io, "Generated by:        " * prepared_by)
        println(io, "Generated date:      " * date)
        println(io, "Counter:             $(model_counter)")
        println(io, "Primary gas atomic #:$(gas_atomic_number)")
        println(io, "Secondary atomic #:  $(secondary_gas_atomic_number)")
        println(io, "Input instructions:  " * input_description)
        println(io, "Output instructions: " * output_description)
        println(io, "Particle counts (densities, flows, etc.) are divided by 1e19.")
        println(
            io,
            "Prepared with the following assumptions about PCS setup (these are hard coded into the PCS and they must match):",
        )
        println(io, "::")
        println(io, "DENSITYCS_N_PREDICTOR_INPUTS  = $(DENSITYCS_N_PREDICTOR_INPUTS)")
        println(io, "DENSITYCS_N_PREDICTOR_OUTPUTS = $(DENSITYCS_N_PREDICTOR_OUTPUTS)")
        println(io, "DENSITYCS_N_PREDICTOR_STATES  = $(DENSITYCS_N_PREDICTOR_STATES)")
        println(io, "DENSITYCS_N_PREDICTOR_HISTORY = $(DENSITYCS_N_PREDICTOR_HISTORY)")
        println(io, "DENSITYCS_N_JMATRIX_COLS      = $(DENSITYCS_N_JMATRIX_COLS)")
        println(io, "DENSITYCS_N_KMATRIX_COLS      = $(DENSITYCS_N_KMATRIX_COLS)")
        println(io, "::")
        println(io, "Number of states in the model      : $(n_states)")
        println(io, "Number of inputs                   : $(n_inputs)")
        println(io, "Number of outputs                  : $(n_outputs)")
        println(io, "Number of history steps (tracking) : $(n_history)")
        println(
            io,
            "Model step time (seconds)          : " *
            Printf.format(ff, model_time_step_sec),
        )
        println(
            io,
            "G_P (1e19 molec/s per physics unit): " *
            zero_pad_str(Gp, DENSITYCS_N_PREDICTOR_INPUTS),
        )
        println(
            io,
            "Tau_I (seconds)                    : " *
            zero_pad_str(tauI, DENSITYCS_N_PREDICTOR_INPUTS),
        )
        println(
            io,
            "Input offsets (10^19 molecules/sec): " *
            zero_pad_str(input_offset, DENSITYCS_N_PREDICTOR_INPUTS),
        )
        println(
            io,
            "Input norm factors                 : " *
            zero_pad_str(input_factor, DENSITYCS_N_PREDICTOR_INPUTS),
        )
        println(
            io,
            "Output offsets (physics units)     : " *
            zero_pad_str(output_offset, DENSITYCS_N_PREDICTOR_OUTPUTS),
        )
        println(
            io,
            "Output norm factors                : " *
            zero_pad_str(output_factor, DENSITYCS_N_PREDICTOR_OUTPUTS),
        )
        println(
            io,
            "Input adaptation G_P (1/in.units)  : " *
            zero_pad_str(adapt_input_Gp, DENSITYCS_N_PREDICTOR_OUTPUTS),
        )
        println(
            io,
            "Input adaptation tau_I (seconds)   : " *
            zero_pad_str(adapt_input_tauI, DENSITYCS_N_PREDICTOR_OUTPUTS),
        )
        println(
            io,
            "Output adaptation G_P (1/out.units): " *
            zero_pad_str(adapt_output_Gp, DENSITYCS_N_PREDICTOR_OUTPUTS),
        )
        println(
            io,
            "Output adaptation tau_I (seconds)  : " *
            zero_pad_str(adapt_output_tauI, DENSITYCS_N_PREDICTOR_OUTPUTS),
        )
        println(
            io,
            "Adaptation factor minimum          : " *
            Printf.format(ff, adapt_factor_min),
        ) #* zero_pad_str(adapt_factor_min, DENSITYCS_N_PREDICTOR_OUTPUTS))
        println(
            io,
            "Adaptation factor maximum          : " *
            Printf.format(ff, adapt_factor_max),
        ) #* zero_pad_str(adapt_factor_max, DENSITYCS_N_PREDICTOR_OUTPUTS))
        println(io, "::")
        println(
            io,
            "A = " *
            zero_pad_str(A, DENSITYCS_N_PREDICTOR_STATES, DENSITYCS_N_PREDICTOR_STATES),
        )
        println(
            io,
            "B = " *
            zero_pad_str(B, DENSITYCS_N_PREDICTOR_STATES, DENSITYCS_N_PREDICTOR_INPUTS),
        )
        println(
            io,
            "C = " * zero_pad_str(
                C,
                DENSITYCS_N_PREDICTOR_OUTPUTS,
                DENSITYCS_N_PREDICTOR_STATES,
            ),
        )
        println(
            io,
            "D = " * zero_pad_str(
                D,
                DENSITYCS_N_PREDICTOR_OUTPUTS,
                DENSITYCS_N_PREDICTOR_INPUTS,
            ),
        )
        println(
            io,
            "Y2x (J) = " *
            zero_pad_str(Y2x, DENSITYCS_N_PREDICTOR_STATES, DENSITYCS_N_JMATRIX_COLS),
        )
        return println(
            io,
            "U2x (K) = " *
            zero_pad_str(U2x, DENSITYCS_N_PREDICTOR_STATES, DENSITYCS_N_KMATRIX_COLS),
        )
    end
end

end # module SOLPS2ctrl
