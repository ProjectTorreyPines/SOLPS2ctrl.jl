using SOLPS2ctrl
using SOLPS2ctrl:
    SOLPS2ctrl, SOLPS2imas, IMAS, IMASggd, EFIT, Unitful, Interpolations, Statistics
using Plots
using Test
using DelimitedFiles: readdlm
using ControlSystemsBase: lsim, ssrand, delay, tf, c2d, ss, pid, margin

"""
    make_test_profile()

Makes a modified tanh profile based on the formula from:
Groebner, et al., Nucl. Fusion 41, 1789 (2001) 10.1088/0029-5515/41/12/306

x: position basis. The output will be a function of this coordinate.
Could be psi_N or R or Z or whatever, as long as sym and hwid are
specified in terms of the same coordinate.
Defaults are given in terms of psi_N.
sym: location of the tanh symmetry point
hwid: half width of the tanh
offset: value at the foot of the tanh / SOL value
pedestal: value at the top of the tanh / pedestal value
alpha: interior slope
"""
function make_test_profile(
    x;
    sym=0.9657,
    hwid=0.05421,
    offset=0.0953,
    pedestal=3.58901,
    alpha=0.005,
)
    z = (sym .- x) ./ hwid
    amplitude = (pedestal - offset) / 2.0
    b = offset + amplitude
    return b .+
           amplitude .* ((1 .+ alpha .* z) .* exp.(z) .- exp.(.-z)) ./
           (exp.(z) .+ exp.(.-z))
end

# function plot_test_profiles()
#     edge_rho = Array(LinRange(0.88, 1.0, 18))
#     edge_quantity = make_test_profile(edge_rho)
#     output_rho = Array(LinRange(0, 1.0, 201))
#     output_quantity = SOLPS2ctrl.extrapolate_core(edge_rho, edge_quantity, output_rho)
#     Plots.plot(output_rho, output_quantity)
#     Plots.plot!(edge_rho, edge_quantity, marker='.')
# end

function define_default_sample_set()
    sample_path =
        splitdir(pathof(SOLPS2ctrl))[1] *
        "/../sample/ITER_Lore_2296_00000/extended_output"
    sample_path2 =
        splitdir(pathof(SOLPS2ctrl))[1] * "/../sample/ITER_Lore_2296_00000/baserun"
    sample_path3 =
        splitdir(pathof(SOLPS2ctrl))[1] * "/../sample/ITER_Lore_2296_00000/run_restart"

    # Requires dvc pull before the full samples will be loaded
    # Sorry, the minimal samples are too minimal for this one.
    file_list = SOLPS2ctrl.find_files_in_allowed_folders(
        sample_path, sample_path2, sample_path3;
        eqdsk_file="thereisntoneyet",
    )
    b2fgmtry, b2time, b2mn, eqdsk = file_list
    eqdsk =
        splitdir(pathof(SOLPS2ctrl))[1] *
        "/../sample/ITER_Lore_2296_00000/EQDSK/g002296.00200"
    return b2fgmtry, b2time, b2mn, eqdsk
end

if isempty(ARGS) || "units" in ARGS
    @testset "Unit conversion utilities" begin
        # Gas unit converter
        flow_tls = 40.63 * Unitful.Torr * Unitful.L / Unitful.s
        flow_pam3 =
            SOLPS2ctrl.gas_unit_converter(flow_tls, "torr L s^-1", "Pa m^3 s^-1")
        flow_pam3_no_unitful =
            SOLPS2ctrl.gas_unit_converter(flow_tls.val, "torr L s^-1", "Pa m^3 s^-1")
        @test flow_pam3.val > 0.0
        @test flow_pam3_no_unitful ==
              Unitful.uconvert(Unitful.Pa * Unitful.m^3 / Unitful.s, flow_pam3).val

        flow_molecules1 = SOLPS2ctrl.gas_unit_converter(
            flow_tls,
            "torr L s^-1",
            "molecules s^-1";
            temperature=293.15 * Unitful.K,
        )
        flow_molecules2 = SOLPS2ctrl.gas_unit_converter(
            flow_tls,
            "torr L s^-1",
            "molecules s^-1";
            temperature=300.0 * Unitful.K,
        )
        flow_molecules3 = SOLPS2ctrl.gas_unit_converter(
            flow_tls.val,
            "torr L s^-1",
            "molecules s^-1";
            temperature=300.0,
        )
        @test flow_molecules1 > flow_molecules2
        @test (Unitful.upreferred(flow_molecules2)).val == flow_molecules3
    end
end

if isempty(ARGS) || "core" in ARGS
    @testset "core_profile_extension" begin
        # Just the basic profile extrapolator ------------------
        edge_rho = Array(LinRange(0.88, 1.0, 18))
        edge_quantity = make_test_profile(edge_rho)
        output_rho = Array(LinRange(0, 1.0, 201))
        output_quantity =
            SOLPS2ctrl.extrapolate_core(edge_rho, edge_quantity, output_rho)
        @test length(output_quantity) == length(output_rho)

        # The full workflow --------------------------------------
        # Setup sample DD
        b2fgmtry, b2time, b2mn, eqdsk = define_default_sample_set()
        println(b2fgmtry)
        println(b2time)
        println(b2mn)
        println(eqdsk)
        # If these files don't exist, complete the DVC sample setup and try again
        @test isfile(b2fgmtry)
        @test isfile(b2time)
        @test isfile(b2mn)
        @test isfile(eqdsk)
        eqdsk_time = parse(Float64, split(eqdsk, ".")[end]) / 1000.0
        dd = SOLPS2imas.solps2imas(b2fgmtry, b2time; b2mn=b2mn)
        SOLPS2ctrl.geqdsk_to_imas!(eqdsk, dd; set_time=eqdsk_time)
        rho = dd.equilibrium.time_slice[1].profiles_1d.rho_tor_norm

        if !SOLPS2ctrl.check_rho_1d(dd; time_slice=1)
            SOLPS2ctrl.add_rho_to_equilibrium!(dd)
            rho = dd.equilibrium.time_slice[1].profiles_1d.rho_tor_norm
            println("Repaired missing rho for core profile test")
        end
        # Sample is ready

        # Test settings
        quantity_name = "electrons.density"
        test_slice_idx = 1

        # Do it
        SOLPS2ctrl.fill_in_extrapolated_core_profile!(dd, quantity_name)

        # Inspect results
        @test length(dd.core_profiles.profiles_1d) > 0
        core_prof = dd.core_profiles.profiles_1d[test_slice_idx]
        tags = split(quantity_name, ".")
        quantity = dd.core_profiles.profiles_1d[test_slice_idx]
        for tag ∈ tags
            quantity = getproperty(quantity, Symbol(tag))
        end
        rho_core = dd.core_profiles.profiles_1d[test_slice_idx].grid.rho_tor_norm
        @test length(quantity) > 0
        @test length(quantity) == length(rho_core)
    end
end

if isempty(ARGS) || "edge" in ARGS
    @testset "edge_profile_extension" begin
        # Test for getting mesh spacing
        b2fgmtry, b2time, b2mn, eqdsk = define_default_sample_set()
        eqdsk_time = parse(Float64, split(eqdsk, ".")[end]) / 1000.0
        dd = SOLPS2imas.solps2imas(b2fgmtry, b2time; b2mn=b2mn)
        SOLPS2ctrl.geqdsk_to_imas!(eqdsk, dd; set_time=eqdsk_time)
        dpsin = SOLPS2ctrl.mesh_psi_spacing(dd)
        @test dpsin > 0.0

        # Extend the mesh
        grid_ggd_idx = 1
        grid_ggd = dd.edge_profiles.grid_ggd[grid_ggd_idx]
        extended_subs = 1:5
        orig_subs = [
            deepcopy(IMASggd.get_grid_subset(grid_ggd, i)) for
            i ∈ extended_subs
        ]
        cfn = SOLPS2ctrl.cached_mesh_extension!(dd, eqdsk, b2fgmtry; clear_cache=true)
        println("cleared ext mesh cache: ", cfn)
        SOLPS2ctrl.cached_mesh_extension!(
            dd,
            eqdsk,
            b2fgmtry;
            grid_ggd_idx=grid_ggd_idx,
        )
        for j ∈ extended_subs
            orig_sub = orig_subs[j]
            std_sub = IMASggd.get_grid_subset(grid_ggd, -j)
            all_sub = IMASggd.get_grid_subset(grid_ggd, j)
            ext_sub = IMASggd.get_grid_subset(grid_ggd, -200 - j)
            orig_indices = [ele.object[1].index for ele ∈ orig_sub.element]
            std_indices = [ele.object[1].index for ele ∈ std_sub.element]
            all_indices = [ele.object[1].index for ele ∈ all_sub.element]
            ext_indices = [ele.object[1].index for ele ∈ ext_sub.element]
            @test std_sub.identifier.index == -j
            @test all(orig_indices .== std_indices)
            all_indices_reconstruct = [std_indices; ext_indices]
            @test all(all_indices .== all_indices_reconstruct)

            # Verify that original and standard are separate and not refs to each other
            arbitrary_change = 5
            arb_el = 2
            orig_sub.element[arb_el].object[1].index += arbitrary_change
            @test orig_sub.element[arb_el].object[1].index !=
                  std_sub.element[arb_el].object[1].index
            orig_sub.element[arb_el].object[1].index -= arbitrary_change
        end

        # Prepare profiles that need to be extended
        n_edge = 47
        n_outer_prof = 13
        quantity_edge =
            [collect(LinRange(2.0, 75.0, 7)); collect(LinRange(75.0, 98.3, 17))]
        quantity_edge = [quantity_edge; reverse(quantity_edge)[2:end]]
        gradient_edge =
            [collect(LinRange(-1, -100, 7)); collect(LinRange(-100, -350, 17))]
        gradient_edge = [gradient_edge; reverse(gradient_edge)[2:end]]
        psin_out = collect(LinRange(1.0, 1.25, n_outer_prof + 1))[2:end]
    end
end

if isempty(ARGS) || "heavy" in ARGS
    @testset "heavy_utilities" begin
        # Test for finding files in allowed folders
        file_list = define_default_sample_set()
        @test length(file_list) == 4
        b2fgmtry, b2time, b2mn, eqdsk = file_list
        @test length(b2fgmtry) > 10
        @test endswith(b2fgmtry, "b2fgmtry_red") | endswith(b2fgmtry, "b2fgmtry")
        @test length(b2time) > 10
        @test endswith(b2time, "b2time_red.nc") | endswith(b2time, "b2time.nc")
        @test length(b2mn) > 10
        @test endswith(b2mn, "b2mn.dat")

        # Test for sweeping 1D core profiles into 2D R,Z
        # (or anyway evaluating them at any R,Z location)
        dd = IMAS.dd()
        eqdsk_file =
            splitdir(pathof(SOLPS2ctrl))[1] * "/../sample/geqdsk_iter_small_sample"
        SOLPS2ctrl.geqdsk_to_imas!(eqdsk_file, dd)
        quantity = "electrons.density"
        prof_time_idx = eq_time_idx = 1
        resize!(dd.core_profiles.profiles_1d, prof_time_idx)
        n = 101
        rho_n = Array(LinRange(0, 1.0, n))
        dd.core_profiles.profiles_1d[prof_time_idx].grid.rho_tor_norm = rho_n
        dd.core_profiles.profiles_1d[prof_time_idx].electrons.density =
            make_test_profile(rho_n)
        r_mag_axis = dd.equilibrium.time_slice[1].global_quantities.magnetic_axis.r
        z_mag_axis = dd.equilibrium.time_slice[1].global_quantities.magnetic_axis.z
        rg = r_mag_axis .* [0.85, 0.90, 0.95, 1.0, 1.05, 1.10, 1.15]
        zg = z_mag_axis .+ (r_mag_axis .* [-0.05, -0.025, 0, 0.025, 0.05])
        points = collect(Base.Iterators.product(rg, zg))
        r = getfield.(points, 1)
        z = getfield.(points, 2)
        if !SOLPS2ctrl.check_rho_1d(dd; time_slice=eq_time_idx)
            SOLPS2ctrl.add_rho_to_equilibrium!(dd)
            println("DD was repaired (rho added) for core 2d utility test")
        end
        density_on_grid =
            IMASggd.interp(
                dd.core_profiles.profiles_1d[1].electrons.density,
                dd.core_profiles.profiles_1d[1],
                dd.equilibrium.time_slice[1],
            ).(
                r,
                z,
            )
        @test size(density_on_grid) == (length(rg), length(zg))
    end
end

if isempty(ARGS) || "repair" in ARGS
    @testset "repair_eq" begin
        # Prepare sample
        dd = IMAS.dd()
        eqdsk = splitdir(pathof(SOLPS2ctrl))[1] * "/../sample/geqdsk_iter_small_sample"
        SOLPS2ctrl.geqdsk_to_imas!(eqdsk, dd)
        # Make sure rho is missing
        nt = length(dd.equilibrium.time_slice)
        for it ∈ 1:nt
            eqt = dd.equilibrium.time_slice[it]
            eqt.profiles_1d.rho_tor_norm = Vector{Float64}()
        end
        # Add rho
        SOLPS2ctrl.add_rho_to_equilibrium!(dd)
        # Check
        rho = dd.equilibrium.time_slice[1].profiles_1d.rho_tor_norm
        println(rho)
        @test length(rho) > 10
        @test maximum(rho) > 0
    end
end

if isempty(ARGS) || "geqdsk" in ARGS
    @testset "geqdsk_to_imas" begin
        sample_files =
            (splitdir(pathof(SOLPS2ctrl))[1] * "/../sample/") .* [
                "g184833.03600", "geqdsk_iter_small_sample",
            ]
        append!(
            sample_files,
            ["$(@__DIR__)/../sample/ITER_Lore_2296_00000/EQDSK/g002296.00200"],
        )
        tslice = 1
        for sample_file ∈ sample_files
            println(sample_file)
            dd = IMAS.dd()
            if endswith(sample_file, "00")
                eqdsk_time = parse(Float64, split(sample_file, ".")[end]) / 1000.0
            else
                eqdsk_time = nothing
            end
            SOLPS2ctrl.geqdsk_to_imas!(
                sample_file,
                dd;
                set_time=eqdsk_time,
                time_index=tslice,
            )
            eqt = dd.equilibrium.time_slice[tslice]

            # global
            gq = eqt.global_quantities
            @test gq.magnetic_axis.r > 0
            @test dd.equilibrium.vacuum_toroidal_field.r0 > 0

            # Summary
            @test dd.equilibrium.time[tslice] == dd.summary.time[tslice]
            @test dd.summary.global_quantities.r0.value > 0

            # 1d
            p1 = eqt.profiles_1d
            nprof = length(p1.psi)
            @test nprof > 10
            @test p1.psi[1] == gq.psi_axis
            @test p1.psi[end] == gq.psi_boundary
            @test length(p1.f) == nprof
            @test length(p1.pressure) == nprof
            @test length(p1.rho_tor_norm) == nprof

            # boundary
            r_bry = eqt.boundary.outline.r
            z_bry = eqt.boundary.outline.z
            @test length(r_bry) == length(z_bry)

            # 2d
            # Did the R and Z (dim1 and dim2) coordinates get written properly?
            p2 = eqt.profiles_2d[1]
            r_eq = p2.grid.dim1
            z_eq = p2.grid.dim2
            @test length(r_eq) > 10
            @test minimum(r_eq) > 0  # R should be dim1
            @test minimum(z_eq) < 0
            @test length(z_eq) > 10
            # Does psi match R and Z?
            println(size(p2.psi), (length(r_eq), length(z_eq)))
            @test size(p2.psi) == (length(r_eq), length(z_eq))
            # Does the psi grid look like it's transposed the right way?
            # Many equilibrium reconstructions have equal numbers of cells in R and Z,
            # so transposing the psi map incorrectly would not be detected by checking
            # its array dimensions. It's also possible to incorrectly associate the psi
            # map with R and Z (dim1 and dim2). So what recourse does that leave us?
            # Well, the points on the boundary outline should be at psi_N = 1.0 to
            # within some small tolerance.
            psin2d = (p2.psi .- gq.psi_axis) ./ (gq.psi_boundary - gq.psi_axis)
            tolerance = 2.0e-3  # It's not always a high res contour so cut some slack
            psin_bry =
                Interpolations.linear_interpolation((r_eq, z_eq), psin2d).(r_bry, z_bry)
            @test maximum(psin_bry) < (1.0 + tolerance)
            @test minimum(psin_bry) > (1.0 - tolerance)

            # derived
            @test gq.q_axis == p1.q[1]

            # X-points
            bx = eqt.boundary.x_point
            bsx = eqt.boundary_separatrix.x_point
            bssx = eqt.boundary_secondary_separatrix.x_point
            nxpt = length(bx)
            nprim = length(bsx)
            nsec = length(bssx)
            @test nxpt >= 1
            @test nxpt >= (nprim + nsec)

            # wall
            limiter = dd.wall.description_2d[1].limiter
            @test length(limiter.unit[1].outline.r) > 10
            @test length(limiter.unit[1].outline.r) == length(limiter.unit[1].outline.z)
        end
    end
end

if isempty(ARGS) || "prep" in ARGS
    @testset "preparation" begin
        eqdsk_file = "geqdsk_iter_small_sample"
        sample_paths = [
            splitdir(pathof(SOLPS2ctrl))[1] * "/../sample",
        ]
        core_method = "simple"
        edge_method = "simple"
        filename = splitdir(pathof(SOLPS2ctrl))[1] * "/../sd_input_data"
        output_format = "json"
        dd = SOLPS2ctrl.preparation(
            eqdsk_file,
            sample_paths;
            core_method=core_method,
            filename=filename,
            output_format=output_format,
        )
        out_file = filename * "." * output_format
        p2 = dd.equilibrium.time_slice[1].profiles_2d[1]
        psirz = p2.psi
        r = p2.grid.dim1
        z = p2.grid.dim2
        @test size(psirz) == (length(r), length(z))
        println(out_file)
        @test isfile(out_file)
        print("imas2json timing: ")
        @time IMAS.imas2json(dd, filename * ".json", strict=true, freeze=false)
    end
end

if isempty(ARGS) || "sysid" in ARGS
    @testset "system_id" begin
        case = "$(@__DIR__)/../sample/D3D_Lore_Time_Dep"
        println("Reading interferometer data...")
        print("json2imas time: ")
        @time ids = IMAS.json2imas(case * "/imas_with_ifo.json")
        # The file above "imas_with_ifo.json" has been created using following code.
        # We jump directly to after the interferometer addition to avoid making
        # FusionSyntheticDiagnostics a dependency of this repo.
        #
        # data_dir = case * "/puff2.5e21_ss_noLat_dribble_308_2d_output"
        # b2fgmtry = case * "/baserun/b2fgmtry"
        # b2time = data_dir * "/b2time.nc"
        # b2mn = data_dir * "/b2mn.dat"
        # fort = (data_dir * "/fort.33", data_dir * "/fort.34", data_dir * "/fort.35")
        # eqdsk = case * "/baserun/gfile"
        # SOLPS2ctrl.geqdsk_to_imas!(eqdsk, ids)
        # ids = SOLPS2imas.solps2imas(b2fgmtry, b2time; b2mn=b2mn, fort=fort)
        # SOLPS2ctrl.fill_in_extrapolated_core_profile!(
        #     ids,
        #     "electrons.density";
        #     method="simple",
        #     cell_subset_idx=-5,
        # )
        # FusionSyntheticDiagnostics.add_interferometer!(
        #     case * "/interferometer.json",
        #     ids;
        #     n_e_gsi=-5,
        #     overwrite=true,
        # )

        ifo_ch = 3
        ne_factor = 1e-19
        ne_offset = ids.interferometer.channel[ifo_ch].n_e_line_average.data[1]
        boltzmann_k = 1.38064852e-23
        Temp = 300.0

        ifo_tt = ids.interferometer.channel[1].n_e_line_average.time
        # Make a uniform linear time series
        tt = collect(ifo_tt[1]:0.001:ifo_tt[end])

        # Input gas data
        input_gas_data = readdlm(case * "/input_gas_atomsps.txt"; comments=true)
        gas_factor = 1e-21
        gas_offset = input_gas_data[1, 2]
        # Intepolate input data along the uniform time
        input_gas =
            Interpolations.linear_interpolation(
                input_gas_data[:, 1] .- 5.0,
                input_gas_data[:, 2],
            ).(
                tt,
            )

        # Output interferometer data
        ne_ifo = ids.interferometer.channel[ifo_ch].n_e_line_average.data
        # Interpolate output data along the uniform time
        ne_uniform = Interpolations.linear_interpolation(ifo_tt, ne_ifo).(tt)

        function inp_cond_step(
            gas_so::Float64,
            previous_gas_so::Float64,
            corr::Float64;
            p::Union{Float64, Vector{Float64}}=-0.17,
        )::Tuple{Float64, Float64}
            if p isa Array
                param = p[1]
            else
                param = p
            end
            if gas_so + corr > previous_gas_so
                corr = corr + (gas_so + corr - previous_gas_so) * param
            end
            return gas_so + corr, corr
        end

        """
            inp_cond(inp_gas; p=0.242)

        Condition the input to have this non-linearity. Only when the gas is increasing, the
        increase amount is changed by a factor of 'p' of the actual change. This gives a
        response of gas input when it is rising vs falling and thus constitute a non-linearity.
        """
        function inp_cond(
            gas_so::Array{Float64};
            p::Union{Float64, Vector{Float64}}=-0.17,
        )::Array{Float64}
            ret_vec = deepcopy(gas_so)
            corr = 0.0
            for ii ∈ 2:length(ret_vec)
                ret_vec[ii], corr = inp_cond_step(ret_vec[ii], ret_vec[ii-1], corr; p)
            end
            return ret_vec
        end

        # verbose only during CI testinf
        verbose = get(ENV, "GITHUB_ACTIONS", "") == "true"

        # Third order linear fit of the model
        linear_plant_3 = SOLPS2ctrl.system_id(
            input_gas, ne_uniform, tt, 3;
            inp_offset=gas_offset, inp_factor=gas_factor,
            out_offset=ne_offset, out_factor=ne_factor,
            verbose,
        )

        # Non-linear input + 3rd order linear fit
        nonlin_plant_3, p_opt = system_id_optimal_inp_cond(
            input_gas, ne_uniform, tt, 3, inp_cond, Dict{Symbol, Any}(:p => -0.2);
            inp_offset=gas_offset, inp_factor=gas_factor,
            out_offset=ne_offset, out_factor=ne_factor,
            verbose,
        )

        lin_out = model_evolve(
            linear_plant_3, input_gas;
            inp_offset=gas_offset, inp_factor=gas_factor,
            out_offset=ne_offset, out_factor=ne_factor,
        )

        nonlin_out = model_evolve(
            nonlin_plant_3, input_gas;
            inp_offset=gas_offset, inp_factor=gas_factor,
            out_offset=ne_offset, out_factor=ne_factor,
            inp_cond=inp_cond, inp_cond_kwargs=p_opt)

        @test sqrt(Statistics.mean((lin_out - ne_uniform) .^ 2)) < 1.2e18

        @test sqrt(Statistics.mean((nonlin_out - ne_uniform) .^ 2)) < 0.4e18
    end
end

if isempty(ARGS) || "state" in ARGS
    @testset "state_pred" begin
        nu = rand(1:10, 1)[1]
        ny = rand(1:nu, 1)[1]
        nstates = rand(1:10, 1)[1]
        hh = 20
        rnd_model = ssrand(ny, nu, nstates;
            proper=false,
            stable=true,
            Ts=0.001,
        )
        LL, MM, NN, OO = SOLPS2ctrl._calculate_LMNO(rnd_model, hh)

        x0 = rand(Float64, size(rnd_model.A, 1))
        test_u = randn((size(rnd_model.B, 2), hh))
        test_y, _, test_x, _ = lsim(rnd_model, test_u; x0)
        final_x = rnd_model.A * test_x[:, end] + rnd_model.B * test_u[:, end]

        test_U = zeros(size(rnd_model.B, 2) * hh)
        test_Y = zeros(size(rnd_model.C, 1) * hh)
        for i ∈ 1:hh
            test_U[(size(rnd_model.B, 2)*(i-1)+1):(size(rnd_model.B, 2)*i)] =
                test_u[:, i]
            test_Y[(size(rnd_model.C, 1)*(i-1)+1):(size(rnd_model.C, 1)*i)] =
                test_y[:, i]
        end
        test_est_Y = LL * x0 + MM * test_U
        test_est_x = NN * x0 + OO * test_U

        # println(sqrt(Statistics.mean((test_est_Y - test_Y) .^ 2)))
        # println(sqrt(Statistics.mean((test_est_x - final_x) .^ 2)))
        @test isapprox(test_est_Y, test_Y)
        @test isapprox(test_est_x, final_x)

        noise_scale = 1e-3

        set_u_off = randn() * noise_scale
        test_U_mod = test_U .+ set_u_off

        set_y_off = randn() * noise_scale
        test_Y_mod = test_Y .+ set_y_off

        YY2x, UU2x = state_prediction_matrices(rnd_model, hh)
        test_est_x2 = YY2x * test_Y_mod + UU2x * test_U_mod
        # println(sqrt(Statistics.mean((test_est_x2 - final_x) .^ 2)))
        @test isapprox(test_est_x2, final_x; rtol=10 * noise_scale)
    end
end

if isempty(ARGS) || "controller" in ARGS
    @testset "controller" begin
        # Random 2-pole plant model with resonance between 100 and 200 Hz
        # and Q value between 1 to 6 and a gain of 3.0
        s = tf('s')
        imp = rand(200.0:0.001:300.0)
        rep = rand(-200:0.001:-50)
        abs2 = rep^2 + imp^2
        den = s^2 - 2 * rep * s + abs2
        Ts = 0.001 # Sampling period of discrete system
        plant_model = c2d(ss(abs2 * 3 / den), Ts)
        delay_steps = 100
        delay_model = delay(delay_steps * Ts, Ts)

        # Define an actuator for the plant model
        # Actuator must be a mutable struct and daughter of Actuator type
        mutable struct TestAct <: Actuator
            latency::DelayedActuator
            gain::Float64
        end
        # And all actuators must be callable
        # Here we make this a simple gain actuator with some latency
        (ta::TestAct)(inp::Vector{Float64}) = ta.latency(inp .* ta.gain)
        # Now create an instance of this actuator with gain 5 and delay of 10 time steps
        act = TestAct(DelayedActuator(delay_steps), 5.0)

        # Finally, design a linear controller using pid
        lc = LinearController(c2d(ss(pid(0.1, 0.1 / 10)), Ts), zeros(Float64, 1))

        # Reduce controller gain to ensure stability taking into account the delay
        w = collect(logrange(0.1, 1e3, 400))
        pm = -200
        gm = -1
        while pm < 30 || gm < 3
            wgm, gm, wpm, pm =
                margin(plant_model * act.gain * delay_model * lc.ctrl_ss, w)
            pm = pm[1, 1, 1]
            gm = gm[1, 1, 1]
            if pm < 30 || gm < 3
                lc.ctrl_ss *= 0.9 # Reduce gain slightly
            end
        end

        # Create PVLC and another instance of actuator for this run
        pvlc = PVLC(c2d(ss(pid(1.0, 1.0 / 10)), Ts), zeros(Float64, 1), plant_model, 20)
        act2 = TestAct(DelayedActuator(delay_steps), 5.0)

        # Tune PVLC controller for no delay
        pm = -200
        gm = -1
        while pm < 30 || gm < 3
            w = collect(logrange(0.1, 1e3, 400))
            wgm, gm, wpm, pm = margin(plant_model * act2.gain * pvlc.ctrl_ss, w)
            pm = pm[1, 1, 1]
            gm = gm[1, 1, 1]
            if pm < 30 || gm < 3
                pvlc.ctrl_ss *= 0.9
            end
        end

        # Set a target waveform
        tt = collect(0:Ts:2)
        set_points_tt = [0.0, 0.1, 0.5, 0.7, 1.3, 1.5, 1.7, 2.0]
        set_points = [0.0, 0.0, 1.0, 1.0, -0.5, -0.5, 0.0, 0.0]
        target = Matrix{Float64}(
            Interpolations.linear_interpolation(set_points_tt, set_points).(tt)',
        )

        # Run closed loop simulation
        res = Dict()
        res["PI"] = run_closed_loop_sim(plant_model, act, lc, target)
        res["PVLC"] = run_closed_loop_sim(plant_model, act2, pvlc, target)

        p1 = plot(
            tt, res["PI"][:target][1, :];
            label="Target", linestyle=:dash, color=:black,
        )
        plot!(
            tt, res["PI"][:plant_out][1, :];
            label="PI", linewidth=2, color=:deepskyblue3,
        )
        plot!(
            tt, res["PVLC"][:plant_out][1, :];
            label="PVLC", ylabel="Plant Output", xformatter=_ -> "",
            linewidth=2, color=:orange,
        )
        p2 = plot(
            tt, res["PI"][:plant_inp][1, :];
            label="PI", linewidth=2, color=:deepskyblue3,
        )
        plot!(
            tt, res["PVLC"][:plant_inp][1, :];
            label="PVLC", ylabel="Plant Input", xformatter=_ -> "",
            linewidth=2, color=:orange,
        )
        p3 = plot(
            tt, res["PI"][:ctrl_out][1, :];
            label="PI", linewidth=2, color=:deepskyblue3,
        )
        plot!(
            tt, res["PVLC"][:ctrl_out][1, :];
            label="PVLC", ylabel="Controller\nOutput", xlabel="Time / s",
            linewidth=2, color=:orange,
        )
        l = @layout [a{0.5h}; b{0.25h}; c{0.25h}]
        plot(
            p1, p2, p3;
            layout=l,
            suptitle=(
                "Closed Loop Simulation Results\n" *
                "Actuator has delay of $(act.latency.delay * Ts)s."
            ),
        )

        savefig("$(@__DIR__)/Closed_Loop_Sim_Results.pdf")

        residual_PVLC = sqrt(
            Statistics.mean(
                (res["PVLC"][:target][1, :] .- res["PVLC"][:plant_out][1, :]) .^ 2,
            ),
        )

        residual_LC = sqrt(
            Statistics.mean(
                (res["PI"][:target][1, :] .- res["PI"][:plant_out][1, :]) .^ 2,
            ),
        )

        @test residual_PVLC < residual_LC
    end
end
