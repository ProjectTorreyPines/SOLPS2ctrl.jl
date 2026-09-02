using Revise
using SOLPS2ctrl
using Test

if isempty(ARGS) || "PCS" in ARGS
    @testset "PCS instruction writing" begin
        n_states = 3
        n_inputs = 1
        n_outputs = 1
        n_history = 80
        type_code = 0  # density
        gas_species_atomic_number = 1
        secondary_gas_atomic_number = 0
        counter = 5
        prepared_by = "Testy McTestFace"
        model_description = "This is a test model that I made ALLLLL by myself. It is very nice."
        input_description = "input 1 is D2 flow in molecules / s"
        output_description = "output 1 is density in 10^19 m^-3"
        model_time_step_sec = 250.0e-6
        controller_gp = [5.39]
        controller_taui = [1.333]
        input_offset = [2.5e2]
        input_factor = [1.0e-3]
        output_offset = [6.17832]
        output_factor = [0.1]
        adapt_input_gp = [0.0]
        adapt_input_taui = [0.0]
        adapt_output_gp = [1.2]
        adapt_output_taui = [1.45]
        adapt_factor_min = 1e-2
        adapt_factor_max = 10.0
        a = [[1.2, 0.398, 0.002] [-0.0289, 0.89278, 0.28] [-0.0001, -0.54, 1.2]]
        b = [[0.99] [0.456] [-0.97]]
        c = reshape([1.23, 0.478, -0.002], 3, 1)
        d = zeros(1, 1)
        y2x = [
            [1.321, 1.2, 0.5]
            [1.22, 1.1, 0.6]
            [1.2, 1.2, 0.4]
            [1.1, 1.065, 0.45]
            [1.09, 1.054, 0.39]
            [1.293, 1.02, 0.38]
        ]
        u2x = [
            [0.978, 1.05, 0.43]
            [1.19, 1.02, 0.7]
            [1.4, 1.3, 0.2]
            [0.85, 1.01, 0.55]
            [0.91, 0.99, 0.25]
            [0.97, 0.74, 0.21]
        ]
        println(typeof(prepared_by))
        SOLPS2ctrl.write_pcs_config(
            a, b, c, d,
            y2x, u2x,
            n_inputs, n_outputs, n_states, n_history,
            type_code, gas_species_atomic_number,
            model_time_step_sec,
            adapt_input_gp,
            adapt_input_taui,
            adapt_output_gp,
            adapt_output_taui,
            adapt_factor_min,
            adapt_factor_max,
            input_offset,
            input_factor,
            input_description,
            output_offset,
            output_factor,
            output_description,
            controller_gp,
            controller_taui,
            model_description,
            counter,
            prepared_by,
            secondary_gas_atomic_number,
        )
    end
end
