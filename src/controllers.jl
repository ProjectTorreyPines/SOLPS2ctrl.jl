import ControlSystemIdentification: PredictionStateSpace
import ControlSystemsBase: StateSpace
import LinearAlgebra: pinv
# import ControlSystemsBase: lsim, ss, pid, d2c_exact, d2c, StateSpace

export lsim_step, state_prediction_matrices

function lsim_step(
    sys::Union{PredictionStateSpace, StateSpace}, u::Vector{Float64};
    x0::Vector{Float64}=zeros(size(sys.A, 1)),
)::Tuple{Vector{Float64}, Vector{Float64}}
    x = sys.A * x0 + sys.B * u
    y = sys.C * x0 + sys.D * u
    return y, x
end

function _calculate_LMNO(
    sys::Union{PredictionStateSpace, StateSpace},
    h::Int=size(sys.A, 1),
)
    A = sys.A
    B = sys.B
    C = sys.C
    D = sys.D
    order = size(A, 1)
    ninps = size(B, 2)
    nouts = size(C, 1)

    L = zeros((h * nouts, order))
    M = zeros((h * nouts, h * ninps))
    O = zeros((order, h * ninps))
    Atoim1 = Matrix{Float64}(I, (order, order))

    for j ∈ 1:h
        # println("Filling M with D: ", (nouts*(j-1)+1):(nouts*j), ", ", (nouts*(j-1)+1):(nouts*j))
        M[(nouts*(j-1)+1):(nouts*j), (ninps*(j-1)+1):(ninps*j)] = D
    end
    for i ∈ 1:h
        # println("i = ", i)
        if i > 1
            Atoim1 = Atoim1 * A
        end
        CAtoim1 = C * Atoim1
        # println("Filling L[$((nouts*(i-1)+1):(nouts*i)), :] = CA^$(i-1)")
        L[(nouts*(i-1)+1):(nouts*i), :] = CAtoim1

        for j ∈ (i+1):h
            # println("Filling M[$((nouts*(j-1)+1):(nouts*j)), $((ninps*(j-1-i)+1):(ninps*(j-i)))] = CA^$(i-1)B")
            M[(nouts*(j-1)+1):(nouts*j), (ninps*(j-1-i)+1):(ninps*(j-i))] = CAtoim1 * B
        end
        # println("Filling O[:, $((ninps*(h-i)+1):(ninps*(h-i+1)))] = A^$(i-1)B")
        O[:, (ninps*(h-i)+1):(ninps*(h-i+1))] = Atoim1 * B
    end
    N = Atoim1 * A
    return L, M, N, O
end

"""
    state_prediction_matrices(
        sys::Union{PredictionStateSpace, StateSpace},
        h::Int=size(sys.A, 1),
    )

Calculate state prediction matrices for a linear system `sys` given `h` steps of input
and output history. This function returns two matrices, `Y2x` and `U2x` that can be
used to calculate least square fitted estimate of the current state vector of the
system.
`Y2x` has size: (State size of `sys`) x (No. of outputs times h)
`U2x` has size: (State size of `sys`) x (No. of inputs times h)
The estimated state vector is obtained by:
Y2x * Y + U2x * U
where Y is all the `h` outputs of system stacked into a single vector.
and U is all the `h` inputs to the system stacked into a single vector.
"""
function state_prediction_matrices(
    sys::Union{PredictionStateSpace, StateSpace},
    h::Int=size(sys.A, 1),
)
    L, M, N, O = _calculate_LMNO(sys, h)
    Linv = pinv(L)
    Y2x = N * Linv
    U2x = (O - Y2x * M)
    return Y2x, U2x
end
