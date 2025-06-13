"""
Tools for repairing/filling out partial equilibrium files.

Some of the added fields may not be totally accurate, so it is recommended to
use this tool mainly for test cases, as a utility. For a real equilibrium,
problems should be fixed properly.
"""

using Contour: Contour
using Statistics: Statistics

export add_rho_to_equilibrium!
export check_rho_1d

"""
    check_rho_1d(
        dd::IMAS.dd;
        time_slice::Int=1,
        throw_on_fail::Bool=false,
    )::Bool

Checks to see if rho exists and is valid in the equilibrium 1d profiles
"""
function check_rho_1d(
    dd::IMAS.dd;
    time_slice::Int=1,
    throw_on_fail::Bool=false,
)::Bool
    rho = dd.equilibrium.time_slice[time_slice].profiles_1d.rho_tor_norm
    if length(rho) < 1
        rho_okay = false
    elseif maximum(rho) == 0.0
        rho_okay = false
    else
        rho_okay = true
    end
    if (throw_on_fail) * (!rho_okay)
        if length(rho) < 1
            reason = "rho is missing"
        else
            reason = "rho is all zeros"
        end
        throw(
            ArgumentError(
                string(
                    "Equilibrium rho profile data at time index ", time_slice,
                    " is invalid. This is usually because rho is missing or all",
                    " zeros in the source file. In this case, ", reason, ". You",
                    " can try using add_rho_to_equilibrium!() to calculate rho ",
                    "and add it.",
                ),
            ),
        )
    end
    return rho_okay
end

"""
    function add_rho_to_equilibrium(dd:IMAS.dd)

Adds equilibrium rho profile to the DD
"""
function add_rho_to_equilibrium!(dd::IMAS.dd)
    nt = length(dd.equilibrium.time_slice)
    if nt < 1
        println("No equilibrium time slices to work with; can't add rho")
        return
    end
    for it ∈ 1:nt
        eqt = dd.equilibrium.time_slice[it]
        b0 = dd.equilibrium.vacuum_toroidal_field.b0[it]
        r0 = dd.equilibrium.vacuum_toroidal_field.r0[it]

        psi = eqt.profiles_1d.psi
        n = length(psi)
        if (length(eqt.profiles_1d.rho_tor_norm) > 0)
            if maximum(eqt.profiles_1d.rho_tor_norm) > 0
                println("Slice #", it, " already has a rho_tor_norm profile; skipping")
                continue
            end
        end
        if (
            if IMAS.ismissing(eqt.profiles_1d, :phi)
                true
            else
                IMAS.isempty(eqt.profiles_1d.phi)
            end
        )
            eqt.profiles_1d.phi = Array{Float64}(undef, n)
            psi2 = eqt.profiles_2d[1].psi
            req = collect(eqt.profiles_2d[1].grid.dim1)
            zeq = collect(eqt.profiles_2d[1].grid.dim2)
            if (length(req), length(zeq)) == size(psi2')
                psi2 = Matrix(psi2')
                println(
                    "transposed psi to make it compatible with r,z prior to contouring",
                )
            end
            if (length(req), length(zeq)) != size(psi2)
                println("Invalid equilibrium data. rho cannot be added.")
                println("2D psi profile does not match 2D grid dimensions:")
                println("    dim1 (R): ", length(req), ", ", size(req))
                println("    dim2 (Z): ", length(zeq), ", ", size(zeq))
                println("    psi2d   : ", size(psi2))
                return
            else
                println(
                    "Eq looks okay ",
                    (length(req), length(zeq)),
                    ", ",
                    size(psi2),
                    ". ",
                    (size(req), size(zeq)),
                )
            end
            for j ∈ 1:n
                contour_level = psi[j]
                if j == n
                    # The last contour has a X point if everything is lined up right,
                    # and that could get weird. We want a contour that only takes the
                    # core boundary part of the separatrix.
                    r = eqt.boundary.outline.r
                    z = eqt.boundary.outline.z
                else
                    c = Contour.contour(req, zeq, psi2, contour_level)
                    clines = Contour.lines(c)
                    line_avg_height = [
                        Statistics.mean([
                            clines[i].vertices[v][2] for
                            v ∈ 1:length(clines[i].vertices)
                        ]) for i ∈ 1:length(clines)
                    ]
                    cl = clines[argmin(abs.(line_avg_height))]
                    # Now just do the 2D integral of B within the area of the contour
                    r = [cl.vertices[v][1] for v ∈ 1:length(cl.vertices)]
                    z = [cl.vertices[v][2] for v ∈ 1:length(cl.vertices)]
                end
                # Get the upper & lower halves of the path, each going from rmin to rmax
                rmin = minimum(r)
                rmax = maximum(r)
                irmin = argmin(r)
                irmax = argmax(r)
                if irmax > irmin
                    r1 = r[irmin:irmax]
                    z1 = z[irmin:irmax]
                    r2 = [r[irmax:end]; r[1:irmin]]
                    z2 = [z[irmax:end]; z[1:irmin]]
                    ii1 = sortperm(r1)
                    ii2 = sortperm(r2)
                    r1 = r1[ii1]
                    z1 = z1[ii1]
                    r2 = r2[ii2]
                    z2 = z2[ii2]
                else
                    r1 = r[irmax:irmin]
                    z1 = z[irmax:irmin]
                    r2 = [r[irmin:end]; r[1:irmax]]
                    z2 = [z[irmin:end]; z[1:irmax]]
                    ii1 = sortperm(r1)
                    ii2 = sortperm(r2)
                    r1 = r1[ii1]
                    z1 = z1[ii1]
                    r2 = r2[ii2]
                    z2 = z2[ii2]
                end
                # Vacuum field simplification: B = R0 B0 / R
                Interpolations.deduplicate_knots!(r1)
                Interpolations.deduplicate_knots!(r2)
                z1i = Interpolations.linear_interpolation(r1, z1)
                z2i = Interpolations.linear_interpolation(r2, z2)
                rr = LinRange(rmin, rmax, 101)
                rc = (rr[1:(end-1)] + rr[2:end]) / 2.0
                integral_part_ = [
                    log(rr[i+1] / rr[i]) * abs(z1i(rc[i]) - z2i(rc[i])) for
                    i ∈ 1:length(rc)
                ]
                phi_ = r0 * b0 .* integral_part_
                eqt.profiles_1d.phi[j] = sum(phi_)
            end
        end
        eqt.profiles_1d.rho_tor = sqrt.(eqt.profiles_1d.phi / (π * b0))
        eqt.profiles_1d.rho_tor_norm =
            eqt.profiles_1d.rho_tor / eqt.profiles_1d.rho_tor[end]
    end
    return println("Rho has been added to the equilibrium")
end
