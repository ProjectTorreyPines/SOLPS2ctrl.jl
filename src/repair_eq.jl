"""
Tools for repairing/filling out partial equilibrium files.

Some of the added fields may not be totally accurate, so it is recommended to
use this tool mainly for test cases, as a utility. For a real equilibrium,
problems should be fixed properly.
"""

import Contour
import Statistics

function add_rho_to_equilibrium!(dd::OMAS.dd)
    nt = length(dd.equilibrium.time_slice)
    if nt < 1
        println("No equilibrium time slices to work with; can't add rho")
        return
    end
    for it = 1:nt
        eqt = dd.equilibrium.time_slice[it]
        b0 = dd.equilibrium.vacuum_toroidal_field.b0[it]
        r0 = dd.equilibrium.vacuum_toroidal_field.r0[it]

        psi = eqt.profiles_1d.psi
        n = length(psi)
        if (length(eqt.profiles_1d.rho_tor_norm) > 0)
            if max(eqt.profiles_1d.rho_tor_norm) > 0
                println("Slice #", it, " already has a rho_tor_norm profile; skipping")
                continue
            end
        end
        if length(eqt.profiles_1d.phi) == 0
            resize!(eqt.profiles_1d.phi, n)
            psi2 = eqt.profiles_2d[1].psi
            req = collect(eqt.profiles_2d[1].grid.dim1)
            zeq = collect(eqt.profiles_2d[1].grid.dim2)
            if (length(req), length(zeq)) == size(psi2')
                psi2 = Matrix(psi2')
                println("transposed psi to make it compatible with r,z prior to contouring")
            end
            if (length(req), length(zeq)) != size(psi2)
                println("Invalid equilibrium data. rho cannot be added.")
                println("2D psi profile does not match 2D grid dimensions:")
                println("    dim1 (R): ", length(req), ", ", size(req))
                println("    dim2 (Z): ", length(zeq), ", ", size(zeq))
                println("    psi2d   : ", size(psi2))
                return
            else
                println("Eq looks okay ", (length(req), length(zeq)), ", ", size(psi2), ". ", (size(req), size(zeq)))
            end
            for j = 1:n
                c = Contour.contour(req, zeq, psi2, psi[j])
                clines = Contour.lines(c)
                line_avg_height = [Statistics.mean([clines[i].vertices[v][2] for v in 1:length(clines[i].vertices)]) for i in 1:length(clines)]
                cl = clines[argmin(abs.(line_avg_height))]
                # Now just do the 2D integral of B within the area of the contour
                r = [cl.vertices[v][1] for v in 1:length(cl.vertices)]
                z = [cl.vertices[v][2] for v in 1:length(cl.vertices)]
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
                z1i = Interpolations.LinearInterpolation(r1, z1)
                z2i = Interpolations.LinearInterpolation(r2, z2)
                rr = LinRange(rmin, rmax, 101)
                rc = (rr[1:end-1] + rr[2:end]) / 2.0
                integral_part_ = [log(rr[i+1]/rr[i]) * abs(z1i(rc[i]) - z2i(rc[i])) for i in 1:length(rc)]
                phi_ = r0 * b0 .* integral_part_
                eqt.profiles_1d.phi[j] = sum(phi_)
            end
        end
        eqt.profiles_1d.rho_tor = sqrt.(eqt.profiles_1d.phi / (π * b0))
        eqt.profiles_1d.rho_tor_norm = eqt.profiles_1d.rho_tor / eqt.profiles_1d.rho_tor[end]
    end
    println("Rho has been added to the equilibrium")
end
