#!! Conversion from UPF !!#

function Psp8Header(file::UpfFile, rgrid=file.mesh.r)
    title = ""
    zatom = Float64(PeriodicTable.elements[Symbol(titlecase(file.header.element))].number)
    zion = file.header.z_valence
    pspd = -1
    pspcod = 8
    pspxc = libxc_to_abinit(libxc_string(file.header))
    lmax = file.header.l_max
    lloc = file.header.l_local == -1 ? file.header.l_max + 1 : file.header.l_local
    mmax = length(rgrid)  # UPF may need to be interpolated onto a linear mesh `rgrid`
    r2well = 0.0
    # Use the original UPF mesh to find rchrg
    rchrg = if file.header.core_correction
        for i in reverse(eachindex(file.nlcc))
            if file.nlcc[i] > 1e-6  # This is the tolerance used by ABINIT for UPF -> PSP8
                file.mesh.r[i]
            end
        end
        last(file.mesh.r)
    else
        0.0
    end
    fchrg = file.header.core_correction ? 1.0 : 0.0
    qchrg = 0.0
    nproj = map(0:(file.header.l_max)) do l
        return count(beta -> beta.angular_momentum == l, file.nonlocal.betas)
    end
    extension_switch = file.header.has_so ? 3 : 1  # 3 (spin-orbit and valence charge), 1 (valence charge)
    nprojso = nothing

    return Psp8Header(title, zatom, zion, pspd, pspcod, pspxc, lmax, lloc, mmax, r2well,
                      rchrg, fchrg, qchrg, nproj, extension_switch, nprojso)
end

function Psp8File(file::UpfFile, rgrid=file.mesh.r)
    if file.header.has_so
        # See: https://github.com/abinit/abinit/blob/master/src/57_iopsp_parser/m_pspheads.F90#L958
        error("Converting UPF with spin-orbit to PSP8 is currently not supported.")
    end
    if file.header.number_of_wfc > 0
        @warn "UPF file contains $(file.header.number_of_wfc) pseudo-waves. PSP8 does not support pseudo-waves, they will be lost."
    end
    if !all(isapprox.(diff(rgrid), rgrid[2] - rgrid[1]))
        error("UPF has non-uniform radial mesh. Provide a uniform mesh `rgrid` to interpolate the UPF file for conversion.")
    end

    identifier = file.identifier * ".psp8"
    header = Psp8Header(file, rgrid)

    v_local = extrapolate(
        interpolate(
            file.mesh.r, file.local_ ./ 2,  # Ry -> Ha
            BSplineOrder(4), Natural()
        ),
        Flat()
    ).(rgrid)

    projectors = [Vector{Float64}[] for _ in 0:(file.header.l_max)]
    ekb = [Float64[] for _ in 0:(file.header.l_max)]
    for l in 0:(file.header.l_max)
        for (i, beta) in enumerate(file.nonlocal.betas)
            if beta.angular_momentum == l
                beta_spline = extrapolate(
                    interpolate(
                        file.mesh.r, beta.beta ./ 2,  # Ry -> Ha; could also be done in D_{ij} => E_{KB}
                        BSplineOrder(4), Natural()
                    ),
                    Flat()
                )
                push!(projectors[l + 1], beta_spline.(rgrid))
                push!(ekb[l + 1], file.nonlocal.dij[i, i])
            end
        end
    end

    projectors_so = nothing  # Unsupported
    ekb_so = nothing         # Unsupported

    if file.header.core_correction
        rhoc_spline = extrapolate(
            interpolate(
                file.mesh.r, file.nlcc .* 4π,  # Add 4π prefactor
                BSplineOrder(6), Natural()
            ),
            Smooth()
        )
        rhoc = (Derivative(0) * rhoc_spline).(rgrid)
        d_rhoc_dr = (Derivative(1) * rhoc_spline).(rgrid)
        d2_rhoc_dr2 = (Derivative(2) * rhoc_spline).(rgrid)
        d3_rhoc_dr3 = (Derivative(3) * rhoc_spline).(rgrid)
        d4_rhoc_dr4 = (Derivative(4) * rhoc_spline).(rgrid)
    else
        rhoc = nothing
        d_rhoc_dr = nothing
        d2_rhoc_dr2 = nothing
        d3_rhoc_dr3 = nothing
        d4_rhoc_dr4 = nothing
    end

    i_begin = file.mesh.r[1] ≈ 0 ? 2 : 1
    rhov_spline = extrapolate(
        interpolate(
            file.mesh.r[i_begin:end], file.rhoatom[i_begin:end] ./ file.mesh.r[i_begin:end] .^ 2,  # Remove r^2 prefactor
            BSplineOrder(4), Natural()
        ),
        Flat()
    )
    rhov = rhov_spline.(rgrid)

    ae_rhov = nothing
    ae_rhoc = nothing

    return Psp8File(identifier, header, rgrid, v_local, projectors, ekb, projectors_so,
                    ekb_so, rhoc, d_rhoc_dr, d2_rhoc_dr2, d3_rhoc_dr3, d4_rhoc_dr4, rhov, ae_rhov, ae_rhoc)
end
