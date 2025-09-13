#!! Conversion from PSP8 !!#

function UpfHeader(header::Psp8Header)
    generated = ""
    author = ""
    date = string(header.pspd)
    comment = header.title
    element = String(PeriodicTable.elements[Int(header.zatom)].symbol)
    pseudo_type = "NC"
    relativistic = header.extension_switch in (2, 3) ? "full" : "scalar"
    is_ultrasoft = false
    is_paw = false
    is_coulomb = false
    has_so = header.extension_switch in (2, 3)
    has_wfc = false
    has_gipaw = false
    paw_as_gipaw = nothing
    core_correction = header.fchrg > 0
    functional = libxc_to_qe_old(libxc_string(header))
    z_valence = header.zion
    total_psenergy = nothing
    wfc_cutoff = nothing
    rho_cutoff = nothing
    l_max = header.lmax
    l_max_rho = nothing
    l_local = header.lloc > header.lmax ? -1 : header.lloc
    mesh_size = header.mmax
    number_of_wfc = 0
    number_of_proj = sum(header.nproj)
    return UpfHeader(generated, author, date, comment, element, pseudo_type, relativistic,
                     is_ultrasoft, is_paw, is_coulomb, has_so, has_wfc, has_gipaw,
                     paw_as_gipaw, core_correction, functional, z_valence, total_psenergy,
                     wfc_cutoff, rho_cutoff, l_max, l_max_rho, l_local, mesh_size,
                     number_of_wfc,
                     number_of_proj)
end
UpfHeader(file::Psp8File) = UpfHeader(file.header)

function UpfMesh(file::Psp8File)
    r = file.rgrid
    # PSP8 always uses a uniform mesh. This is not strictly allowed by the UPF specification,
    # but in practice, it works.
    rab = fill(r[2] - r[1], length(r))
    mesh = length(r)
    rmax = maximum(r)
    dx = nothing
    xmin = nothing
    zmesh = nothing
    return UpfMesh(r, rab, mesh, rmax, dx, xmin, zmesh)
end

function UpfNonlocal(file::Psp8File)
    n_projectors = sum(file.header.nproj)
    betas = UpfBeta[]
    dij = zeros(Float64, n_projectors, n_projectors)
    index = 1
    # Non-local projectors must be listed strictly in order of increasing angular momentum.
    for l_index in 1:(file.header.lmax)
        l = l_index - 1
        for (projector, ekb) in zip(file.projectors[l_index], file.ekb[l_index])
            push!(betas,
                  UpfBeta(projector .* 2,  # Ha -> Ry; could also be done in E_{KB} => D_{ij}
                          index,
                          l,
                          file.header.mmax,
                          last(file.rgrid),
                          nothing,
                          nothing,
                          nothing))
            dij[index, index] = ekb
            index += 1
        end
    end
    augmentation = nothing
    return UpfNonlocal(betas, dij, augmentation)
end

function UpfFile(file::Psp8File)
    if file.header.extension_switch in (2, 3)
        # See: https://github.com/abinit/abinit/blob/master/src/57_iopsp_parser/m_pspheads.F90#L958
        error("Converting PSP8 with spin-orbit to UPF is currently not supported.")
    end
    identifier = file.identifier * ".upf"
    version = "(FROM PSP8)"
    info = "Converted from PSP8 format by PseudoPotentialIO.jl $(PseudoPotentialIO.VERSION)"
    header = UpfHeader(file)
    mesh = UpfMesh(file)
    nlcc = file.rhoc ./ 4π  # Remove 4π prefactor
    local_ = file.v_local .* 2  # Ha -> Ry
    nonlocal = UpfNonlocal(file)
    pswfc = nothing
    full_wfc = nothing
    rhoatom = file.rgrid .^ 2 .* file.rhov  # Add r^2 prefactor
    spin_orb = nothing
    paw = nothing
    gipaw = nothing

    return UpfFile(identifier, version, info, header, mesh, nlcc, local_, nonlocal, pswfc, full_wfc,
                   rhoatom, spin_orb, paw, gipaw)
end
