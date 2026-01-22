function upf1_parse_psp(io::IO; identifier="")
    version = "1.old"
    info = upf1_parse_info(io)
    header = upf1_parse_header(io)
    mesh = upf1_parse_mesh(io, header.mesh_size)
    if header.core_correction
        nlcc = upf1_parse_nlcc(io, header.mesh_size)
    else
        nlcc = nothing
    end
    if !header.is_coulomb
        local_ = upf1_parse_local(io, header.mesh_size)
    else
        local_ = nothing
    end
    nonlocal = upf1_parse_nonlocal(io, header.pseudo_type, header.mesh_size,
                                   header.number_of_proj, header.l_max)
    pswfc = upf1_parse_pswfc(io, header.mesh_size, header.number_of_wfc)
    full_wfc = nothing  # Not supported by UPF v1
    rhoatom = upf1_parse_rhoatom(io, header.mesh_size)
    if has_tag(io, "<PP_ADDINFO>")
        spinorb = upf1_parse_addinfo(io, header.number_of_proj, header.number_of_wfc)
    else
        spinorb = nothing
    end
    paw = nothing  # Not supported by UPF v1
    gipaw = nothing  # Not supported by UPF v1
    return UpfFile(identifier, version, info, header, mesh, nlcc, local_,
                   nonlocal, pswfc, full_wfc, rhoatom, spinorb, paw, gipaw)
end

"""
Read the human-readable `PP_INFO` (comment) block describing the pseudopotential.
"""
function upf1_parse_info(io::IO)
    pos = position(io)
    info = String[]
    read_until(io, "<PP_INFO>")
    # Read the first line of PP_INFO then read until the end tag
    # TODO: this will break if there are no lines in PP_INFO
    line = string(strip(readline(io)))
    push!(info, line)
    while !occursin("</PP_INFO>", line)
        line = string(strip(readline(io)))
        push!(info, line)
    end
    seek(io, pos)
    return join(info, '\n')
end

function upf1_parse_header(io::IO)
    read_until(io, "<PP_HEADER>")

    # Line 1 contains the "Version Number" of the pseudopotential (not the file format)
    s = split(readline(io))
    _ = parse(Int, s[1])

    # Line 2 contains the IUPAC "Element" symbol (e.g. Br for Bromine)
    s = split(readline(io))
    element = strip(s[1])

    # Line 3 contains a code for the pseudopotential formalism: `US` or `USPP` for
    # ultrasoft, `PAW` for projector-augmented wave, `NC` for norm-convserving, `1/r` for
    # coulombic, or `SL` for semi-local
    s = split(readline(io))
    pseudo_type = s[1]
    is_ultrasoft = pseudo_type in ("PAW", "US", "USPP")
    is_paw = pseudo_type == "PAW"
    is_coulomb = pseudo_type == "1/r"
    pseudo_type == "SL" && error("Semilocal pseudopotentials are not supported")

    # Line 4 contains a boolean flag in the form of `T` or `F` for the presence of
    # non-linear core correction
    s = split(readline(io))
    core_correction = occursin("T", uppercase(s[1])) ? true : false

    # Line 5 contains the QuantumESPRESSO-style string detailing the exchange-correlation
    # functional.
    # The line looks something like:
    #     SLA  PW   PBX  PBC    PBE  Exchange-Correlation functional
    # We look at the first four words:
    #     SLA PW PBX PBC
    # which correspond to the LDA exchange, LDA correlation, GGA exchange, and GGA
    # correlation functionals.
    # We check that each word is a valid UPF functional string then add it to the cleaned
    # functional string.
    s = strip.(split(strip(readline(io))))
    functional = []
    functional_dictionaries = (
        # The first word could be LDA exchange or a short name
        [UPF_LDA_EXCHANGE..., UPF_SHORT_NAMES...],
        UPF_LDA_CORRELATION,
        UPF_GGA_EXCHANGE,
        UPF_GGA_CORRELATION
    )
    for (word, dictionary) in zip(s, functional_dictionaries)
        if in(word, [entry["name"] for entry in dictionary])
            push!(functional, lowercase(word))
        end
    end
    functional = join(functional, ' ')

    functional_words = filter(s_i -> length(s_i) <= 6, s)
    functional_words = functional_words[1:min(4, length(functional_words))]
    functional = join(functional_words, ' ')

    # Line 6 contains the pseudo-ionic valence charge
    s = split(readline(io))
    z_valence = parse(Float64, s[1])

    # Line 7 contains the total energy of the pseudo-ion (in Rydberg)
    s = split(readline(io))
    total_psenergy = parse(Float64, s[1])

    # Line 8 contains the suggested kinetic-energy cutoffs for the wavefunctions and
    # charge density (in Rydberg)
    s = split(readline(io))
    wfc_cutoff = parse(Float64, s[1])
    rho_cutoff = parse(Float64, s[2])

    # Line 9 contains the maximum angular momentum channel
    s = split(readline(io))
    l_max = parse(Int, s[1])

    # Line 10 contains the number of points on the radial mesh
    s = split(readline(io))
    mesh_size = parse(Int, s[1])

    # Line 11 contains the number of chi-functions and the number of KB
    # non-local projectors
    s = split(readline(io))
    number_of_wfc = parse(Int, s[1])
    number_of_proj = parse(Int, s[2])

    # These values are given by UPF v2 but not by v1. They are set to `nothing` or
    # a meaningful value depending on what UPF v1 supports
    generated = nothing
    author = nothing
    date = nothing
    comment = nothing
    paw_as_gipaw = nothing
    l_max_rho = nothing
    l_local = nothing
    has_so = has_tag(io, "<PP_ADDINFO>")
    has_gipaw = false
    relativistic = has_so ? "full" : "scalar"
    has_wfc = false

    return UpfHeader(generated, author, date, comment, element, pseudo_type, relativistic,
                     is_ultrasoft, is_paw, is_coulomb, has_so, has_wfc, has_gipaw,
                     paw_as_gipaw,
                     core_correction, functional, z_valence, total_psenergy, wfc_cutoff,
                     rho_cutoff, l_max, l_max_rho, l_local, mesh_size, number_of_wfc,
                     number_of_proj)
end

function upf1_parse_mesh(io::IO, mesh_size::Int)
    pos = position(io)
    read_until(io, "<PP_R>")
    r = read_mesh_data(Float64, io, mesh_size)

    read_until(io, "<PP_RAB>")
    rab = read_mesh_data(Float64, io, mesh_size)
    seek(io, pos)

    # These values are given by UPF v2 but not by v1. They could be found by fitting
    # a mesh function to `r` and `rab` but are currently left undetermined.
    mesh = nothing
    rmax = nothing
    dx = nothing
    xmin = nothing
    zmesh = nothing

    return UpfMesh(r, rab, mesh, rmax, dx, xmin, zmesh)
end

function upf1_parse_nlcc(io::IO, mesh_size::Int)
    pos = position(io)
    read_until(io, "<PP_NLCC>")
    nlcc = read_mesh_data(Float64, io, mesh_size)
    seek(io, pos)
    return nlcc
end

function upf1_parse_local(io::IO, mesh_size::Int)
    pos = position(io)
    read_until(io, "<PP_LOCAL>")
    local_ = read_mesh_data(Float64, io, mesh_size)
    seek(io, pos)
    return local_
end

function upf1_parse_betas(io::IO, number_of_proj::Int, mesh_size::Int)
    pos = position(io)
    betas = Vector{UpfBeta}(undef, number_of_proj)
    for index in 1:number_of_proj
        read_until(io, "<PP_BETA>")

        # The first line of each block contains the index of the block and the angular
        # momentum of the KB projector in the block
        s = split(readline(io))
        index = parse(Int, s[1])
        angular_momentum = parse(Int, s[2])

        # The second line contains the index of the cutoff radius on the radial grid.
        s = readline(io)
        cutoff_radius_index = parse(Int, s)

        beta = read_mesh_data(Float64, io, cutoff_radius_index)

        # make beta size consistent with mesh_size
        beta = [beta; zeros(Float64, mesh_size - cutoff_radius_index)]

        # These values are (sometimes) given by UPF v2 but not by v1. The cutoff radius
        # could be found by looking up the value in `r` corresponding to
        # `cutoff_radius_index`, but it is currently left undetermined.
        cutoff_radius = nothing
        norm_conserving_radius = nothing
        ultrasoft_cutoff_radius = nothing
        label = nothing

        beta = UpfBeta(beta, index, angular_momentum, cutoff_radius_index, cutoff_radius,
                       norm_conserving_radius, ultrasoft_cutoff_radius, label)
        betas[index] = beta
    end
    seek(io, pos)
    return betas
end

function upf1_parse_dij(io::IO, number_of_proj::Int)
    pos = position(io)
    read_until(io, "<PP_DIJ>")

    dij = zeros(Float64, number_of_proj, number_of_proj)

    # The first line contains the number of non-zero projector coupling constants, i.e.
    # the number of data lines in this block
    s = split(readline(io))
    nd = parse(Int, s[1])  # Number of non-zero Dij components

    for _ in 1:nd
        # Each line contains three values: the index of the first projector, the index
        # of the second projector, and the value of the coupling constant between them
        s = split(readline(io))
        i, j = parse.(Int, s[1:2])
        d = parse(Float64, s[3])
        dij[i, j] = d
        dij[j, i] = d
    end
    seek(io, pos)
    return dij
end

function upf1_parse_augmentation(io::IO, mesh_size::Int, number_of_proj::Int, l_max::Int)
    pos = position(io)
    read_until(io, "<PP_QIJ>")

    # The first line contains the number of augmentation functions (`Number of Q Functions`)
    s = split(readline(io))
    nqf = parse(Int, s[1])

    if nqf == 0
        return nothing
    end

    # Read the block containing the inner radial cutoffs for each m quantum number of
    # the maximum angular momentum channel. Inside these radii, the augmentation functions
    # are determined by a set of polynomial coefficients `qfcoef`
    read_until(io, "<PP_RINNER>")
    rinner = Vector{Float64}(undef, 2l_max + 1)
    for i in 1:(2l_max + 1)
        s = split(readline(io))
        rinner[i] = parse(Float64, s[2])
    end
    read_until(io, "</PP_RINNER>")

    q = zeros(Float64, number_of_proj, number_of_proj)
    qijs = UpfQij[]
    qfcoefs = UpfQfcoef[]

    # For each pair of projectors, read the augmentation function block
    for i in 1:number_of_proj, _ in i:number_of_proj
        # The first line contains the indices of the two projectors and the angular momentum
        # of the projector corresponding to the second index
        s = split(readline(io))
        first_index, second_index, second_index_angular_momentum = parse.(Int, s[1:3])

        # The second line contains the integral of the augmentation charge
        s = split(readline(io))
        Q_int = parse(Float64, s[1])  # integral of augmentation charge

        qij = read_mesh_data(Float64, io, mesh_size)

        # Read the `nqf` polynomial coefficients for each magnetic quantum number of the
        # maximum angular momentum channel
        read_until(io, "<PP_QFCOEF>")
        qfcoef = read_mesh_data(Float64, io, nqf * (2l_max + 1))
        read_until(io, "</PP_QFCOEF>")

        # These values are given by UPF v2 but not by UPF v1
        is_null = nothing
        composite_index = nothing

        q[first_index, second_index] = Q_int
        q[second_index, first_index] = Q_int
        push!(qijs, UpfQij(qij, first_index, second_index, composite_index, is_null))
        push!(qfcoefs, UpfQfcoef(qfcoef, first_index, second_index, composite_index))
    end

    # These values are (sometimes) given by UPF v2 but not by UPF v1
    multipoles = nothing
    qijls = nothing
    q_with_l = false
    nqlc = nothing
    shape = nothing
    iraug = nothing
    raug = nothing
    l_max_aug = nothing
    augmentation_epsilon = nothing
    cutoff_r = nothing
    cutoff_r_index = nothing

    augmentation = UpfAugmentation(q, multipoles, qfcoefs, rinner, qijs, qijls, q_with_l,
                                   nqf, nqlc, shape, iraug, raug, l_max_aug,
                                   augmentation_epsilon, cutoff_r, cutoff_r_index)
    seek(io, pos)
    return augmentation
end

function upf1_parse_nonlocal(io::IO, pseudo_type::String, mesh_size::Int,
                             number_of_proj::Int, l_max::Int)
    # The `PP_NONLOCAL` superblock contains the KB non-local projectors, projector coupling
    # coefficients, and (sometimes) augmentation charges
    betas = upf1_parse_betas(io, number_of_proj, mesh_size)
    dij = upf1_parse_dij(io, number_of_proj)
    if pseudo_type in ("US", "USPP", "PAW")  # Ultrasoft and PAW have augmentation charges
        augmentation = upf1_parse_augmentation(io, mesh_size, number_of_proj, l_max)
    else  # All other formalisms do not have augmentation charges
        augmentation = nothing
    end
    return UpfNonlocal(betas, dij, augmentation)
end

function upf1_parse_pswfc(io::IO, mesh_size::Int, number_of_wfc::Int)
    pos = position(io)
    read_until(io, "<PP_PSWFC>")
    pswfc = Vector{UpfChi}(undef, number_of_wfc)
    for index in 1:number_of_wfc
        # The first line of each block contains a label (e.g. `3S`), the angular momentum,
        # and the occupation (in the range 0.0 - 2.0) of the pseudo-atomic wavefunction
        s = split(readline(io))
        label = s[1]
        l = parse(Int, s[2])
        occupation = parse(Float64, s[3])

        chi = read_mesh_data(Float64, io, mesh_size)

        # These values are (sometimes) given by UPF v2 but not by UPF v2
        n = nothing
        pseudo_energy = nothing
        cutoff_radius = nothing
        ultrasoft_cutoff_radius = nothing

        pswfc[index] = UpfChi(chi, l, occupation, index, label, n, pseudo_energy,
                              cutoff_radius, ultrasoft_cutoff_radius)
    end
    seek(io, pos)
    return pswfc
end

function upf1_parse_rhoatom(io::IO, mesh_size::Int)
    pos = position(io)
    read_until(io, "<PP_RHOATOM>")
    rhoatom = read_mesh_data(Float64, io, mesh_size)
    seek(io, pos)
    return rhoatom
end

function upf1_parse_addinfo(io::IO, number_of_proj::Int, number_of_wfc::Int)
    # The PP_ADDINFO block contains optional spin-orbit coupling information. It has two
    # sub-blocks, in order: pseudo-atomic wavfunctions and KB non-local projectors
    pos = position(io)
    read_until(io, "<PP_ADDINFO>")
    relwfcs = Vector{UpfRelWfc}(undef, number_of_wfc)
    for index in 1:number_of_wfc
        s = split(readline(io))
        els = strip(s[1])  # Label
        nn = parse(Int, s[2])  # Principal quantum number
        lchi = parse(Int, s[3])  # Angular momentum quantum number
        jchi = parse(Float64, s[4])  # Total angular momentum quantum number
        oc = parse(Float64, s[5])  # Occupation
        relwfcs[index] = UpfRelWfc(jchi, index, els, nn, lchi, oc)
    end
    relbetas = Vector{UpfRelBeta}(undef, number_of_proj)
    for index in 1:number_of_proj
        s = split(readline(io))
        lll = parse(Int, s[1])  # Angular momentum quantum number
        jjj = parse(Float64, s[2])  # Total angular momentum quantum number
        relbetas[index] = UpfRelBeta(index, jjj, lll)
    end
    seek(io, pos)
    return UpfSpinOrb(relwfcs, relbetas)
end

"""
Read through the `io` line-by-line until a line contains the `tag` then return. The `io`
read head will then be at the line following the `tag`.

Example:
    >    <tag1>
            some content
        <\tag1>
        <tag2>
            some content
        <\tag2>
        EOF

    read_until(io, "tag2")

         <tag1>
             some content
         <\tag1>
         <tag2>
    >        some content
         <\tag2>
         EOF
"""
function read_until(io::IO, tag::AbstractString)
    while true
        line = readline(io; keep=true)
        if isempty(line)
            throw(EOFError())
        end
        if occursin(tag, line)
            return
        end
    end
end

"""
Check if the `tag` occurs in the file contents of `io`. On return, the `io` read haed will
be returned to the same position as at call time.
"""
function has_tag(io::IO, tag::AbstractString)::Bool
    pos = position(io)
    try
        read_until(io, tag)
        return true
    catch e
        !isa(e, EOFError) && rethrow(e)
        return false
    finally
        seek(io, pos)
    end
end

"""
Read and parse vector data on an `n`-length mesh as type `T` from `io` regardless of whether
the data are wrapped (i.e. separated by spaces/tabs _or_ newlines). Note that the data are
read line-by-line, so partially reading lines is _not_ supported.

Example:
    >    1.0 2.0 3.0
         4.0 5.0
         6.0 7.0 8.0
         EOF

    read_mesh_data(Int, io, 5)

    [1, 2, 3, 4, 5]

         1.0 2.0 3.0
         4.0 5.0
    >    6.0 7.0 8.0
         EOF
"""
function read_mesh_data(T::Type, io::IO, n::Integer)
    mesh_data = T[]
    while length(mesh_data) < n
        line_data = parse.(T, split(readline(io)))
        append!(mesh_data, line_data)
    end
    length(mesh_data) != n && error("Too many values.")
    return mesh_data
end
