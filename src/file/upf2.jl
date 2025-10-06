function upf2_parse_psp(io::IO; identifier="")
    text = read(io, String)
    # Remove end-of-file junk (input data, etc.)
    text = string(split(text, "</UPF>")[1], "</UPF>")
    # Clean any errant `&` characters
    text = replace(text, "&" => "")
    doc = parsexml(text)

    root_node = root(doc)
    version = get_attr(String, root_node, "version")
    #* PP_INFO
    info_node = findfirst("PP_INFO", root_node)
    info = isnothing(info_node) ? nothing : nodecontent(info_node)
    #* PP_HEADER
    header = upf2_parse_header(doc)
    #* PP_MESH
    mesh = upf2_parse_mesh(doc)
    #* PP_NLCC
    nlcc_node = findfirst("PP_NLCC", root_node)
    if isnothing(nlcc_node)
        nlcc = nothing
    else
        nlcc = parse_nodecontent(Float64, nlcc_node)
    end
    #* PP_LOCAL
    local_node = findfirst("PP_LOCAL", root_node)
    if isnothing(local_node) | header.is_coulomb
        local_ = nothing
    else
        local_ = parse_nodecontent(Float64, local_node)
    end
    #* PP_NONLOCAL
    nonlocal = upf2_parse_nonlocal(doc)
    #* PP_PSWFC
    pswfc_node = findfirst("PP_PSWFC", root_node)
    pswfc = [upf2_parse_chi(n)
             for n in eachnode(pswfc_node) if occursin("PP_CHI.", nodename(n))]
    #* PP_FULL_WFC
    if isnothing(findfirst("PP_FULL_WFC", root_node))
        full_wfc = nothing
    else
        full_wfc = upf2_parse_full_wfc(doc)
    end
    #* PP_RHOATOM
    rhoatom_node = findfirst("PP_RHOATOM", root_node)
    rhoatom = parse_nodecontent(Float64, rhoatom_node)
    #* PP_SPINORB
    if isnothing(findfirst("PP_SPIN_ORB", root_node))
        spinorb = nothing
    else
        spinorb = upf2_parse_spin_orb(doc)
    end
    #* PP_PAW
    if isnothing(findfirst("PP_PAW", root_node))
        paw = nothing
    else
        paw = upf2_parse_paw(doc)
    end
    #* PP_GIPAW
    if isnothing(findfirst("PP_GIPAW", root_node))
        gipaw = nothing
    else
        gipaw = upf2_parse_gipaw(doc)
    end

    return UpfFile(identifier, version, info, header, mesh, nlcc, local_,
                   nonlocal, pswfc, full_wfc, rhoatom, spinorb, paw, gipaw)
end

function upf2_parse_header(node::EzXML.Node)
    functional_dlm(x) = isspace(x) || x == '-'

    generated = get_attr(String, node, "generated")
    author = get_attr(String, node, "author")
    date = get_attr(String, node, "date")
    comment = get_attr(String, node, "comment")
    element = get_attr(String, node, "element")
    pseudo_type = get_attr(String, node, "pseudo_type")
    relativistic = get_attr(String, node, "relativistic")
    is_ultrasoft = get_attr(Bool, node, "is_ultrasoft")
    is_paw = get_attr(Bool, node, "is_paw")
    is_coulomb = get_attr(Bool, node, "is_coulomb")
    has_so = get_attr(Bool, node, "has_so")
    has_wfc = get_attr(Bool, node, "has_wfc")
    has_gipaw = get_attr(Bool, node, "has_gipaw")
    paw_as_gipaw = get_attr(Bool, node, "paw_as_gipaw")
    core_correction = get_attr(Bool, node, "core_correction")
    functional = join(filter(!isempty, split(get_attr(String, node, "functional"), functional_dlm)), ' ')
    z_valence = get_attr(Float64, node, "z_valence")
    total_psenergy = get_attr(Float64, node, "total_psenergy")
    wfc_cutoff = get_attr(Float64, node, "wfc_cutoff")
    rho_cutoff = get_attr(Float64, node, "rho_cutoff")
    l_max = get_attr(Int, node, "l_max")
    l_max_rho = get_attr(Int, node, "l_max_rho")
    l_local = get_attr(Int, node, "l_local")
    mesh_size = get_attr(Int, node, "mesh_size")
    number_of_wfc = get_attr(Int, node, "number_of_wfc")
    number_of_proj = get_attr(Int, node, "number_of_proj")

    pseudo_type == "SL" && error("Semilocal pseudopotentials are not supported")

    return UpfHeader(generated, author, date, comment, element, pseudo_type,
                     relativistic, is_ultrasoft, is_paw, is_coulomb, has_so, has_wfc,
                     has_gipaw, paw_as_gipaw, core_correction, functional, z_valence,
                     total_psenergy, wfc_cutoff, rho_cutoff, l_max, l_max_rho, l_local,
                     mesh_size, number_of_wfc, number_of_proj)
end

function upf2_parse_header(doc::EzXML.Document)
    return upf2_parse_header(findfirst("PP_HEADER", root(doc)))
end

function upf2_parse_mesh(node::EzXML.Node)
    # Metadata
    dx = get_attr(Float64, node, "dx")
    mesh = get_attr(Int, node, "mesh")
    xmin = get_attr(Float64, node, "xmin")
    rmax = get_attr(Float64, node, "rmax")
    zmesh = get_attr(Float64, node, "zmesh")
    # PP_R
    r_node = findfirst("PP_R", node)
    if isnothing(mesh)
        mesh = get_attr(Int, r_node, "size")
    end
    r = parse_nodecontent(Float64, r_node)  # Bohr
    # PP_RAB
    rab_node = findfirst("PP_RAB", node)
    rab = parse_nodecontent(Float64, rab_node)
    return UpfMesh(r, rab, mesh, rmax, dx, xmin, zmesh)
end
upf2_parse_mesh(doc::EzXML.Document) = upf2_parse_mesh(findfirst("PP_MESH", root(doc)))

function upf2_parse_qij(node::EzXML.Node)
    # Metadata
    first_index = get_attr(Int, node, "first_index")
    second_index = get_attr(Int, node, "second_index")
    composite_index = get_attr(Int, node, "composite_index")
    is_null = get_attr(Bool, node, "is_null")
    # PP_QIJ.$i.$j
    qij = parse_nodecontent(Float64, node)
    return UpfQij(qij, first_index, second_index, composite_index, is_null)
end

function upf2_parse_qijl(node::EzXML.Node)
    # Metadata
    angular_momentum = get_attr(Int, node, "angular_momentum")
    first_index = get_attr(Int, node, "first_index")
    second_index = get_attr(Int, node, "second_index")
    composite_index = get_attr(Int, node, "composite_index")
    is_null = get_attr(Bool, node, "is_null")
    # PP_QIJL.$i.$j
    qijl = parse_nodecontent(Float64, node)
    return UpfQijl(qijl, angular_momentum, first_index, second_index, composite_index,
                   is_null)
end

function upf2_parse_augmentation(node::EzXML.Node)
    q_with_l = get_attr(Bool, node, "q_with_l")
    nqf = get_attr(Int, node, "nqf")
    nqlc = get_attr(Float64, node, "nqlc")
    shape = get_attr(String, node, "shape")
    iraug = get_attr(Int, node, "iraug")
    raug = get_attr(Float64, node, "raug")
    l_max_aug = get_attr(Float64, node, "l_max_aug")
    augmentation_epsilon = get_attr(Float64, node, "augmentation_epsilon")
    cutoff_r = get_attr(Float64, node, "cutoff_r")
    cutoff_r_index = get_attr(Float64, node, "cutoff_r_index")

    q_node = findfirst("PP_Q", node)
    q_vector = parse_nodecontent(Float64, q_node)
    q_size = get_attr(Int, q_node, "size")
    nq = Int(sqrt(q_size))
    q = reshape(q_vector, nq, nq)

    multipoles_node = findfirst("PP_MULTIPOLES", node)
    if isnothing(multipoles_node)
        multipoles = nothing
    else
        multipoles = parse_nodecontent(Float64, multipoles_node)
    end

    qfcoef_node = findfirst("PP_QFCOEF", node)
    if isnothing(qfcoef_node)
        qfcoefs = nothing
    else
        error("Cannot parse UPF v2.0.1 with PP_QFCOEF")
        # qfcoefs = parse_nodecontent(Float64, qfcoef_node)
    end

    rinner_node = findfirst("PP_RINNER", node)
    if isnothing(rinner_node)
        rinner = nothing
    else
        error("Cannot parse UPF v2.0.1 with PP_RINNER")
        # rinner = parse_nodecontent(Float64, rinner_node)
    end

    qij_nodes = [n for n in eachnode(node) if occursin("PP_QIJ.", nodename(n))]
    if isempty(qij_nodes)
        qijs = nothing
    else
        error("Cannot parse UPF v2.0.1 with PP_QIJ.i.j")
        # qijs = upf2_parse_qij.(qij_nodes)
    end

    qijl_nodes = [n for n in eachnode(node) if occursin("PP_QIJL.", nodename(n))]
    if isempty(qijl_nodes)
        qijls = nothing
    else
        qijls = upf2_parse_qijl.(qijl_nodes)
    end

    @assert !isnothing(qijs) | !isnothing(qijls)

    return UpfAugmentation(q, multipoles, qfcoefs, rinner, qijs, qijls, q_with_l, nqf, nqlc,
                           shape, iraug, raug, l_max_aug, augmentation_epsilon, cutoff_r,
                           cutoff_r_index)
end
function upf2_parse_augmentation(doc::EzXML.Document)
    return upf2_parse_augmentation(findfirst("PP_NONLOCAL/PP_AUGMENTATION", root(doc)))
end

function upf2_parse_beta(node::EzXML.Node)
    # Metadata
    name = nodename(node)
    index = get_attr(String, node, "index")
    if (index == "*") | isnothing(index)
        # If two digits, will be written as "*" because of Fortran printing, and sometimes
        # it is missing entirely. In both cases, we parse the index from the node name
        # which has the form "BETA.(\d+)"
        index = parse(Int, split(name, ".")[2])
    else
        index = parse(Int, index)
    end
    angular_momentum = get_attr(Int, node, "angular_momentum")
    cutoff_radius_index = get_attr(Int, node, "cutoff_radius_index")
    cutoff_radius = get_attr(Float64, node, "cutoff_radius")
    norm_conserving_radius = get_attr(Float64, node, "norm_conserving_radius")
    ultrasoft_cutoff_radius = get_attr(Float64, node, "ultrasoft_cutoff_radius")
    label = get_attr(String, node, "label")
    # PP_BETA.$i
    #* Note: all the data are parsed, not just up until the cutoff radius index
    #* Note: it is the _user's_ responsibility to cut off the projector data
    beta = parse_nodecontent(Float64, node)  # [1:cutoff_radius_index]
    return UpfBeta(beta, index, angular_momentum, cutoff_radius_index, cutoff_radius,
                   norm_conserving_radius, ultrasoft_cutoff_radius, label)
end

function upf2_parse_nonlocal(node::EzXML.Node)
    beta_nodes = [n for n in eachnode(node) if occursin("PP_BETA.", nodename(n))]
    betas = upf2_parse_beta.(beta_nodes)

    dij_node = findfirst("PP_DIJ", node)
    dij = parse_nodecontent(Float64, dij_node)
    dij = reshape(dij, length(betas), length(betas))

    augmentation_node = findfirst("PP_AUGMENTATION", node)
    if isnothing(augmentation_node)
        augmentation = nothing
    else
        augmentation = upf2_parse_augmentation(augmentation_node)
    end

    return UpfNonlocal(betas, dij, augmentation)
end
function upf2_parse_nonlocal(doc::EzXML.Document)
    return upf2_parse_nonlocal(findfirst("PP_NONLOCAL", root(doc)))
end

function upf2_parse_chi(node::EzXML.Node)
    # Metadata
    l = get_attr(Int, node, "l")
    occupation = get_attr(Float64, node, "occupation")
    index = get_attr(Int, node, "index")
    label = get_attr(String, node, "label")
    n = get_attr(Int, node, "n")
    pseudo_energy = get_attr(Float64, node, "pseudo_energy")
    cutoff_radius = get_attr(Float64, node, "cutoff_radius")
    ultrasoft_cutoff_radius = get_attr(Float64, node, "ultrasoft_cutoff_radius")
    # PP_CHI.$i
    chi = parse_nodecontent(Float64, node)
    return UpfChi(chi, l, occupation, index, label, n, pseudo_energy, cutoff_radius,
                  ultrasoft_cutoff_radius)
end

function upf2_parse_relwfc(node::EzXML.Node)
    jchi = get_attr(Float64, node, "jchi")
    index = get_attr(Int, node, "index")
    els = get_attr(String, node, "els")
    nn = get_attr(Int, node, "nn")
    lchi = get_attr(Int, node, "lchi")
    oc = get_attr(Float64, node, "oc")
    return UpfRelWfc(jchi, index, els, nn, lchi, oc)
end

function upf2_parse_relbeta(node::EzXML.Node)
    index = get_attr(Int, node, "index")
    jjj = get_attr(Float64, node, "jjj")
    lll = get_attr(Int, node, "lll")
    return UpfRelBeta(index, jjj, lll)
end

function upf2_parse_spin_orb(node::EzXML.Node)
    relwfc_nodes = [n for n in eachnode(node) if occursin("PP_RELWFC.", nodename(n))]
    relwfcs = upf2_parse_relwfc.(relwfc_nodes)

    relbeta_nodes = [n for n in eachnode(node) if occursin("PP_RELBETA.", nodename(n))]
    relbetas = upf2_parse_relbeta.(relbeta_nodes)

    return UpfSpinOrb(relwfcs, relbetas)
end
function upf2_parse_spin_orb(doc::EzXML.Document)
    return upf2_parse_spin_orb(findfirst("PP_SPIN_ORB", root(doc)))
end

function upf2_parse_wfc(node::EzXML.Node)
    index = get_attr(Int, node, "index")
    # Sometimes the `index` attribute is missng, so we parse it from the node name which is
    # of the form "PSWFC.(\d+)"
    if isnothing(index)
        index = parse(Int, split(nodename(node), '.')[end])
    end
    l = get_attr(Int, node, "l")
    label = get_attr(String, node, "label")
    wfc = parse_nodecontent(Float64, node)
    return UpfWfc(wfc, index, l, label)
end

function upf2_parse_full_wfc(node::EzXML.Node)
    aewfc_nodes = [n for n in eachnode(node) if occursin("PP_AEWFC", nodename(n))]
    aewfcs = upf2_parse_wfc.(aewfc_nodes)

    pswfc_nodes = [n for n in eachnode(node) if occursin("PP_PSWFC", nodename(n))]
    pswfcs = upf2_parse_wfc.(pswfc_nodes)

    return UpfFullWfc(aewfcs, pswfcs)
end
function upf2_parse_full_wfc(doc::EzXML.Document)
    return upf2_parse_full_wfc(findfirst("PP_FULL_WFC", root(doc)))
end

function upf2_parse_paw(node::EzXML.Node)
    paw_data_format = get_attr(Int, node, "paw_data_format")
    core_energy = get_attr(Float64, node, "core_energy")

    occupations_node = findfirst("PP_OCCUPATIONS", node)
    occupations = parse_nodecontent(Float64, occupations_node)

    ae_nlcc_node = findfirst("PP_AE_NLCC", node)
    ae_nlcc = parse_nodecontent(Float64, ae_nlcc_node)

    ae_vloc_node = findfirst("PP_AE_VLOC", node)
    ae_vloc = parse_nodecontent(Float64, ae_vloc_node)

    aewfc_nodes = [n for n in eachnode(node) if occursin("PP_AEWFC", nodename(n))]
    aewfcs = upf2_parse_wfc.(aewfc_nodes)

    pswfc_nodes = [n for n in eachnode(node) if occursin("PP_PSWFC", nodename(n))]
    pswfcs = upf2_parse_wfc.(pswfc_nodes)

    return UpfPaw(paw_data_format, core_energy, occupations, ae_nlcc, ae_vloc, aewfcs,
                  pswfcs)
end
upf2_parse_paw(doc::EzXML.Document) = upf2_parse_paw(findfirst("PP_PAW", root(doc)))

function upf2_parse_gipaw_core_orbital(node::EzXML.Node)
    index = get_attr(Int, node, "index")
    label = get_attr(String, node, "label")
    # Sometimes these integers are printed as floats
    n = Int(get_attr(Float64, node, "n"))
    l = Int(get_attr(Float64, node, "l"))
    core_orbital = parse_nodecontent(Float64, node)
    return UpfGipawCoreOrbital(index, label, n, l, core_orbital)
end

function upf2_parse_gipaw(node::EzXML.Node)
    gipaw_data_format = get_attr(Int, node, "gipaw_data_format")
    core_orbitals_node = findfirst("PP_GIPAW_CORE_ORBITALS", node)
    core_orbital_nodes = [n
                          for n in eachnode(core_orbitals_node)
                          if occursin("PP_GIPAW_CORE_ORBITAL.", nodename(n))]
    core_orbitals = upf2_parse_gipaw_core_orbital.(core_orbital_nodes)
    return UpfGipaw(gipaw_data_format, core_orbitals)
end
upf2_parse_gipaw(doc::EzXML.Document) = upf2_parse_gipaw(findfirst("PP_GIPAW", root(doc)))

parse_bool(s::AbstractString)::Bool = occursin("T", uppercase(s)) ? true : false
parse_bool(s::Char)::Bool = uppercase(s) == 'T' ? true : false

function get_attr(::Type{T}, node::EzXML.Node, key;
                  default=nothing)::Union{Nothing,T} where {T<:Number}
    if haskey(node, key)
        value = strip(node[key])
        value = replace(uppercase(value), "D" => "E")
        attr = parse(T, value)
    else
        attr = default
    end
    return attr
end

function get_attr(::Type{T}, node::EzXML.Node, key;
                  default=nothing)::Union{Nothing,T} where {T<:AbstractString}
    return haskey(node, key) ? T(strip(node[key])) : default
end

function get_attr(::Type{Bool}, node::EzXML.Node, key; default=nothing)::Union{Nothing,Bool}
    return haskey(node, key) ? parse_bool(strip(node[key])) : default
end

function parse_nodecontent(::Type{T}, node::EzXML.Node) where {T}
    words = split(strip(nodecontent(node)))
    data = Vector{T}(undef, length(words))
    data .= parse.(T, words)
    return data
end
