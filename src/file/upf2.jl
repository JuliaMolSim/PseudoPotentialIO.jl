function upf2_parse_psp(io::IO; identifier="")
    text = read(io, String)
    # Remove end-of-file junk (input data, etc.)
    text = string(split(text, "</UPF>")[1], "</UPF>")
    text = escape_ampersand(text)

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

function upf2_dump_psp(upffile::UpfFile)::EzXML.Node
    root_node = ElementNode("UPF")
    set_attr!(root_node, "version", "2.0.1")

    #* PP_INFO
    if !isnothing(upffile.info)
        addelement!(root_node, "PP_INFO", upffile.info)
    end
    # PP_HEADER
    link!(root_node, upf2_dump_header(upffile.header))
    # PP_MESH
    link!(root_node, upf2_dump_mesh(upffile.mesh))
    # PP_NLCC
    if !isnothing(upffile.nlcc)
        addelement!(root_node, "PP_NLCC", array_to_text(upffile.nlcc))
    end
    # PP_LOCAL
    if !isnothing(upffile.local_)
        addelement!(root_node, "PP_LOCAL", array_to_text(upffile.local_))
    end
    # PP_NONLOCAL
    link!(root_node, upf2_dump_nonlocal(upffile.nonlocal))
    # PP_PSWFC
    if !isnothing(upffile.pswfc)
        pswfc_node = ElementNode("PP_PSWFC")
        for chi in upffile.pswfc
            link!(pswfc_node, upf2_dump_chi(chi))
        end
        link!(root_node, pswfc_node)
    end
    # PP_FULL_WFC
    if !isnothing(upffile.full_wfc)
        link!(root_node, upf2_dump_full_wfc(upffile.full_wfc))
    end
    # PP_RHOATOM
    addelement!(root_node, "PP_RHOATOM", array_to_text(upffile.rhoatom))
    # PP_SPINORB
    if !isnothing(upffile.spin_orb)
        link!(root_node, upf2_dump_spin_orb(upffile.spin_orb))
    end

    # PP_PAW
    if !isnothing(upffile.paw)
        link!(root_node, upf2_dump_paw(upffile.paw))
    end

    # PP_GIPAW
    if !isnothing(upffile.gipaw)
        link!(root_node, upf2_dump_gipaw(upffile.gipaw))
    end

    return root_node
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

function upf2_dump_header(h::UpfHeader)::EzXML.Node
    node = ElementNode("PP_HEADER")

    for n in fieldnames(UpfHeader)
        set_attr!(node, n, getfield(h, n))
    end

    return node
end

function upf2_parse_mesh(node::EzXML.Node)
    # Metadata
    dx = get_attr(Float64, node, "dx")
    mesh = get_attr(Int, node, "mesh")
    xmin = get_attr(Float64, node, "xmin")
    rmax = get_attr(Float64, node, "rmax")
    zmesh = get_attr(Float64, node, "zmesh")

    # Parse from isolated node
    if isnothing(node.document)
        doc = XMLDocument()
        setroot!(doc, node)
    end

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

function upf2_dump_mesh(m::UpfMesh)::EzXML.Node
    node = ElementNode("PP_MESH")

    for n in [n for n in fieldnames(UpfMesh) if n != :r && n != :rab]
        set_attr!(node, n, getfield(m, n))
    end

    # PP_R
    addelement!(node, "PP_R", array_to_text(m.r))

    # PP_RAB
    addelement!(node, "PP_RAB", array_to_text(m.rab))

    return node
end

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

function upf2_dump_qij(qij::UpfQij)::EzXML.Node
    node = ElementNode("PP_QIJ")
    set_attr!(node, "first_index", qij.first_index)
    set_attr!(node, "second_index", qij.second_index)
    set_attr!(node, "composite_index", qij.composite_index)
    set_attr!(node, "is_null", qij.is_null)
    text = array_to_text(qij.qij)
    link!(node, TextNode(text))

    return node
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

function upf2_dump_qijl(qijl::UpfQijl)::EzXML.Node
    node = ElementNode("PP_QIJL.$(qijl.first_index).$(qijl.second_index).$(qijl.angular_momentum)")

    for n in fieldnames(UpfQijl)
        set_attr!(node, n, getfield(qijl, n))
    end

    text = array_to_text(qijl.qijl)
    link!(node, TextNode(text))
    return node
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

function upf2_dump_augmentation(aug::UpfAugmentation)::EzXML.Node
    node = ElementNode("PP_AUGMENTATION")

    for n in [n for n in fieldnames(UpfAugmentation) if
        n != :qij && n != :qijl && n != :qfcoefs && n != :rinner && n != :multipoles && n != :q]
        set_attr!(node, n, getfield(aug, n))
    end

    # PP_Q
    q_node = ElementNode("PP_Q")
    set_attr!(q_node, "size", length(aug.q))
    text = array_to_text(aug.q)
    link!(q_node, TextNode(text))
    link!(node, q_node)

    # PP_MULTIPOLES
    if !isnothing(aug.multipoles)
        addelement!(node, "PP_MULTIPOLES", array_to_text(aug.multipoles))
    end

    # PP_QFCOEF
    if !isnothing(aug.qfcoefs)
        # addelement!(node, "PP_QFCOEF", array_to_text(aug.qfcoefs))
        @warn "Cannot dump UPF v2.0.1 with PP_QFCOEF"
    end

    # PP_RINNER
    if !isnothing(aug.rinner)
        # addelement!(node, "PP_RINNER", array_to_text(aug.rinner))
        @warn "Cannot dump UPF v2.0.1 with PP_RINNER"
    end

    # PP_QIJ
    if !isnothing(aug.qijs)
        #for qij in aug.qijs
        #    link!(node, upf2_dump_qij(qij))
        #end
        @warn "Cannot dump UPF v2.0.1 with PP_QIJ.i.j"
    end

    # PP_QIJL
    if !isnothing(aug.qijls)
        for qijl in aug.qijls
            link!(node, upf2_dump_qijl(qijl))
        end
    end

    return node
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

function upf2_dump_beta(beta::UpfBeta)::EzXML.Node
    node = ElementNode("PP_BETA.$(beta.index)")

    for n in [n for n in fieldnames(UpfBeta) if n != :beta]
        set_attr!(node, n, getfield(beta, n))
    end

    text = array_to_text(beta.beta)
    link!(node, TextNode(text))
    return node
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

function upf2_dump_nonlocal(nl::UpfNonlocal)::EzXML.Node
    node = ElementNode("PP_NONLOCAL")

    # PP_BETA
    for beta in nl.betas
        link!(node, upf2_dump_beta(beta))
    end

    # PP_DIJ
    addelement!(node, "PP_DIJ", array_to_text(nl.dij))

    # PP_AUGMENTATION
    if !isnothing(nl.augmentation)
        link!(node, upf2_dump_augmentation(nl.augmentation))
    end

    return node
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

function upf2_dump_chi(chi::UpfChi)::EzXML.Node
    node = ElementNode("PP_CHI.$(chi.index)")

    for n in [n for n in fieldnames(UpfChi) if n != :chi]
        set_attr!(node, n, getfield(chi, n))
    end

    text = array_to_text(chi.chi)
    link!(node, TextNode(text))
    return node
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

function upf2_dump_relwfc(relwfc::UpfRelWfc)::EzXML.Node
    node = ElementNode("PP_RELWFC")

    for n in fieldnames(UpfRelWfc)
        set_attr!(node, n, getfield(relwfc, n))
    end

    return node
end

function upf2_parse_relbeta(node::EzXML.Node)
    index = get_attr(Int, node, "index")
    jjj = get_attr(Float64, node, "jjj")
    lll = get_attr(Int, node, "lll")
    return UpfRelBeta(index, jjj, lll)
end

function upf2_dump_relbeta(relbeta::UpfRelBeta)::EzXML.Node
    node = ElementNode("PP_RELBETA")
    set_attr!(node, "index", relbeta.index)
    set_attr!(node, "jjj", relbeta.jjj)
    set_attr!(node, "lll", relbeta.lll)
    return node
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

function upf2_dump_spin_orb(so::UpfSpinOrb)::EzXML.Node
    node = ElementNode("PP_SPIN_ORB")
    for relwfc in so.relwfcs
        link!(node, upf2_dump_relwfc(relwfc))
    end
    for relbeta in so.relbetas
        link!(node, upf2_dump_relbeta(relbeta))
    end
    return node
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

function upf2_dump_wfc(wfc::UpfWfc)::EzXML.Node
    node = ElementNode("PP_PSWFC.$(wfc.index)")
    set_attr!(node, "index", wfc.index)
    set_attr!(node, "l", wfc.l)
    set_attr!(node, "label", wfc.label)
    text = array_to_text(wfc.wfc)
    link!(node, TextNode(text))
    return node
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

function upf2_dump_full_wfc(fwfc::UpfFullWfc)::EzXML.Node
    node = ElementNode("PP_FULL_WFC")
    for aewfc in fwfc.aewfcs
        link!(node, upf2_dump_wfc(aewfc))
    end
    for pswfc in fwfc.pswfcs
        link!(node, upf2_dump_wfc(pswfc))
    end
    return node
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

function upf2_dump_paw(paw::UpfPaw)::EzXML.Node
    node = ElementNode("PP_PAW")
    set_attr!(node, "paw_data_format", paw.paw_data_format)
    set_attr!(node, "core_energy", paw.core_energy)

    # PP_OCCUPATIONS
    addelement!(node, "PP_OCCUPATIONS", array_to_text(paw.occupations))

    # PP_AE_NLCC
    addelement!(node, "PP_AE_NLCC", array_to_text(paw.ae_nlcc))

    # PP_AE_VLOC
    addelement!(node, "PP_AE_VLOC", array_to_text(paw.ae_vloc))

    # PP_AEWFC
    for aewfc in paw.aewfcs
        link!(node, upf2_dump_wfc(aewfc))
    end

    # PP_PSWFC
    for pswfc in paw.pswfcs
        link!(node, upf2_dump_wfc(pswfc))
    end

    return node
end

function upf2_parse_gipaw_core_orbital(node::EzXML.Node)
    index = get_attr(Int, node, "index")
    label = get_attr(String, node, "label")
    # Sometimes these integers are printed as floats
    n = Int(get_attr(Float64, node, "n"))
    l = Int(get_attr(Float64, node, "l"))
    core_orbital = parse_nodecontent(Float64, node)
    return UpfGipawCoreOrbital(index, label, n, l, core_orbital)
end

function upf2_dump_gipaw_core_orbital(core_orbital::UpfGipawCoreOrbital)::EzXML.Node
    node = ElementNode("PP_GIPAW_CORE_ORBITAL.$(core_orbital.index)")

    for n in [n for n in fieldnames(UpfGipawCoreOrbital) if n != :core_orbital]
        set_attr!(node, n, getfield(core_orbital, n))
    end

    text = array_to_text(core_orbital.core_orbital)
    link!(node, TextNode(text))
    return node
end

function upf2_parse_gipaw_orbital(node::EzXML.Node)
    index = get_attr(Int, node, "index")
    label = get_attr(String, node, "label")
    l = Int(get_attr(Float64, node, "l"))
    cutoff_radius = get_attr(Float64, node, "cutoff_radius")
    ultrasoft_cutoff_radius = get_attr(Float64, node, "ultrasoft_cutoff_radius")
    wf_ae_node = findfirst("PP_GIPAW_WFS_AE", node)
    wf_ae = parse_nodecontent(Float64, wf_ae_node)
    wf_ps_node = findfirst("PP_GIPAW_WFS_PS", node)
    wf_ps = parse_nodecontent(Float64, wf_ps_node)
    return UpfGipawOrbital(index, label, l, cutoff_radius, ultrasoft_cutoff_radius, wf_ae, wf_ps)
end

function upf2_dump_gipaw_orbital(orbital::UpfGipawOrbital)::EzXML.Node
    node = ElementNode("PP_GIPAW_ORBITAL.$(orbital.index)")
    for n in [n for n in fieldnames(UpfGipawOrbital) if n != :wf_ae && n != :wf_ps]
        value = getfield(orbital, n)
        if !isnothing(value)
            set_attr!(node, n, value)
        end
    end
    node_ae = ElementNode("PP_GIPAW_WFS_AE")
    text_ae = array_to_text(orbital.wf_ae)
    link!(node_ae, TextNode(text_ae))
    link!(node, node_ae)
    node_ps = ElementNode("PP_GIPAW_WFS_PS")
    text_ps = array_to_text(orbital.wf_ps)
    link!(node_ps, TextNode(text_ps))
    link!(node, node_ps)
    return node
end

function upf2_parse_gipaw_vlocal(node::EzXML.Node)
    vlocal_ae_node = findfirst("PP_GIPAW_VLOCAL_AE", node)
    vlocal_ae = parse_nodecontent(Float64, vlocal_ae_node)
    vlocal_ps_node = findfirst("PP_GIPAW_VLOCAL_PS", node)
    vlocal_ps = parse_nodecontent(Float64, vlocal_ps_node)
    return UpfGipawVlocal(vlocal_ae, vlocal_ps)
end

function upf2_dump_gipaw_vlocal(vlocal::UpfGipawVlocal)::EzXML.Node
    node = ElementNode("PP_GIPAW_VLOCAL")
    node_ae = ElementNode("PP_GIPAW_VLOCAL_AE")
    text_ae = array_to_text(vlocal.vlocal_ae)
    link!(node_ae, TextNode(text_ae))
    link!(node, node_ae)
    node_ps = ElementNode("PP_GIPAW_VLOCAL_PS")
    text_ps = array_to_text(vlocal.vlocal_ps)
    link!(node_ps, TextNode(text_ps))
    link!(node, node_ps)
    return node
end

function upf2_parse_gipaw(node::EzXML.Node)
    gipaw_data_format = get_attr(Int, node, "gipaw_data_format")
    core_orbitals_node = findfirst("PP_GIPAW_CORE_ORBITALS", node)
    core_orbital_nodes = [n
                          for n in eachnode(core_orbitals_node)
                          if occursin("PP_GIPAW_CORE_ORBITAL.", nodename(n))]
    core_orbitals = upf2_parse_gipaw_core_orbital.(core_orbital_nodes)
    orbitals_node = findfirst("PP_GIPAW_ORBITALS", node)
    if isnothing(orbitals_node)
        orbitals = nothing
    else
        orbital_nodes = [n
                         for n in eachnode(orbitals_node)
                         if occursin("PP_GIPAW_ORBITAL.", nodename(n))]
        orbitals = upf2_parse_gipaw_orbital.(orbital_nodes)
    end
    vlocal_node = findfirst("PP_GIPAW_VLOCAL", node)
    if isnothing(vlocal_node)
        vlocal = nothing
    else
        vlocal = upf2_parse_gipaw_vlocal(vlocal_node)
    end
    return UpfGipaw(gipaw_data_format, core_orbitals, orbitals, vlocal)
end
upf2_parse_gipaw(doc::EzXML.Document) = upf2_parse_gipaw(findfirst("PP_GIPAW", root(doc)))

function upf2_dump_gipaw(gipaw::UpfGipaw)::EzXML.Node
    node = ElementNode("PP_GIPAW")
    set_attr!(node, "gipaw_data_format", gipaw.gipaw_data_format)

    # PP_GIPAW_CORE_ORBITALS
    core_orbitals_node = ElementNode("PP_GIPAW_CORE_ORBITALS")
    set_attr!(core_orbitals_node, "number_of_orbitals", length(gipaw.core_orbitals))
    for core_orbital in gipaw.core_orbitals
        link!(core_orbitals_node, upf2_dump_gipaw_core_orbital(core_orbital))
    end
    link!(node, core_orbitals_node)

    # PP_GIPAW_ORBITALS
    if !isnothing(gipaw.orbitals)
        orbitals_node = ElementNode("PP_GIPAW_ORBITALS")
        for orbital in gipaw.orbitals
            link!(orbitals_node, upf2_dump_gipaw_orbital(orbital))
        end
        link!(node, orbitals_node)
    end

    # PP_GIPAW_VLOCAL
    if !isnothing(gipaw.vlocal)
        link!(node, upf2_dump_gipaw_vlocal(gipaw.vlocal))
    end

    return node
end

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

"""
convert a vector to a string with line break every `width` elements
"""
function array_to_text(v::AbstractVector; width::Int=4)::AbstractString
    s = "\n"
    for i in eachindex(v)
        s *= @sprintf("%20.14E", v[i]) * ' '
        if i % width == 0
            s *= '\n'
        end
    end
    s *= '\n'

    return s
end

"""
convert a matrix to a string
"""
function array_to_text(m::AbstractMatrix)::AbstractString
    s = "\n"
    for row in eachrow(m)
        for el in row
            s *= @sprintf("%20.14E", el) * ' '
        end
        s *= '\n'
    end

    s
end

function set_attr!(node::EzXML.Node, key::String, value::Any)
    # Convert Bool to Fortran-style T/F
    if isa(value, Bool)
        value = value ? "T" : "F"
    end

    # Convert Float64 to Fortran-style D
    if isa(value, Float64)
        value = replace(string(value), "E" => "D")
    end

    # skip if value is nothing
    if !isnothing(value)
        node[key] = value
    end
end

function set_attr!(node::EzXML.Node, key::Symbol, value::Any)
    set_attr!(node, string(key), value)
end

function parse_nodecontent(::Type{T}, node::EzXML.Node) where {T}
    words = split(strip(nodecontent(node)))
    data = Vector{T}(undef, length(words))
    data .= parse.(T, words)
    return data
end
