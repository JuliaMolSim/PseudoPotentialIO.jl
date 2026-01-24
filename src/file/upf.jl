"""
UPF `<PP_HEADER>`
"""
struct UpfHeader <: PsPFile
    "Generation code"
    generated::Union{Nothing,String}
    author::Union{Nothing,String}
    "Generation date"
    date::Union{Nothing,String}
    comment::Union{Nothing,String}
    "A valid chemical symbol: `{H, He, Li, ..., Og}`"
    element::String
    "A valid type of pseudopotential: `{NC, SL, 1/r, US, PAW, USPP}`"
    pseudo_type::String
    "A valid relativistic treatment: `{scalar, full, relativistic}`"
    relativistic::Union{Nothing,String}
    is_ultrasoft::Bool
    is_paw::Bool
    "True of the pseudopotential is just a bare Coulomb potential (all-electron)"
    is_coulomb::Bool
    "True if fully-relativistic with spin-orbit terms"
    has_so::Bool
    "True if pseudo-atomic wavefunctions present"
    has_wfc::Bool
    "True if data for GIPAW reconstruction is present"
    has_gipaw::Bool
    "True if data for GIPAW reconstruction is present"
    paw_as_gipaw::Union{Nothing,Bool}
    "True if non-linear core correction is included"
    core_correction::Bool
    "QuantumEspresso exchange-correlation identifiers"
    functional::String
    "Pseudo-atomic charge"
    z_valence::Float64
    "Total pseudo-valence energy of the pseudopotential"
    total_psenergy::Union{Nothing,Float64}
    "Suggested plane wave cutoff for expansion of Kohn-Sham orbitals"
    wfc_cutoff::Union{Nothing,Float64}
    "Suggested plane wave cutoff for expansion of charge density"
    rho_cutoff::Union{Nothing,Float64}
    "Maximum angular momentum channel in the pseudopotential"
    l_max::Int
    "Maximum angular momentum channel in the atomic charge density (PAW only)"
    l_max_rho::Union{Nothing,Int}
    "Angular momentum chosen to be the local potential (-1 if none)"
    l_local::Union{Nothing,Int}
    "Number of points in the radial grid"
    mesh_size::Int
    "Number of chi-functions"
    number_of_wfc::Int
    "Number of Kleinman-Bylander nonlocal projectors"
    number_of_proj::Int
end

function Base.show(io::IO, ::MIME"text/plain", header::UpfHeader)
    println(io, typeof(header))
    if !isnothing(header.generated)
        @printf "%032s: %s\n" "generated" header.generated
    end
    if !isnothing(header.author)
        @printf "%032s: %s\n" "author" header.author
    end
    if !isnothing(header.date)
        @printf "%032s: %s\n" "date" header.date
    end
    @printf "%032s: %s\n" "element" header.element
    @printf "%032s: %s\n" "pseudopotential type" header.pseudo_type
    if !isnothing(header.relativistic)
        @printf "%032s: %s\n" "relativistic treatment" header.relativistic
    end
    @printf "%032s: %s\n" "is ultrasoft" header.is_ultrasoft
    @printf "%032s: %s\n" "is PAW" header.is_paw
    @printf "%032s: %s\n" "is bare Coulomb" header.is_coulomb
    @printf "%032s: %s\n" "has spin-orbit" header.has_so
    @printf "%032s: %s\n" "has pseudo-atomic wavefunctions" header.has_wfc
    @printf "%032s: %s\n" "has GIPAW data" header.has_gipaw
    if !isnothing(header.paw_as_gipaw)
        @printf "%032s: %s\n" "PAW as GIPAW" header.paw_as_gipaw
    end
    @printf "%032s: %s\n" "has non-linear core correction" header.core_correction
    @printf "%032s: %s\n" "functional (QuantumESPRESSO)" header.functional
    @printf "%032s: %.6f\n" "valence charge" header.z_valence
    if !isnothing(header.total_psenergy)
        @printf "%032s: %.6f\n" "total pseudo energy" header.total_psenergy
    end
    if !isnothing(header.wfc_cutoff)
        @printf "%032s: %.6f Ry\n" "suggested wavefunction cutoff" header.wfc_cutoff
    end
    if !isnothing(header.rho_cutoff)
        @printf "%032s: %.6f Ry\n" "suggested charge density cutoff" header.rho_cutoff
    end
    @printf "%032s: %d\n" "max angular momentum" header.l_max
    if !isnothing(header.l_max_rho)
        @printf "%032s: %d\n" "max angular momentum in charge density (PAW)" header.l_max_rho
    end
    @printf "%032s: %d\n" "local potential angular momentum" header.l_local
    @printf "%032s: %d\n" "radial mesh size" header.mesh_size
    @printf "%032s: %d\n" "number of pseudo wavefunctions" header.number_of_wfc
    @printf "%032s: %d\n" "number of KB projectors" header.number_of_proj
end

"""
UPF `<PP_MESH>`
"""
struct UpfMesh
    "Radial mesh"
    r::Vector{Float64}
    "Integration factors for integrating quantities on the radial mesh"
    rab::Vector{Float64}
    "Number of points in the radial mesh"
    mesh::Union{Nothing,Int}
    "Maximum value of the radial mesh"
    rmax::Union{Nothing,Float64}
    # Mesh generation parameters
    dx::Union{Nothing,Float64}
    xmin::Union{Nothing,Float64}
    zmesh::Union{Nothing,Float64}
end

"""
UPF `<PP_NONLOCAL/PP_AUGMENTATION/PP_QIJ.[i].[j]>`
"""
struct UpfQij
    qij::Vector{Float64}
    first_index::Union{Nothing,Int}
    second_index::Union{Nothing,Int}
    composite_index::Union{Nothing,Int}
    is_null::Union{Nothing,Bool}
end

"""
UPF `<PP_NONLOCAL/PP_AUGMENTATION/PP_QIJL.[i].[j].[l]>`
"""
struct UpfQijl
    qijl::Vector{Float64}
    angular_momentum::Int
    first_index::Union{Nothing,Int}
    second_index::Union{Nothing,Int}
    composite_index::Union{Nothing,Int}
    is_null::Union{Nothing,Bool}
end

struct UpfQfcoef
    qfcoef::Vector{Float64}
    angular_momentum::Union{Nothing,Int}
    first_index::Union{Nothing,Int}
    second_index::Union{Nothing,Int}
end

"""
UPF `<PP_NONLOCAL/PP_AUGMENTATION>`
"""
struct UpfAugmentation
    """Integrals of the augmentation functions 4π ∫ Qij(r) r^2 dr.
    NB: `q = 0` does _not_ guarantee that the corresponding augmentation function is
    zero."""
    q::Matrix{Float64}
    """(PAW) Electronic multipoles of the corresponding augmentation channel. If the
    absolute value of a multipole is less than `augmentation_epsilon`, the corresponding
    augmentation function should be considered zero"""
    multipoles::Union{Nothing,Vector{Float64}}
    """Coefficients used to perform a Taylor expansion of the augmentation functions at
    small radii (NB: compulsory if `nqf > 0`, ignored otherwise)"""
    qfcoefs::Union{Nothing,Vector{UpfQfcoef}}
    rinner::Union{Nothing,Vector{Float64}}
    "If `q_with_l` is false, the augmentation functions for `i,j in 1:number_of_proj`"
    qijs::Union{Nothing,Vector{UpfQij}}
    """If `q_with_l` is true, the angular-momentum dependent augmentation functions for
    `i,j in 1:number_of_proj` and `l in 0:l_max_rho`"""
    qijls::Union{Nothing,Vector{UpfQijl}}
    "True if augmentation charge functions are decomposed into angular momentum components"
    q_with_l::Bool
    """Number of expansion coefficients for analytical expansion of the augmentation
    charge at small radius."""
    nqf::Int
    "Number of angular momenta terms in the augmentation charge, unused if `nqf = 0`"
    nqlc::Union{Nothing,Float64}
    """(UNUSED) (PAW) Shape of the augmentation function: `{PSQ, GAUSS, BESSEL}`, could
    be used for analyical overlap of PAW augmentation charge"""
    shape::Union{Nothing,String}
    "(PAW) Radial grid index beyond which augmentation charge is zero"
    iraug::Union{Nothing,Int}
    "(PAW) Radial distance beyond which augmentation charge is zero"
    raug::Union{Nothing,Float64}
    "(PAW): Maximum angular momentum appearing in augmentation charge"
    l_max_aug::Union{Nothing,Float64}
    """(PAW): Augmentation functions whose norms are less than `augmentation_epsilon` are
    considered zero"""
    augmentation_epsilon::Union{Nothing,Float64}
    "(DEPRECATED?)"
    cutoff_r::Union{Nothing,Float64}
    "(DEPRECATED?)"
    cutoff_r_index::Union{Nothing,Float64}
end

"""
UPF `<PP_NONLOCAL/PP_BETA.[i]>`
"""
struct UpfBeta
    "Kleinman-Bylander nonlocal projector multiplied by the radial mesh, on the radial mesh"
    beta::Vector{Float64}  # misnomer: actually r*β(r)!
    index::Union{Nothing,Int}
    angular_momentum::Int
    cutoff_radius_index::Union{Nothing,Int}
    cutoff_radius::Union{Nothing,Float64}
    norm_conserving_radius::Union{Nothing,Float64}
    ultrasoft_cutoff_radius::Union{Nothing,Float64}
    label::Union{Nothing,String}
end

"""
UPF `<PP_NONLOCAL>`
"""
struct UpfNonlocal
    """Kleinman-Bylander nonlocal projectors"""
    betas::Vector{UpfBeta}
    "Kleinman-Bylander energies"
    dij::Matrix{Float64}
    "Agumentation data for ultrasoft and PAW pseudopotentials"
    augmentation::Union{Nothing,UpfAugmentation}
end

"""
UPF `<PP_PSWFC/PP_CHI>`
"""
struct UpfChi
    "Pseudo-atomic valence wavefunction on the radial mesh"
    chi::Vector{Float64}  # misnomer: actually r*χ(r)!
    "Angular momentum"
    l::Int
    occupation::Float64
    index::Union{Nothing,Int}
    label::Union{Nothing,String}
    "Principal quantum number"
    n::Union{Nothing,Int}
    pseudo_energy::Union{Nothing,Float64}
    cutoff_radius::Union{Nothing,Float64}
    ultrasoft_cutoff_radius::Union{Nothing,Float64}
end

"""
UPF `<PP_SPIN_ORB/PP_RELWFC.[i]>`
"""
struct UpfRelWfc
    jchi::Float64
    index::Union{Nothing,Int}
    els::Union{Nothing,String}
    nn::Union{Nothing,Int}
    lchi::Union{Nothing,Int}
    oc::Union{Nothing,Float64}
end

"""
UPF `<PP_SPIN_ORB/PP_RELBETA.[i]>`
"""
struct UpfRelBeta
    index::Union{Nothing,Int}
    jjj::Float64
    lll::Union{Nothing,Int}
end

"""
UPF `<PP_SPIN_ORB>`
"""
struct UpfSpinOrb
    relwfcs::Vector{UpfRelWfc}
    relbetas::Vector{UpfRelBeta}
end

"""
UPF `<//PP_(AE|PS)WFC.[i]>`
"""
struct UpfWfc
    wfc::Vector{Float64}
    index::Int
    l::Int
    label::Union{Nothing,String}
end

"""
UPF `<PP_FULL_WFC>`
"""
struct UpfFullWfc
    aewfcs::Vector{UpfWfc}
    pswfcs::Vector{UpfWfc}
end

"""
UPF `<PP_PAW>`
"""
struct UpfPaw
    paw_data_format::Union{Nothing,Int}
    core_energy::Union{Nothing,Float64}
    occupations::Vector{Float64}
    ae_nlcc::Vector{Float64}
    ae_vloc::Vector{Float64}
    aewfcs::Vector{UpfWfc}
    pswfcs::Vector{UpfWfc}
end

"""
UPF `<PP_GIPAW/PP_GIPAW_CORE_ORBITALS/PP_GIPAW_CORE_ORBITAL.[i]>`
"""
struct UpfGipawCoreOrbital
    index::Int
    label::Union{Nothing,String}
    "Principal quantum number"
    n::Int
    "Angular momentum"
    l::Int
    core_orbital::Vector{Float64}
end

"""
UPF `<PP_GIPAW>`
"""
struct UpfGipaw
    gipaw_data_format::Int
    core_orbitals::Vector{UpfGipawCoreOrbital}
end

"""
Universal Pseudopotential Format file contents.
"""
struct UpfFile <: PsPFile
    "Identifier"
    identifier::String
    "UPF format version"
    version::String
    "Optional general information about the pseudopotential, often generation input"
    info::Union{Nothing,String}
    "Various pseudopotential metadata"
    header::UpfHeader
    "Radial mesh, mesh integration factors, and other mesh information"
    mesh::UpfMesh
    "Pseudized core charge on the radial grid, (ignored if `core_correction` is false)"
    nlcc::Union{Nothing,Vector{Float64}}
    "Local part of the pseudopotential on the radial grid (ignored if `is_coulomb`)"
    local_::Union{Nothing,Vector{Float64}}
    "Nonlocal part of the pseudopotential"
    nonlocal::UpfNonlocal
    "Pseudo-atomic valence wavefunctions"
    pswfc::Vector{UpfChi}
    "All-electron wavefunctions"
    full_wfc::Union{Nothing,UpfFullWfc}
    "Pseudo-atomic valence charge density on the radial grid with prefactor 4πr²"
    rhoatom::Vector{Float64}
    "Spin-orbit coupling data, (ignored if `has_so` is false)"
    spin_orb::Union{Nothing,UpfSpinOrb}
    "PAW data, (ignored if `is_paw` is false)"
    paw::Union{Nothing,UpfPaw}
    "GIPAW data"
    gipaw::Union{Nothing,UpfGipaw}
end

function UpfFile(path::AbstractString; identifier="")
    identifier = isempty(identifier) ? splitpath(path)[end] : identifier
    open(path, "r") do io
        return UpfFile(io; identifier)
    end
end

function UpfFile(io::IO; identifier="")
    version = _get_upf_version(io)
    if version == 1
        return upf1_parse_psp(io; identifier)
    end
    if version == 2
        return upf2_parse_psp(io; identifier)
    end
    error("Unknown UPF version.")
end

"""
Write UPF file to `io`. Note that the file is written the UPF version passed
as a string and this value takes preference over the version stored in `psp`.
"""
function save_psp(io::IO, psp::UpfFile; format="UPF v2.0.1")
    if format == "UPF v2.0.1"
        root_node = upf2_dump_psp(psp)
        prettyprint(io, root_node)
    elseif format == "UPF v1.old"
        error("$format not supported yet.")
    else
        error("$format is unknown. Only 'UPF v2.0.1' and 'UPF v1.old' are known.")
    end
end

function _get_upf_version(io::IO)::Int
    pos = position(io)
    line1 = readline(io)
    line2 = readline(io)
    seek(io, pos)
    if occursin("<PP_INFO>", line1)
        # Old UPF files start with the `<PP_INFO>` section
        return 1
    elseif occursin("UPF version=\"2.0.1\"", line1)
        # New UPF files with schema are in XML and start with a version tag
        return 2
    elseif occursin("xml version=\"1.0\"", line1) && occursin("UPF version=\"2.0.1\"", line2)
        # dumped UPF files start with a xml declaration
        return 2
    else
        error("Unknown UPF version")
    end
end

function _get_upf_version(path::AbstractString)::Int
    open(path, "r") do io
        return _get_upf_version(io)
    end
end

identifier(psp::UpfFile)::String = psp.identifier
format(file::UpfFile)::String = "UPF v$(file.version)"
function libxc_string(header::UpfHeader)
    functional = header.functional
    upf_codes = lowercase.(split(functional))

    if length(upf_codes) == 1  # Short code
        entry_index = findfirst(e -> e["name"] == upf_codes[1], UPF_SHORT_NAMES)
        long_code = UPF_SHORT_NAMES[entry_index]["full_name"]
        upf_codes = split(long_code, '+')
    end

    @assert length(upf_codes) == 4

    if upf_codes[3] == "nogx"
        exc_index = findfirst(e -> e["name"] == upf_codes[1], UPF_LDA_EXCHANGE)
        exc = UPF_LDA_EXCHANGE[exc_index]["libxc"]
    else
        exc_index = findfirst(e -> e["name"] == upf_codes[3], UPF_GGA_EXCHANGE)
        exc = UPF_GGA_EXCHANGE[exc_index]["libxc"]
    end

    if upf_codes[4] == "nogc"
        corr_index = findfirst(e -> e["name"] == upf_codes[2], UPF_LDA_CORRELATION)
        corr = UPF_LDA_CORRELATION[corr_index]["libxc"]
    else
        corr_index = findfirst(e -> e["name"] == upf_codes[4], UPF_GGA_CORRELATION)
        corr = UPF_GGA_CORRELATION[corr_index]["libxc"]
    end

    return join([exc, corr], ' ')
end
functional(file::UpfFile)::Vector{Symbol} = map(s -> Symbol(lowercase(s)), split(libxc_string(file.header)))
element(file::UpfFile)::PeriodicTable.Element = PeriodicTable.elements[Symbol(file.header.element)]
is_norm_conserving(file::UpfFile)::Bool = file.header.pseudo_type == "NC"
is_ultrasoft(file::UpfFile)::Bool = file.header.pseudo_type in ("US", "USPP")
is_paw(file::UpfFile)::Bool = file.header.pseudo_type == "PAW"
has_spin_orbit(file::UpfFile)::Bool = file.header.has_so
has_model_core_charge_density(file::UpfFile)::Bool = file.header.core_correction
valence_charge(file::UpfFile)::Int = file.header.z_valence

function libxc_to_qe_libxc(libxc_string::AbstractString)::String
    # https://gitlab.com/QEF/q-e/-/blob/develop/Modules/funct.f90?ref_type=heads#L301
    # ! NOTE FOR LIBXC USERS: to use libxc functionals you must enforce them from input (use
    # ! 'input_dft' in &system) and write their IDs in the input string. The only notation
    # ! now allowed (v7.0) for input DFTs containing Libxc terms is:
    # ! XC-000i-000i-000i-000i-000i-000i
    # ! where you put the functional IDs instead of the zeros and an 'L' instead of
    # ! 'i' if the functional is from Libxc. The order is the usual one:
    # ! LDAexch - LDAcorr - GGAexch - GGAcorr - MGGAexch - MGGAcorr
    # ! however QE will automatically adjust it if needed. You can skip zero tails (e.g.
    # ! you don't need GGA/MGGA slots if the dft is LDA only and so on.
    # ! You can use combinations of qe and libxc functionals, when they are compatible.
    # ! You can also add vdW terms after it, for example, sla+pw+rw86+vdw2 is:
    # ! input_dft='XC-001i-004i-013i-vdw2'.
    # ! For more details see the user_guide (in 'Doc' folder).
    # https://www.quantum-espresso.org/wp-content/uploads/2022/03/user_guide.pdf
    # Some functionals in libxc incorporate the exchange part and the correlation one into one term
    # only (e.g. the ones that include the ‘ xc’ kind label in their name). In these cases the whole
    # functional is formally treated as ‘correlation only’ by Quantum ESPRESSO. This does not
    # imply any loss of information.
    terms = zeros(Int, 6)
    for substring in split(libxc_string)
        substring_parts = split(substring, '_')
        if substring_parts[1] == "LDA"
            if substring_parts[2] == "X"
                terms[1] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            elseif substring_parts[2] in ("C", "XC")
                terms[2] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            else
                error("Unsupported LibXC functional $(libxc_string)")
            end
        elseif substring_parts[1] == "GGA"
            if substring_parts[2] == "X"
                terms[3] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            elseif substring_parts[2] in ("C", "XC")
                terms[4] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            else
                error("Unsupported LibXC functional $(libxc_string)")
            end
        elseif subsrting_parts[1] == "MGGA"
            if substring_parts[2] == "X"
                terms[5] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            elseif substring_parts[2] in ("C", "XC")
                terms[6] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            else
                error("Unsupported LibXC functional $(libxc_string)")
            end
        else
            error("Unsupported LibXC functional $(libxc_string)")
        end
    end
    return "XC" * prod(terms) do term
        iszero(term) && return @sprintf "-%03dI" term
        return @sprintf "-%0dL" term
    end
end

function libxc_to_qe_old(libxc_string::AbstractString, sep::Char=' ')::String
    luts = Dict(
        "LDA" => Dict("X" => UPF_LDA_EXCHANGE, "C" => UPF_LDA_CORRELATION),
        "GGA" => Dict("X" => UPF_GGA_EXCHANGE, "C" => UPF_GGA_CORRELATION),
    )
    labels_indices = Dict(
        "LDA" => Dict("X" => 1, "C" => 2),
        "GGA" => Dict("X" => 3, "C" => 4)
    )
    qe_labels = ["nox", "noc", "nogx", "nogc"]
    for libxc_name in split(libxc_string)
        family, x_or_c = split(libxc_name, '_')[1:2]

        family_lut = get(luts, family, nothing)
        isnothing(family_lut) && error("Could not convert $(libxc_string) to QE complete name.")

        lut = get(family_lut, x_or_c, nothing)
        isnothing(lut) && return error("Could not convert $(libxc_string) to QE complete name.")

        labels_idx = labels_indices[family][x_or_c]

        lut_idx = findfirst(entry -> entry["libxc"] == libxc_name, lut)
        qe_labels[labels_idx] = isnothing(lut_idx) ? nothing : lut[lut_idx]["name"]
    end
    # LibXC GGA functionals include LDA exchange and correlation: Slater exchange, Perdew-Wang correlation.
    # Lee-Yang-Parr functionals are an exception, but they don't have QE old equivalents in the lookup tables,
    # so we don't handle them here.
    if qe_labels[3] != "nogx"

        qe_labels[1] = "sla"
    end
    if qe_labels[4] != "nogc"
        qe_labels[2] = "pw"
    end
    any(isnothing.(qe_labels)) && error("Could not convert $(libxc_string) to QE complete name.")
    return join(qe_labels, sep)
end
