struct UpfArray{T,N}
    content::Array{T,N}
    size::Int
    type::String
    columns::Int
end

UpfVector{T} = UpfArray{T,1}
UpfMatrix{T} = UpfArray{T,2}

struct UpfPswfc
    chi::Vector{UpfChi}
end

struct UpfChi
    label::String
    l::Int
    content::Vector{Float64}
    pseudo_energy::Float64
    size::Int
    occupation::Float64
    type::String
    columns::Int
    index::Int
    n::Union{Nothing,Int}
    cutoff_radius::Union{Nothing,Float64}
    ultrasoft_cutoff_radius::Union{Nothing,Float64}
end

struct UpfRelWfc
    jchi::Float64
    index::Union{Nothing,Int}
    els::Union{Nothing,String}
    nn::Union{Nothing,Int}
    lchi::Union{Nothing,Int}
    oc::Union{Nothing,Float64}
end

struct UpfRelBeta
    index::Union{Nothing,Int}
    jjj::Float64
    lll::Union{Nothing,Int}
end

struct UpfSpinOrb
    relwfcs::Vector{UpfRelWfc}
    relbetas::Vector{UpfRelBeta}
end

struct UpfWfc
    content::Vector{Float64}
    index::Int
    l::Int
    label::Union{Nothing,String}
end

struct UpfFullWfc
    aewfcs::Vector{UpfChi}
    pswfcs::Vector{UpfChi}
end

struct UpfHeader
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
    "True if data for meta-GGA functionals is present"
    has_metagga::Bool
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

struct UpfMesh
    r::UpfVector{Float64}
    rab::UpfVector{Float64}
    dx::Union{Nothing,Float64}
    mesh::Union{Nothing,Int}
    rmax::Union{Nothing,Float64}
    xmin::Union{Nothing,Float64}
    zmesh::Union{Nothing,Float64}
end

struct UpfInfo
    inputfile::String
end

struct UpfAugmentation
    nqf::Int
    nqlc::Union{Nothing,Int}
    q_with_l::Bool
    qfcoef::Union{Nothing,Vector{UpfQfcoef}}
    rinner::Union{Nothing,UpfRinner}
    qij::Union{Nothing,Vector{UpfQij}}
    qijl::Union{Nothing,Vector{UpfQijl}}
    q::UpfVector{Float64}
end

struct UpfNonLocal
    dij::UpfVector{Float64}
    beta::Vector{UpfBeta}
    augmentation::Union{Nothing,UpfAugmentation}
end

struct UpfBeta
    angular_momentum::Int
    content::Vector{Float64}
    cutoff_radius_index::Int
    cutoff_radius::Float64
    size::Int
    type::String
    columns::Int
    index::Int
    norm_conserving_radius::Union{Nothing,Float64}
    ultrasoft_cutoff_radius::Union{Nothing,Float64}
    label::Union{Nothing,String}
end

struct UpfPaw
    paw_data_format::Union{Nothing,Int}
    core_energy::Union{Nothing,Float64}
    occupations::Vector{Float64}
    ae_nlcc::Vector{Float64}
    ae_vloc::Vector{Float64}
    aewfcs::Vector{UpfWfc}
    pswfcs::Vector{UpfWfc}
end

struct UpfGipawCoreOrbital
    index::Int
    label::Union{Nothing,String}
    n::Int
    l::Int
    content::Vector{Float64}
end

struct UpfGipaw
    gipaw_data_format::Int
    core_orbitals::Vector{UpfGipawCoreOrbital}
end

@serde @de_name @ser_name struct UpfFile
    info::UpfInfo | "info" | "info"
    header::UpfHeader | "header" | "header"
    mesh::UpfMesh | "mesh" | "mesh"
    local_::UpfVector{Float64} | "local" | "local"
    nonlocal::UpfNonLocal | "nonlocal" | "nonlocal"
    nlcc::Union{Nothing,UpfVector{Float64}} | "nlcc" | "nlcc"
    taumod::Union{Nothing,UpfVector{Float64}} | "taumod" | "taumod"
    rhoatom::Union{Nothing,UpfVector{Float64}} | "rhoatom" | "rhoatom"
    spin_orb::Union{Nothing,UpfSpinOrb} | "spin_orb" | "spin_orb"
    paw::Union{Nothing,UpfPaw} | "paw" | "paw"
    gipaw::Union{Nothing,UpfGipaw} | "gipaw" | "gipaw"
    tauatom::Union{Nothing,UpfVector{Float64}} | "tauatom" | "tauatom"
    pswfc::Union{Nothing,UpfPswfc} | "pswfc" | "pswfc"
    full_wfc::Union{Nothing,UpfFullWfc} | "full_wfc" | "full_wfc"
end
