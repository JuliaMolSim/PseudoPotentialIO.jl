abstract type PsmlField end
Serde.deser(::Type{<:PsmlField}, ::Type{Bool}, text::String) = lowercase(strip(text)) == "yes"
function Serde.deser(::Type{<:PsmlField}, ::Type{Vector{Float64}}, text::String)
    return parse.(Float64, strip.(split(text)))
end

struct PsmlInputFile <: PsmlField
    name::String
    _::String
end

struct PsmlProvenance <: PsmlField
    record_number::Union{Int64,Nothing}  # Optional record number
    creator::String
    date::String
    annotation::Union{Dict{String,Any},Nothing}  # Optional annotation
    input_file::Union{Vector{PsmlInputFile},Nothing}  # Zero or more input files
end
Serde.custom_name(::Type{PsmlProvenance}, ::Val{:record_number}) = "record-number"
Serde.custom_name(::Type{PsmlProvenance}, ::Val{:input_file}) = "input-file"
function Serde.deser(::Type{PsmlProvenance}, ::Type{Vector{PsmlInputFile}}, data::Dict{String,Any})
    return [Serde.deser(PsmlInputFile, data)]
end

struct PsmlFunctional <: PsmlField
    id::Int64
    name::String
    weight::Union{Float64,Nothing}
    # exchange | correlation | exchange-correlation |
    # XC_EXCHANGE | XC_CORRELATION | XC_EXCHANGE_CORREALTION
    type::Union{String,Nothing}
end

struct PsmlLibxcInfo <: PsmlField
    number_of_functionals::Int64
    functional::Vector{PsmlFunctional}
end
Serde.custom_name(::Type{PsmlLibxcInfo}, ::Val{:number_of_functionals}) = "number-of-functionals"

struct PsmlExchangeCorrelation <: PsmlField
    annotation::Union{Dict{String,Any},Nothing}
    libxc_info::PsmlLibxcInfo
end
Serde.custom_name(::Type{PsmlExchangeCorrelation}, ::Val{:libxc_info}) = "libxc-info"

struct PsmlShell <: PsmlField
    l::String
    n::Int64
    occupation::Float64
    occupation_up::Union{Float64,Nothing}
    occupation_down::Union{Float64,Nothing}
end
Serde.custom_name(::Type{PsmlShell}, ::Val{:occupation_up}) = "occupation-up"
Serde.custom_name(::Type{PsmlShell}, ::Val{:occupation_down}) = "occupation-down"

struct PsmlValenceConfiguration <: PsmlField
    total_valence_charge::Float64
    annotation::Union{Dict{String,Any},Nothing}
    shell::Vector{PsmlShell}
end
Serde.custom_name(::Type{PsmlValenceConfiguration}, ::Val{:total_valence_charge}) = "total-valence-charge"

struct PsmlCoreConfiguration <: PsmlField
    total_core_charge::Float64
    annotation::Union{Dict{String,Any},Nothing}
    shell::Vector{PsmlShell}
end
Serde.custom_name(::Type{PsmlCoreConfiguration}, ::Val{:total_core_charge}) = "total-core-charge"

struct PsmlPseudoAtomSpec <: PsmlField
    atomic_label::String
    atomic_number::Float64
    z_pseudo::Float64
    core_corrections::Bool
    meta_gga::Union{Bool,Nothing}
    relativity::String
    spin_dft::Union{Bool,Nothing}
    flavor::Union{String,Nothing}
    annotation::Union{Dict{String,Any},Nothing}
    exchange_correlation::PsmlExchangeCorrelation
    valence_configuration::PsmlValenceConfiguration
    core_configuration::Union{PsmlCoreConfiguration,Nothing}
end
Serde.custom_name(::Type{PsmlPseudoAtomSpec}, ::Val{:atomic_label}) = "atomic-label"
Serde.custom_name(::Type{PsmlPseudoAtomSpec}, ::Val{:atomic_number}) = "atomic-number"
Serde.custom_name(::Type{PsmlPseudoAtomSpec}, ::Val{:z_pseudo}) = "z-pseudo"
Serde.custom_name(::Type{PsmlPseudoAtomSpec}, ::Val{:core_corrections}) = "core-corrections"
Serde.custom_name(::Type{PsmlPseudoAtomSpec}, ::Val{:meta_gga}) = "meta-gga"
Serde.custom_name(::Type{PsmlPseudoAtomSpec}, ::Val{:spin_dft}) = "spin-dft"
Serde.custom_name(::Type{PsmlPseudoAtomSpec}, ::Val{:exchange_correlation}) = "exchange-correlation"
Serde.custom_name(::Type{PsmlPseudoAtomSpec}, ::Val{:valence_configuration}) = "valence-configuration"
Serde.custom_name(::Type{PsmlPseudoAtomSpec}, ::Val{:core_configuration}) = "core-configuration"

struct PsmlGridData <: PsmlField
    _::Vector{Float64}
end

struct PsmlGrid <: PsmlField
    npts::Int64
    annotation::Dict{String,Any}
    grid_data::PsmlGridData
end
Serde.custom_name(::Type{PsmlGrid}, ::Val{:grid_data}) = "grid-data"

struct PsmlData <: PsmlField
    npts::Union{Int64,Nothing}
    _::Vector{Float64}
end

struct PsmlRadfunc <: PsmlField
    grid::Union{PsmlGrid,Nothing}
    data::PsmlData
end

struct PsmlValenceCharge <: PsmlField
    total_charge::Float64
    is_unscreening_charge::Union{Bool,Nothing}
    rescaled_to_z_pseudo::Union{Bool,Nothing}
    annotation::Union{Dict{String,Any},Nothing}
    radfunc::PsmlRadfunc
end
Serde.custom_name(::Type{PsmlValenceCharge}, ::Val{:total_charge}) = "total-charge"
Serde.custom_name(::Type{PsmlValenceCharge}, ::Val{:is_unscreening_charge}) = "is-unscreening-charge"
Serde.custom_name(::Type{PsmlValenceCharge}, ::Val{:rescaled_to_z_pseudo}) = "rescaled-to-z-pseudo"

struct PsmlCoreCharge <: PsmlField
    matching_radius::Union{Float64,Nothing}
    number_of_continuous_derivatives::Int64
    annotation::Union{Dict{String,Any},Nothing}
    radfunc::PsmlRadfunc
end
Serde.custom_name(::Type{PsmlCoreCharge}, ::Val{:matching_radius}) = "matching-radius"
Serde.custom_name(::Type{PsmlCoreCharge}, ::Val{:number_of_continuous_derivatives}) = "number-of-continuous-derivatives"

struct PsmlValenceKineticDensity <: PsmlField
    is_unscreening_tau::Union{Bool,Nothing}
    annotation::Union{Dict{String,Any},Nothing}
    radfunc::PsmlRadfunc
end
Serde.custom_name(::Type{PsmlValenceKineticDensity}, ::Val{:is_unscreening_tau}) = "is-unscreening-tau"

struct PsmlCoreKineticDensity <: PsmlField
    matching_radius::Union{Float64,Nothing}
    number_of_continuous_derivatives::Int64
    annotation::Union{Dict{String,Any},Nothing}
    radfunc::PsmlRadfunc
end
Serde.custom_name(::Type{PsmlCoreKineticDensity}, ::Val{:number_of_continuous_derivatives}) = "number-of-continuous-derivatives"

struct PsmlPotential <: PsmlField
    flavor::Union{String,Nothing}
    l::String
    j::Union{Float64,Nothing}
    n::Int64
    rc::Float64
    eref::Union{Float64,Nothing}
    radfunc::PsmlRadfunc
end

struct PsmlSemilocalPotentials <: PsmlField
    set::String
    flavor::Union{String,Nothing}
    annotation::Union{Dict{String,Any},Nothing}
    grid::Union{PsmlGrid,Nothing}
    slps::Vector{PsmlPotential}
end

struct PsmlLocalCharge <: PsmlField
    radfunc::PsmlRadfunc
end

struct PsmlLocalPotential <: PsmlField
    type::String
    annotation::Union{Dict{String,Any},Nothing}
    grid::Union{PsmlGrid,Nothing}
    radfunc::PsmlRadfunc
    local_charge::Union{PsmlLocalCharge,Nothing}
end
Serde.custom_name(::Type{PsmlLocalPotential}, ::Val{:local_charge}) = "local-charge"

struct PsmlProjector <: PsmlField
    ekb::Float64
    eref::Union{Float64,Nothing}
    l::String
    j::Union{Float64,Nothing}
    seq::Int64
    type::String
    radfunc::PsmlRadfunc
end

struct PsmlNonlocalProjectors <: PsmlField
    set::String
    annotation::Union{Dict{String,Any},Nothing}
    grid::Union{PsmlGrid,Nothing}
    proj::Vector{PsmlProjector}
end

struct PsmlPseudoWf <: PsmlField
    l::String
    j::Union{Float64,Nothing}
    n::Int64
    energy_level::Union{Float64,Nothing}
    radfunc::PsmlRadfunc
end
Serde.custom_name(::Type{PsmlPseudoWf}, ::Val{:energy_level}) = "energy-level"

struct PsmlPseudoWaveFunctions <: PsmlField
    set::String
    annotation::Union{Dict{String,Any},Nothing}
    grid::Union{PsmlGrid,Nothing}
    pswf::Vector{PsmlPseudoWf}
end

struct PsmlFile <: PsPFile
    version::String
    energy_unit::String
    length_unit::String
    uuid::String
    # One or more provenance elements
    provenance::Vector{PsmlProvenance}
    pseudo_atom_spec::PsmlPseudoAtomSpec
    # Optional top-level grid
    grid::Union{PsmlGrid,Nothing}
    valence_charge::PsmlValenceCharge
    # Optional pseudo-core charge
    pseudocore_charge::Union{PsmlCoreCharge,Nothing}
    # Optional valence-kinetic energy density for MGGA
    valence_kinetic_density::Union{PsmlValenceKineticDensity,Nothing}
    # Optional core-kinetic energy density for MGGA
    core_kinetic_density::Union{PsmlCoreKineticDensity,Nothing}
    # One or more semilocal potential groups and optionally (local potential and zero or more fully nonlocal groups)
    # OR zero or more semilocal potential groups and (local potential and zero or more fully nonlocal groups)
    semilocal_potentials::Union{Vector{PsmlSemilocalPotentials},Nothing}
    local_potential::Union{PsmlLocalPotential,Nothing}
    nonlocal_projectors::Union{Vector{PsmlNonlocalProjectors},Nothing}
    # Zero or more Pseudo Wavefunction groups
    pseudo_wave_functions::Union{Vector{PsmlPseudoWaveFunctions},Nothing}
end
Serde.custom_name(::Type{PsmlFile}, ::Val{:pseudo_atom_spec}) = "pseudo-atom-spec"
Serde.custom_name(::Type{PsmlFile}, ::Val{:valence_charge}) = "valence-charge"
Serde.custom_name(::Type{PsmlFile}, ::Val{:pseudocore_charge}) = "pseudocore-charge"
Serde.custom_name(::Type{PsmlFile}, ::Val{:valence_kinetic_density}) = "valence-kinetic-density"
Serde.custom_name(::Type{PsmlFile}, ::Val{:core_kinetic_density}) = "core-kinetic-density"
Serde.custom_name(::Type{PsmlFile}, ::Val{:semilocal_potentials}) = "semilocal-potentials"
Serde.custom_name(::Type{PsmlFile}, ::Val{:local_potential}) = "local-potential"
Serde.custom_name(::Type{PsmlFile}, ::Val{:nonlocal_projectors}) = "nonlocal-projectors"
function Serde.deser(::Type{PsmlFile}, ::Type{Vector{PsmlProvenance}}, data::Dict{String,Any})
    return [Serde.deser(PsmlProvenance, data)]
end
function Serde.deser(::Type{PsmlFile}, ::Type{Vector{PsmlSemilocalPotentials}}, data::Dict{String,Any})
    return [Serde.deser(PsmlSemilocalPotentials, data)]
end
function Serde.deser(::Type{PsmlFile}, ::Type{Vector{PsmlNonlocalProjectors}}, data::Dict{String,Any})
    return [Serde.deser(PsmlNonlocalProjectors, data)]
end
function Serde.deser(::Type{PsmlFile}, ::Type{Vector{PsmlPseudoWaveFunctions}}, data::Dict{String,Any})
    return [Serde.deser(PsmlPseudoWaveFunctions, data)]
end

identifier(file::PsmlFile) = file.uuid
format(::PsmlFile) = "PSML"
element(file::PsmlFile) = file.pseudo_atom_spec.atomic_label
functional(file::PsmlFile) = [Symbol(LIBXC_FUNCTIONALS_BY_ID[func.id]) for func in file.pseudo_atom_spec.exchange_correlation.libxc_info.functional]
valence_charge(file::PsmlFile) = file.pseudo_atom_spec.valence_configuration.total_valence_charge
is_norm_conserving(::PsmlFile) = true
is_ultrasoft(::PsmlFile) = false
is_paw(::PsmlFile) = false
has_spin_orbit(file::PsmlFile) = file.pseudo_atom_spec.relativity == "dirac"
has_model_core_charge_density(file::PsmlFile) = file.pseudo_atom_spec.core_corrections
