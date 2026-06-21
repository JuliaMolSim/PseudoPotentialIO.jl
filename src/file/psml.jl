using Serde.SerXml

"""
Abstract type for PSML fields. Used for `Bool` and `Vector{Float64}` (de)serialization.
"""
abstract type PsmlField end

# "yes" | "no" -> true | false
function Serde.deser(::Type{<:PsmlField}, ::Type{Bool}, text::String)::Bool
    if text == "yes"
        return true
    elseif text == "no"
        return false
    else
        error("Invalid serialized boolean value $(text)")
    end
end

# Deserialize numeric vectors
function Serde.deser(::Type{<:PsmlField}, ::Type{Vector{T}}, text::String)::Vector{T} where {T<:Number}
    return parse.(T, strip.(split(text)))
end

# true | false -> "yes" | "no"
function Serde.ser_type(::Type{T}, b::Bool) where {T<:PsmlField}
    b ? "yes" : "no"
end

# Monkey patch to properly serialize the quantity arrays
function Serde.SerXml._node_content(node::T)::String where {T<:PsmlField}
    if hasfield(T, Symbol(Serde.SerXml.CONTENT_WORD))
        value = getfield(node, Symbol(Serde.SerXml.CONTENT_WORD))
        io = IOBuffer()
        write(io, "\n")
        for (i, v) in enumerate(value)
            @printf io " %19.12E" v
            mod(i, 4) == 0 && write(io, "\n")
        end
        write(io, "\n")
        seekstart(io)
        return read(io, String)
    end
    return ""
end

# Monkey patch to escape string XML attributes
function Serde.SerXml._xml_attributes_string(node::T) where {T<:PsmlField}
    attributes = Serde.SerXml._xml_attributes(node)
    attributes = Dict(
        n => typeof(v) <: AbstractString ? escape_for_xml(v) : v for (n, v) in attributes
    )
    return join([" $n=\"$v\"" for (n, v) in attributes])
end

function Base.:(==)(f1::T, f2::T) where {T<:PsmlField}
    for fn in fieldnames(T)
        x1 = getfield(f1, fn)
        x2 = getfield(f2, fn)
        if typeof(x1) <: AbstractDict && typeof(x2) <: AbstractDict
            for (k1, v1) in x1
                !haskey(x2, k1) && return false
                v1 != x2[k1] && return false
            end
        elseif typeof(x1) <: AbstractArray && typeof(x2) <: AbstractArray
            !all(x1 .== x2) && return false
        else
            x1 != x2 && return false
        end
    end
    return true
end

"""
Input file text container.
"""
struct PsmlInputFile <: PsmlField
    "No spaces or commas allowed"
    name::String
    _::String
end
# More pirating to wrap the input file in a CDATA tag
# so that we don't have to explicitly escape its contents.
function Serde.SerXml._node_content(node::PsmlInputFile)::String
    "<![CDATA[" * node._ * "]]>"
end

"""
Provenance entry.
"""
struct PsmlProvenance <: PsmlField
    "Optional record number"
    record_number::Union{Int64, Nothing}
    creator::String
    date::String
    "Optional annotation"
    annotation::Union{Dict{String, Any}, Nothing}
    "Zero or more input files"
    input_file::Union{Vector{PsmlInputFile}, Nothing}
end
function Serde.deser(::Type{PsmlProvenance}, ::Type{Vector{PsmlInputFile}}, data::Dict{String, Any})
    return [Serde.deser(PsmlInputFile, data)]
end

"""
LibXC-compatible functional description.
"""
struct PsmlFunctional <: PsmlField
    "LibXC ID"
    id::Int64
    "Functional name (no convention)"
    name::String
    "Optional weight coefficient"
    weight::Union{Float64, Nothing}
    "Optional canonical name or LibXC symbol: exchange | correlation | exchange-correlation | XC_EXCHANGE | XC_CORRELATION | XC_EXCHANGE_CORRELATION"
    type::Union{String, Nothing}
end

"""
LibXC-compatible functional container.
"""
struct PsmlLibxcInfo <: PsmlField
    number_of_functionals::Int64
    functional::Vector{PsmlFunctional}
end

"""
Exchange-correlation container.
"""
struct PsmlExchangeCorrelation <: PsmlField
    annotation::Union{Dict{String, Any}, Nothing}
    libxc_info::PsmlLibxcInfo
end

"""
Atomic shell description.
"""
struct PsmlShell <: PsmlField
    l::String
    n::Int64
    occupation::Float64
    occupation_up::Union{Float64, Nothing}
    occupation_down::Union{Float64, Nothing}
end

"""
Valence configuration.
"""
struct PsmlValenceConfiguration <: PsmlField
    total_valence_charge::Float64
    annotation::Union{Dict{String, Any}, Nothing}
    shell::Vector{PsmlShell}
end

"""
Core configuration.
"""
struct PsmlCoreConfiguration <: PsmlField
    total_core_charge::Float64
    annotation::Union{Dict{String, Any}, Nothing}
    shell::Vector{PsmlShell}
end

"""
Pseudo-atom specification.
"""
struct PsmlPseudoAtomSpec <: PsmlField
    atomic_label::String
    atomic_number::Float64
    z_pseudo::Float64
    core_corrections::Bool
    meta_gga::Union{Bool, Nothing}
    relativity::String
    spin_dft::Union{Bool, Nothing}
    flavor::Union{String, Nothing}
    annotation::Union{Dict{String, Any}, Nothing}
    exchange_correlation::PsmlExchangeCorrelation
    valence_configuration::PsmlValenceConfiguration
    core_configuration::Union{PsmlCoreConfiguration, Nothing}
end

"""
Grid vector container.
"""
struct PsmlGridData <: PsmlField
    _::Vector{Float64}
end
Base.:(==)(gd1::PsmlGridData, gd2::PsmlGridData) = all(gd1._ .== gd2._)

"""
Grid information
"""
struct PsmlGrid <: PsmlField
    npts::Int64
    annotation::Dict{String, Any}
    grid_data::PsmlGridData
end

"""
Radial function container.
"""
struct PsmlData <: PsmlField
    npts::Union{Int64, Nothing}
    _::Vector{Float64}
end
Base.:(==)(d1::PsmlData, d2::PsmlData) = (d1.npts == d2.npts) && all(d1._ .== d2._)
Base.:(≈)(d1::PsmlData, d2::PsmlData) = (d1.npts == d2.npts) && all(d1._ .≈ d2._)

"""
Radial function.
"""
struct PsmlRadfunc <: PsmlField
    grid::Union{PsmlGrid, Nothing}
    data::PsmlData
end

"""
Valence charge density.
"""
struct PsmlValenceCharge <: PsmlField
    total_charge::Float64
    is_unscreening_charge::Union{Bool, Nothing}
    rescaled_to_z_pseudo::Union{Bool, Nothing}
    annotation::Union{Dict{String, Any}, Nothing}
    radfunc::PsmlRadfunc
end

"""
Model core charge density.
"""
struct PsmlCoreCharge <: PsmlField
    matching_radius::Union{Float64, Nothing}
    number_of_continuous_derivatives::Int64
    annotation::Union{Dict{String, Any}, Nothing}
    radfunc::PsmlRadfunc
end

"""
Valence kinetic energy density.
"""
struct PsmlValenceKineticDensity <: PsmlField
    is_unscreening_tau::Union{Bool, Nothing}
    annotation::Union{Dict{String, Any}, Nothing}
    radfunc::PsmlRadfunc
end

"""
Model core kinetic energy density.
"""
struct PsmlCoreKineticDensity <: PsmlField
    matching_radius::Union{Float64, Nothing}
    number_of_continuous_derivatives::Int64
    annotation::Union{Dict{String, Any}, Nothing}
    radfunc::PsmlRadfunc
end

"""
Semi-local potential.
"""
struct PsmlPotential <: PsmlField
    flavor::Union{String, Nothing}
    l::String
    j::Union{Float64, Nothing}
    n::Int64
    rc::Float64
    eref::Union{Float64, Nothing}
    radfunc::PsmlRadfunc
end

"""
Semi-local potentials container.
"""
struct PsmlSemilocalPotentials <: PsmlField
    set::String
    flavor::Union{String, Nothing}
    annotation::Union{Dict{String, Any}, Nothing}
    grid::Union{PsmlGrid, Nothing}
    slps::Vector{PsmlPotential}
end

"""
Local charge density container.
"""
struct PsmlLocalCharge <: PsmlField
    radfunc::PsmlRadfunc
end

"""
Local potential.
"""
struct PsmlLocalPotential <: PsmlField
    type::String
    annotation::Union{Dict{String, Any}, Nothing}
    grid::Union{PsmlGrid, Nothing}
    radfunc::PsmlRadfunc
    local_charge::Union{PsmlLocalCharge, Nothing}
end

"""
Nonlocal projector.
"""
struct PsmlProjector <: PsmlField
    ekb::Float64
    eref::Union{Float64, Nothing}
    l::String
    j::Union{Float64, Nothing}
    seq::Int64
    type::String
    radfunc::PsmlRadfunc
end

"""
Nonlocal projectors / potential container.
"""
struct PsmlNonlocalProjectors <: PsmlField
    set::String
    annotation::Union{Dict{String, Any}, Nothing}
    grid::Union{PsmlGrid, Nothing}
    proj::Vector{PsmlProjector}
end

"""
Pseudo-atomic wavefunction.
"""
struct PsmlPseudoWf <: PsmlField
    l::String
    j::Union{Float64, Nothing}
    n::Int64
    energy_level::Union{Float64, Nothing}
    radfunc::PsmlRadfunc
end

"""
Pseudo-atomic wavefunctions container.
"""
struct PsmlPseudoWaveFunctions <: PsmlField
    set::String
    annotation::Union{Dict{String, Any}, Nothing}
    grid::Union{PsmlGrid, Nothing}
    pswf::Vector{PsmlPseudoWf}
end

"""
PSeudopotential Markup Language file contents. Information on file format and the meaning
of the quantities within the file can be found at the
[libPSML website](https://siesta-project.github.io/psml-docs/page/index.html)
or in the [PSML format paper](https://doi.org/10.1016/j.cpc.2018.02.011).
"""
struct PsmlFile <: PsPFile
    "PSML format version"
    version::String
    "Energy units"
    energy_unit::String
    "Length units"
    length_unit::String
    "Pseudopotential UUID"
    uuid::String
    "One or more provenance elements"
    provenance::Vector{PsmlProvenance}
    "Pseudo-atom specification"
    pseudo_atom_spec::PsmlPseudoAtomSpec
    "Optional top-level grid"
    grid::Union{PsmlGrid, Nothing}
    "Valence charge density"
    valence_charge::PsmlValenceCharge
    "Optional model core charge density"
    pseudocore_charge::Union{PsmlCoreCharge, Nothing}
    "Optional valence-kinetic energy density for MGGA"
    valence_kinetic_energy_density::Union{PsmlValenceKineticDensity, Nothing}
    "Optional model core-kinetic energy density for MGGA"
    pseudocore_kinetic_energy_density::Union{PsmlCoreKineticDensity, Nothing}
    "Optionally, one or more semi-local potential groups"
    semilocal_potentials::Union{Vector{PsmlSemilocalPotentials}, Nothing}
    "Local potential (optional if at least one semilocal potential group is present)"
    local_potential::Union{PsmlLocalPotential, Nothing}
    "Zero or more groups of non-local projectors"
    nonlocal_projectors::Union{Vector{PsmlNonlocalProjectors}, Nothing}
    "Optional groups of pseudo-atomic wavefunctions"
    pseudo_wave_functions::Union{Vector{PsmlPseudoWaveFunctions}, Nothing}
end
function PsmlFile(path::String)
    open(path, "r") do io
        PsmlFile(io)
    end
end

function PsmlFile(io::IO)
    Serde.deser_xml(PsmlFile, read(io, String))
end

function save_psp(io, file::PsmlFile, args...; kwargs...)
    write(io, Serde.to_xml(file; key="psml"))
end

# Ensure Vector values in deserialization for Vector-type elements even if only
# a single occurrence exists in the file.
function Serde.deser(::Type{PsmlFile}, ::Type{Vector{T}}, data::Dict{String, Any}) where {T}
    return [Serde.deser(T, data)]
end

# Generate custom_name (deserialization) and ser_name (serialization) methods so that
# the PSML convention of using "-" as a separator in the XML tags is converted to an
# internal convention of using "_".
const dash_underscore_exceptions = ("_", "energy_unit", "length_unit", "energy_level")
for type in [subtypes(PsmlField)..., PsmlFile]
    for fn in fieldnames(type)
        fn_str = string(fn)
        if occursin("_", fn_str) && !(fn_str in dash_underscore_exceptions)
            @eval function Serde.custom_name(
                ::Type{$type},
                ::Val{$Symbol($fn_str)}
            )
                join(split($fn_str, "_"), "-")
            end
            @eval function Serde.ser_name(
                ::Type{$type},
                ::Val{$Symbol($fn_str)}
            )
                Symbol(join(split($fn_str, "_"), "-"))
            end
        end
    end
end

identifier(file::PsmlFile) = file.uuid
format(::PsmlFile) = "PSML"
element(file::PsmlFile) = file.pseudo_atom_spec.atomic_label
function functional(file::PsmlFile)
    functionals = file.pseudo_atom_spec.exchange_correlation.libxc_info.functional
    return [Symbol(LIBXC_FUNCTIONALS_BY_ID[func.id]) for func in functionals]
end
valence_charge(file::PsmlFile) = file.pseudo_atom_spec.valence_configuration.total_valence_charge
is_norm_conserving(::PsmlFile) = true
is_ultrasoft(::PsmlFile) = false
is_paw(::PsmlFile) = false
has_spin_orbit(file::PsmlFile) = file.pseudo_atom_spec.relativity == "dirac"
has_model_core_charge_density(file::PsmlFile) = file.pseudo_atom_spec.core_corrections
