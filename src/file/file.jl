@doc raw"""
Abstract type representing a pseudopotential file.

The structure of the data should closely mirror the format of the file, and the values
of quantities should be exactly those found in the file.
"""
abstract type PsPFile end

"""
Identifying data.
"""
function identifier(file::PsPFile)::String end

"""
Pseudopotential file format.
"""
function format(file::PsPFile)::String end

"""
The element which the pseudopotential was constructed to reproduce.
"""
function element(file::PsPFile)::PeriodicTable.Element end

"""
Exchange-correlation functional in LibXC format.
"""
function functional(file::PsPFile)::Vector{Functional} end

"""
Pseudo-atomic core charge.
"""
function valence_charge(file::PsPFile)::Int end

"""
Whether the pseudopotential is of the norm-conserving kind.
"""
function is_norm_conserving(file::PsPFile)::Bool end

"""
Whether the pseudopotential is of the ultrasoft kind.
"""
function is_ultrasoft(file::PsPFile)::Bool end

"""
Whether the pseudopotential is of the plane-augmented wave kind.
"""
function is_paw(file::PsPFile)::Bool end

"""
Whether the pseudopotential contains relativistic spin-orbit coupling data.
"""
function has_spin_orbit(file::PsPFile)::Bool end

"""
Whether the pseudopotential contains a model core charge density used in non-linear core correction.
"""
function has_model_core_charge_density(file::PsPFile)::Bool end

"""
Formalism of the pseudopotential.
"""
function formalism(file::PsPFile)::String
    # Must be done in this order because some UPF PAW files also have is_ultrasoft == true
    is_paw(file) && return "projector-augmented wave"
    is_ultrasoft(file) && return "ultrasoft"
    is_norm_conserving(file) && return "norm-conserving"
end

Base.Broadcast.broadcastable(file::PsPFile) = Ref(file)

function Base.show(io::IO, file::PsPFile)
    typename = string(typeof(file))
    el = element(file)
    z = valence_charge(file)
    functional = map(f -> String(f.identifier), libxc_functional(file))
    formalism = formalism(file)
    nlcc = has_model_core_charge_density(file)
    return print(io, "$typename(element=$el, z_valence=$z, xc=$functional, formalism=$formalism, nlcc=$nlcc)")
end

function Base.show(io::IO, ::MIME"text/plain", file::PsPFile)
    println(io, typeof(file))
    @printf "%032s: %s\n" "identifier" identifier(file)
    @printf "%032s: %s\n" "format" format(file)
    @printf "%032s: %s\n" "element" element(file)
    @printf "%032s: %s\n" "exchange-correlation (Libxc)" join(map(f -> String(f.identifier), functional(file)), '+')
    @printf "%032s: %d\n" "valence charge" valence_charge(file)
    @printf "%032s: %s\n" "spin-orbit coupling" has_spin_orbit(file)
    @printf "%032s: %s\n" "model core charge density" has_model_core_charge_density(file)
    @printf "%032s: %s\n" "formalism" formalism(file)
end
