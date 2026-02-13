module PseudoPotentialIO
using EzXML
using Printf
using PeriodicTable
using Statistics

using PeriodicTable: PeriodicTable
import Base.Broadcast.broadcastable

## DocStringExtensions Templates
using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(TYPEDSIGNATURES)
                                         $(DOCSTRING)
                                         $(METHODLIST)
                                         """

@template TYPES = """
                  $(TYPEDEF)
                  $(DOCSTRING)
                  $(TYPEDFIELDS)
                  """
## Data
include("data/libxc_functionals.jl")
include("data/psp8_functionals.jl")
include("data/upf_functionals.jl")

## Common utilities
include("common/mesh.jl")
include("common/xml.jl")

## File datastructures and interface
# These are only 'public' to avoid clashes in the ecosystem with functions
# of similar names, e.g. in Libxc.jl and AtomsBase.jl
export PsPFile
#= TODO: Can be changed to this once 1.11 is hard enforced
public is_norm_conserving
public is_ultrasoft
public is_paw
public has_spin_orbit
public has_model_core_charge_density
public identifier
public format
public element
public functional
public valence_charge
public formalism
=#
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("public identifier, format, element, functional, valence_charge, formalism, is_norm_conserving, is_ultrasoft, is_paw, has_spin_orbit, has_model_core_charge_density"))
end
include("file/file.jl")

export UpfFile
include("file/upf.jl")
include("file/upf1.jl")
include("file/upf2.jl")

export Psp8File
include("file/psp8.jl")

export HghFile
include("file/hgh.jl")

## Loading/listing functions
export load_psp_file
include("load.jl")

## Save to file
export save_psp_file
include("save.jl")

# Conversion
include("conversion/to_upf.jl")

end
