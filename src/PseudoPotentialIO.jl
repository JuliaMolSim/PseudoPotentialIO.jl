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

## File datastructures and interface
export PsPFile
export identifier
export format
export element
export functional
export valence_charge
export is_norm_conserving
export is_ultrasoft
export is_paw
export has_spin_orbit
export has_model_core_charge_density
export formalism
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
