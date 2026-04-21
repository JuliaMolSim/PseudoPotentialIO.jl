module PseudoPotentialIO
using EzXML
using Printf
using PeriodicTable
using Statistics
using Serde
using PrecompileTools: @setup_workload, @compile_workload

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
include("file/upf2_serde.jl")
include("file/psml_serde.jl")

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

@setup_workload begin
    text = """
    <UPF version="2.0.1">
        <PP_INFO>text<PP_INPUTFILE>text</PP_INPUTFILE></PP_INFO>
        <PP_HEADER
            generated="" author="" date="" comment="" element="H"
            pseudo_type="PAW" relativistic="fully-relativistic"
            is_ultrasoft="t" is_paw="t" is_coulomb="f" has_so="f"
            has_wfc="t" has_gipaw="t"
            core_correction="t" functional="" z_valence="1.0"
            total_psenergy="1.0" rho_cutoff="1.0" l_max="1" l_local="-1"
            mesh_size="1" number_of_wfc="1" number_of_proj="1"
        />
        <PP_MESH>
            <PP_R>0.0</PP_R>
            <PP_RAB>0.0</PP_RAB>
        </PP_MESH>
        <PP_LOCAL>0.0</PP_LOCAL>
        <PP_NONLOCAL>
            <PP_BETA.1>0.0</PP_BETA.1>
            <PP_DIJ>1.0</PP_DIJ>
            <PP_AUGMENTATION q_with_l="f" nqf="1" nqlc="1">
                <PP_Q>1.0</PP_Q>
                <PP_QIJ.1.1 first_index="1" second_index="1" composite_index="1">
                    1.0
                </PP_QIJ.1.1>
            </PP_AUGMENTATION>
        </PP_NONLOCAL>
        <PP_PSWFC>
            <PP_CHIS.1>0.0</PP_CHIS.1>
        </PP_PSWFC>
        <PP_FULL_WFC number_of_wfc="1">
            <PP_AEWFC.1>0.0</PP_AEWFC.1>
            <PP_PSWFC.1>0.0</PP_PSWFC.1>
        </PP_FULL_WFC>
        <PP_PAW paw_data_format="2">
            <PP_OCCUPATIONS>0.0</PP_OCCUPATIONS>
            <PP_AE_NLCC>0.0</PP_AE_NLCC>
            <PP_AE_VLOC>0.0</PP_AE_VLOC>
        </PP_PAW>
        <PP_SPINORB>
            <PP_RELBETA.1 index="1" lll="0" jjj="0.5"/>
            <PP_RELWFC.1 index="1" lchi="0" jchi="0.5" nn="1"/>
        </PP_SPINORB>
        <PP_GIPAW gipaw_data_format="2">
            <PP_GIPAW_CORE_ORBITALS>
                <PP_GIPAW_CORE_ORBITAL.1>0.0</PP_GIPAW_CORE_ORBITAL.1>
            </PP_GIPAW_CORE_ORBITALS>
            <PP_GIPAW_ORBITALS>
                <PP_GIPAW_ORBITAL.1>
                    <PP_GIPAW_WFS_AE>0.0</PP_GIPAW_WFS_AE>
                    <PP_GIPAW_WFS_PS>0.0</PP_GIPAW_WFS_PS>
                </PP_GIPAW_ORBITAL.1>
            </PP_GIPAW_ORBITALS>
            <PP_GIPAW_VLOCAL>
                <PP_GIPAW_VLOCAL_AE>0.0</PP_GIPAW_VLOCAL_AE>
                <PP_GIPAW_VLOCAL_PS>0.0</PP_GIPAW_VLOCAL_PS>
            </PP_GIPAW_VLOCAL>
        </PP_GIPAW>
    </UPF>
    """
    @compile_workload begin
        upf = deser_xml(UPFSerde, text);
    end
end

end
