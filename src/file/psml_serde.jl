
struct InputFile
    name::String
    _::String
end

struct Provenance
    creator::String
    date::String
    annoatation::Dict{String,Any}
end

struct PseudoAtomSpec
    atomic_label::String
    atomic_number::Int
    z_pseudo::Int
    flavor::String
    relativity::String
    spin_dft::Bool
    core_correction::Bool
end

struct PSML
    psml::PSML
    provenance::Provenance
    pseudo_atom_spec::PseudoAtomSpec
    grid::Grid
    valence_charge::ValenceCharge
    pseudocore_charge::PseudocoreCharge
    semilocal_potentials::Vector{SemilocalPotential}
    local_potential::LocalPotential
    nonlocal_projectors::Vector{NonlocalProjector}
end
