function load_upf(filepath::AbstractString)
    open(filepath, "r") do io
        return load_upf(io)
    end
end

function load_upf(io::IO)
    @warn "load_upf is deprecated, use load_psp_file and the new PspFile structures"
    format_version = _get_upf_version(io)
    new_upf = UpfFile(io)

    old_upf = Dict()

    q_with_l = isnothing(new_upf.nonlocal.augmentation) ? false :
               new_upf.nonlocal.augmentation.q_with_l
    old_upf["header"] = Dict("format_version" => format_version,
                             "has_wfc" => new_upf.header.has_wfc,
                             "number_of_wfc" => new_upf.header.number_of_wfc,
                             "number_of_proj" => new_upf.header.number_of_proj,
                             "ecutwfc" => new_upf.header.wfc_cutoff,
                             "has_so" => new_upf.header.has_so,
                             "q_with_l" => q_with_l,
                             "has_gipaw" => new_upf.header.has_gipaw,
                             "core_correction" => new_upf.header.core_correction,
                             "element" => new_upf.header.element,
                             "is_ultrasoft" => new_upf.header.is_ultrasoft,
                             "version" => 0,
                             "is_paw" => new_upf.header.is_paw,
                             "is_coulomb" => new_upf.header.is_coulomb,
                             "l_max" => new_upf.header.l_max,
                             "z_valence" => new_upf.header.z_valence,
                             "total_psenergy" => new_upf.header.total_psenergy,
                             "pseudo_type" => new_upf.header.pseudo_type,
                             "mesh_size" => new_upf.header.mesh_size,
                             "ecutrho" => new_upf.header.rho_cutoff)
    if old_upf["header"]["is_paw"]
        @warn "PAW in UPF v$(format_version) is not implemented."
    elseif old_upf["header"]["is_ultrasoft"]
        @warn "Ultrasoft in UPF v$(format_version) is not implemented."
    end

    old_upf["radial_grid"] = new_upf.mesh.r
    old_upf["radial_grid_derivative"] = new_upf.mesh.rab
    (mesh_type, mesh_a, mesh_b) = guess_mesh_type(old_upf["radial_grid"],
                                                  old_upf["radial_grid_derivative"])
    if mesh_type == "unknown"
        raise(ExceptionError("Unknown mesh type"))
    end
    old_upf["radial_grid_parameters"] = Dict("dx" => 0.0,
                                             "mesh" => length(old_upf["radial_grid"]),
                                             "xmin" => 0.0,
                                             "rmax" => maximum(old_upf["radial_grid"]),
                                             "zmesh" => 0.0,
                                             "a" => mesh_a,
                                             "b" => mesh_b,
                                             "mesh_type" => mesh_type)

    if !isnothing(new_upf.nlcc)
        old_upf["core_charge_density"] = new_upf.nlcc
    end
    old_upf["local_potential"] = new_upf.local_

    old_upf["beta_projectors"] = [Dict("label" => "",
                                       "index" => beta.index,
                                       "angular_momentum" => beta.angular_momentum,
                                       "cutoff_radius_index" => beta.cutoff_radius_index,
                                       "radial_function" => beta.beta[1:beta.cutoff_radius_index],
                                       "cutoff_radius" => 0.0)
                                  for beta in new_upf.nonlocal.betas]
    old_upf["D_ion"] = new_upf.nonlocal.dij

    if isnothing(new_upf.pswfc)
        old_upf["atomic_wave_functions"] = []
    else
        old_upf["atomic_wave_functions"] = [Dict("label" => chi.label,
                                                 "angular_momentum" => chi.l,
                                                 "occupation" => chi.occupation,
                                                 "radial_function" => chi.chi)
                                            for chi in new_upf.pswfc]
    end

    old_upf["total_charge_density"] = new_upf.rhoatom

    if !isnothing(new_upf.spin_orb)
        for relwfc in new_upf.spin_orb.relwfcs
            old_upf["atomic_wave_functions"][relwfc.index]["label"] = relwfc.els
            old_upf["atomic_wave_functions"][relwfc.index]["principal_quantum_number"] = relwfc.nn
            old_upf["atomic_wave_functions"][relwfc.index]["total_angular_momentum"] = relwfc.jchi
            old_upf["atomic_wave_functions"][relwfc.index]["angular_momentum"] = relwfc.lchi
            old_upf["atomic_wave_functions"][relwfc.index]["occupation"] = isnothing(relwfc.oc) ?
                                                                           0.0 : relwfc.oc
        end
        for relbeta in new_upf.spin_orb.relbetas
            old_upf["beta_projectors"][relbeta.index]["angular_momentum"] = relbeta.lll
            old_upf["beta_projectors"][relbeta.index]["total_angular_momentum"] = relbeta.jjj
        end
    end

    return old_upf
end
