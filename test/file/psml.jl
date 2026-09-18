@testset "PSML" begin
    @testset "Internal data consistency" begin
        # Test after a round-trip
        @testset "$(case)" for (case, path) in PSML_CASE_FILEPATHS
            file1 = load_psp_file(path)
            file2 = mktempdir() do tmp
                save_psp_file(joinpath(tmp, "tmp.psml"), file1)
                load_psp_file(joinpath(tmp, "tmp.psml"))
            end
            # Strict equality between the input/output.
            # We can do this because essentially all PSML files have numerical arrays
            # written as %19.12E, including ours. If this is ever not true,
            # strict equality might fail, and we may need to loosen this to approximation.
            @test file1.version == file2.version
            @test file1.energy_unit == file2.energy_unit
            @test file1.length_unit == file2.length_unit
            @test file1.uuid == file2.uuid
            @test length(file1.provenance) == length(file2.provenance)
            @test all(file1.provenance .== file2.provenance)
            spec1 = file1.pseudo_atom_spec
            spec2 = file2.pseudo_atom_spec
            @test spec1.atomic_label == spec2.atomic_label
            @test spec1.atomic_number == spec2.atomic_number
            @test spec1.z_pseudo == spec2.z_pseudo
            @test spec1.core_corrections == spec2.core_corrections
            @test spec1.meta_gga == spec2.meta_gga
            @test spec1.relativity == spec2.relativity
            @test spec1.spin_dft == spec2.spin_dft
            @test spec1.flavor == spec2.flavor
            @test isnothing(spec1.annotation) == isnothing(spec2.annotation)
            if !isnothing(spec1.annotation)
                for (k1, v1) in spec1.annotation
                    @test haskey(spec2.annotation, k1)
                    @test v1 == spec2.annotation[k1]
                end
            end
            @test spec1.exchange_correlation == spec2.exchange_correlation
            @test spec1.valence_configuration == spec2.valence_configuration
            @test spec1.core_configuration == spec2.core_configuration
            @test file1.grid == file2.grid
            @test file1.valence_charge == file2.valence_charge
            @test file1.pseudocore_charge == file2.pseudocore_charge
            @test file1.valence_kinetic_energy_density == file2.valence_kinetic_energy_density
            @test file1.pseudocore_kinetic_energy_density == file1.pseudocore_kinetic_energy_density
            @test isnothing(file1.semilocal_potentials) == isnothing(file2.semilocal_potentials)
            if !isnothing(file1.semilocal_potentials)
                @test all(file1.semilocal_potentials .== file2.semilocal_potentials)
            end
            @test file1.local_potential == file1.local_potential
            @test isnothing(file1.nonlocal_projectors) == isnothing(file2.nonlocal_projectors)
            if !isnothing(file1.nonlocal_projectors)
                @test all(file1.nonlocal_projectors .== file1.nonlocal_projectors)
            end
            @test isnothing(file1.pseudo_wave_functions) == isnothing(file2.pseudo_wave_functions)
            if !isnothing(file1.pseudo_wave_functions)
                @test all(file1.pseudo_wave_functions .== file1.pseudo_wave_functions)
            end
        end
    end
    @testset "Contents" begin
        @testset "14_Si_UPF_r.psml" begin
            filename = "14_Si_UPF_r.psml"
            file = load_psp_file(PSML_CASE_FILEPATHS[filename])

            @test file.version == "1.1"
            @test file.energy_unit == "hartree"
            @test file.length_unit == "bohr"
            @test file.uuid == "a9c2bc50-7302-11e7-4671-34606fcfdcae"

            @test length(file.provenance) == 1
            provenance = first(file.provenance)
            @test provenance.creator == "ONCVPSP-3.3.0+psml-3.3.0-73 (fully-relativistic)"
            @test provenance.date == "2017-07-27"
            @test provenance.annotation["action"] == "semilocal-pseudopotential-generation"
            @test provenance.annotation["action-cont"] == "projectors-generation"
            @test length(provenance.input_file) == 1

            @test provenance.input_file[1].name == "oncvpsp-input"

            spec = file.pseudo_atom_spec
            @test spec.atomic_label == "Si"
            @test spec.atomic_number == 14
            @test spec.z_pseudo == 4.0
            @test spec.flavor == "Hamann's oncvpsp"
            @test spec.relativity == "dirac"
            @test !spec.spin_dft
            @test spec.core_corrections
            @test spec.annotation["pseudo-energy"] == "-3.79205493419"

            exc = spec.exchange_correlation
            @test exc.annotation["oncvpsp-xc-code"] == "3"
            @test exc.annotation["oncvpsp-xc-type"] == "LDA -- Ceperley-Alder Perdew-Zunger"
            @test exc.libxc_info.number_of_functionals == 2
            @test length(exc.libxc_info.functional) == 2
            @test exc.libxc_info.functional[1].name == "Slater exchange (LDA)"
            @test exc.libxc_info.functional[1].type == "exchange"
            @test exc.libxc_info.functional[1].id == 1
            @test exc.libxc_info.functional[2].name == "Perdew & Zunger (LDA)"
            @test exc.libxc_info.functional[2].type == "correlation"
            @test exc.libxc_info.functional[2].id == 9

            val = spec.valence_configuration
            @test val.total_valence_charge == 4
            @test length(val.shell) == 2
            @test val.shell[1].n == 3
            @test val.shell[1].l == "s"
            @test val.shell[1].occupation == 2.0
            @test val.shell[2].n == 3
            @test val.shell[2].l == "p"
            @test val.shell[2].occupation == 2.0

            grid = file.grid
            @test grid.npts == 457
            @test grid.annotation["type"] == "sampled from oncvpsp log grid"
            @test grid.annotation["recipe"] == "r(i:1..N) = r1*exp(a*(i-1)) + r=0; resampled"
            @test grid.annotation["recipe-cont"] == "r1: scale; a: step"
            @test grid.annotation["scale"] == "0.357142857143E-04"
            @test grid.annotation["step"] == "0.119285708653E-01"
            @test grid.annotation["delta"] == "0.500000000000E-02"
            @test grid.annotation["rmax"] == "44.1697227368"
            @test length(grid.grid_data._) == 457
            @test grid.grid_data._[1] == 0.0
            @test grid.grid_data._[2] == 5.043762840785E-03
            @test grid.grid_data._[5] == 2.012306398724E-02
            @test grid.grid_data._[9] == 4.067657200651E-02
            @test grid.grid_data._[end] == 4.416972273678E+01

            rhoval = file.valence_charge
            @test rhoval.total_charge == 4.0
            @test rhoval.is_unscreening_charge
            @test !rhoval.rescaled_to_z_pseudo
            @test rhoval.radfunc.data.npts == 457
            @test length(rhoval.radfunc.data._) == 457
            @test rhoval.radfunc.data._[1] == 7.349279715428E-02
            @test rhoval.radfunc.data._[2] == 7.350892183566E-02
            @test rhoval.radfunc.data._[5] == 7.374948026829E-02
            @test rhoval.radfunc.data._[9] == 7.454184602293E-02
            @test rhoval.radfunc.data._[end] == 0.000000000000E+00

            rhocore = file.pseudocore_charge
            @test !isnothing(rhocore)
            @test rhocore.matching_radius == 1.40590689067
            @test rhocore.number_of_continuous_derivatives == 4
            @test rhocore.annotation["model-charge-form"] == "Polynomial"
            @test rhocore.radfunc.data.npts == 275
            @test length(rhocore.radfunc.data._) == 275
            @test rhocore.radfunc.data._[1] == 7.024804422741E+00
            @test rhocore.radfunc.data._[2] == 7.024804399941E+00
            @test rhocore.radfunc.data._[5] == 7.024798209454E+00
            @test rhocore.radfunc.data._[9] == 7.024702449567E+00
            @test rhocore.radfunc.data._[end] == 9.586604490682E-11

            @test length(file.semilocal_potentials) == 1
            vsl = first(file.semilocal_potentials)
            @test vsl.set == "lj"
            @test length(vsl.slps) == 5
            @test vsl.slps[1].n == 3
            @test vsl.slps[1].l == "s"
            @test vsl.slps[1].j == 0.5
            @test vsl.slps[1].rc == 2.70947249846
            @test vsl.slps[1].eref == -0.399979628917
            @test vsl.slps[1].radfunc.data.npts == 275
            @test length(vsl.slps[1].radfunc.data._) == 275
            @test vsl.slps[1].radfunc.data._[1] == 5.658531131228E+00
            @test vsl.slps[1].radfunc.data._[2] == 5.657736450274E+00
            @test vsl.slps[1].radfunc.data._[5] == 5.645895916023E+00
            @test vsl.slps[1].radfunc.data._[9] == 5.607093899264E+00
            @test vsl.slps[1].radfunc.data._[end] == -7.939420485318E-01

            vloc = file.local_potential
            @test !isnothing(vloc)
            @test vloc.type == "oncv-fit"
            @test vloc.radfunc.data.npts == 275
            @test length(vloc.radfunc.data._) == 275
            @test vloc.radfunc.data._[1] == -2.739597733669E+00
            @test vloc.radfunc.data._[2] == -2.739587948184E+00
            @test vloc.radfunc.data._[5] == -2.739442186089E+00
            @test vloc.radfunc.data._[9] == -2.738965169301E+00
            @test vloc.radfunc.data._[end] == -7.939420485318E-01

            @test length(file.nonlocal_projectors) == 1
            vnl = first(file.nonlocal_projectors)
            @test vnl.set == "lj"
            @test length(vnl.proj) == 10
            @test vnl.proj[1].l == "s"
            @test vnl.proj[1].j == 0.5
            @test vnl.proj[1].seq == 1
            @test vnl.proj[1].ekb == -0.119802984280
            @test vnl.proj[1].eref == -0.399979628917
            @test vnl.proj[1].type == "oncv"
            @test vnl.proj[1].radfunc.data.npts == 275
            @test length(vnl.proj[1].radfunc.data._) == 275

            @test length(file.pseudo_wave_functions) == 1
            wfn = first(file.pseudo_wave_functions)
            @test wfn.set == "lj"
            @test length(wfn.pswf) == 3
            @test wfn.pswf[1].n == 3
            @test wfn.pswf[1].l == "s"
            @test wfn.pswf[1].j == 0.5
            @test wfn.pswf[1].energy_level == -0.399979628917
        end
        @testset "56_Ba_m.psml" begin
            filename = "56_Ba_m.psml"
            file = load_psp_file(PSML_CASE_FILEPATHS[filename])

            @test file.version == "1.2"
            @test file.energy_unit == "hartree"
            @test file.length_unit == "bohr"
            @test file.uuid == "c1d00520-3d8b-11f1-7c86-6936d773923a"

            @test length(file.provenance) == 1
            provenance = first(file.provenance)
            @test provenance.creator == "METAPSP-1.0.2+psml-4.0.1-81 (non-relativistic)"
            @test provenance.date == "2026-04-21"
            @test provenance.annotation["action-cont"] == "projectors-generation"
            @test length(provenance.input_file) == 1

            @test provenance.input_file[1].name == "oncvpsp-input"

            spec = file.pseudo_atom_spec
            @test spec.atomic_label == "Ba"
            @test spec.atomic_number == 56
            @test spec.z_pseudo == 10.0
            @test spec.flavor == "Hamann`s oncvpsp"
            @test spec.relativity == "no"
            @test !spec.spin_dft
            @test spec.core_corrections
            @test spec.meta_gga
            @test spec.annotation["pseudo-energy"] == "-14.0606365385"
            @test spec.annotation["cutoff_hint_low"] == "19"
            @test spec.annotation["cutoff_hint_normal"] == "24"
            @test spec.annotation["cutoff_hint_high"] == "28"

            exc = spec.exchange_correlation
            @test exc.annotation["oncvpsp-xc-code"] == "5"
            @test exc.annotation["oncvpsp-xc-type"] == "METAGGA -- Perdew"
            @test exc.libxc_info.number_of_functionals == 2
            @test length(exc.libxc_info.functional) == 2
            @test exc.libxc_info.functional[1].name == "R2SCAN01 (MGGA)"
            @test exc.libxc_info.functional[1].type == "exchange"
            @test exc.libxc_info.functional[1].id == 645
            @test exc.libxc_info.functional[2].name == "R2SCAN01 (MGGA)"
            @test exc.libxc_info.functional[2].type == "correlation"
            @test exc.libxc_info.functional[2].id == 642

            val = spec.valence_configuration
            @test val.total_valence_charge == 10
            @test length(val.shell) == 3
            @test val.shell[1].n == 5
            @test val.shell[1].l == "s"
            @test val.shell[1].occupation == 2.0
            @test val.shell[2].n == 5
            @test val.shell[2].l == "p"
            @test val.shell[2].occupation == 6.0
            @test val.shell[3].n == 6
            @test val.shell[3].l == "s"
            @test val.shell[3].occupation == 2.0

            grid = file.grid
            @test grid.npts == 1363
            @test grid.annotation["type"] == "sampled from oncvpsp log grid"
            @test grid.annotation["recipe"] == "r(i:1..N) = r1*exp(a*(i-1)) + r=0; resampled"
            @test grid.annotation["recipe-cont"] == "r1: scale; a: step"
            @test grid.annotation["scale"] == "0.892857142857E-05"
            @test grid.annotation["step"] == "0.299550897980E-02"
            @test grid.annotation["delta"] == "0.500000000000E-02"
            @test grid.annotation["rmax"] == "44.8629642481"
            @test length(grid.grid_data._) == 1363
            @test grid.grid_data._[1] == 0.0
            @test grid.grid_data._[2] == 5.007874647428E-03
            @test grid.grid_data._[5] == 2.004404817718E-02
            @test grid.grid_data._[9] == 4.028124004324E-02
            @test grid.grid_data._[end] == 4.486296424810E+01

            rhoval = file.valence_charge
            @test rhoval.total_charge == 10.0
            @test rhoval.is_unscreening_charge
            @test !rhoval.rescaled_to_z_pseudo
            @test rhoval.radfunc.data.npts == 1363
            @test length(rhoval.radfunc.data._) == 1363
            @test rhoval.radfunc.data._[1] == 1.274602136897E-01
            @test rhoval.radfunc.data._[2] == 1.275262822970E-01
            @test rhoval.radfunc.data._[5] == 1.285187002847E-01
            @test rhoval.radfunc.data._[9] == 1.317358880843E-01
            @test rhoval.radfunc.data._[end] == 0.000000000000E+00

            tauval = file.valence_kinetic_energy_density
            @test !isnothing(tauval)
            @test tauval.is_unscreening_tau
            @test tauval.radfunc.data.npts == 1363
            @test length(tauval.radfunc.data._) == 1363
            @test tauval.radfunc.data._[1] == 2.818357349572E+00
            @test tauval.radfunc.data._[2] == 2.818401001816E+00
            @test tauval.radfunc.data._[5] == 2.819055534715E+00
            @test tauval.radfunc.data._[9] == 2.821162451665E+00
            @test tauval.radfunc.data._[end] == 0.000000000000E+00

            rhocore = file.pseudocore_charge
            @test !isnothing(rhocore)
            @test rhocore.matching_radius == 2.70154864719
            @test rhocore.number_of_continuous_derivatives == 4
            @test rhocore.annotation["model-charge-form"] == "MGGA function adjusted by input data fcfact and rcfact"
            @test rhocore.radfunc.data.npts == 692
            @test length(rhocore.radfunc.data._) == 692
            @test rhocore.radfunc.data._[1] == 1.005787209781E+01
            @test rhocore.radfunc.data._[2] == 1.005768106996E+01
            @test rhocore.radfunc.data._[5] == 1.005481212583E+01
            @test rhocore.radfunc.data._[9] == 1.004551796741E+01
            @test rhocore.radfunc.data._[end] == 6.141529677639E-12

            taucore = file.pseudocore_kinetic_energy_density
            @test !isnothing(taucore)
            @test taucore.matching_radius == 2.70154864719
            @test taucore.number_of_continuous_derivatives == 4
            @test taucore.annotation["model-tau-form"] == "MGGA function adjusted by input data fcfact and rcfact"
            @test taucore.radfunc.data.npts == 692
            @test length(taucore.radfunc.data._) == 692
            @test taucore.radfunc.data._[1] == 1.179060763006E+01
            @test taucore.radfunc.data._[2] == 1.179043750827E+01
            @test taucore.radfunc.data._[5] == 1.178788247121E+01
            @test taucore.radfunc.data._[9] == 1.177960433965E+01
            @test taucore.radfunc.data._[end] == 2.213427807804E-11

            @test isnothing(file.semilocal_potentials)

            vloc = file.local_potential
            @test vloc.type == "oncv-fit"
            @test vloc.radfunc.data.npts == 692
            @test length(vloc.radfunc.data._) == 692
            @test vloc.radfunc.data._[1] == -9.747863845348E+00
            @test vloc.radfunc.data._[2] == -9.747818019644E+00
            @test vloc.radfunc.data._[5] == -9.747146985483E+00
            @test vloc.radfunc.data._[9] == -9.744972378241E+00
            @test vloc.radfunc.data._[end] == -1.663558550343E+00

            @test length(file.nonlocal_projectors) == 1
            vnl = first(file.nonlocal_projectors)
            @test vnl.set == "non_relativistic"
            @test length(vnl.proj) == 9
            @test vnl.proj[1].l == "s"
            @test isnothing(vnl.proj[1].j)
            @test vnl.proj[1].seq == 1
            @test vnl.proj[1].ekb == 6.59722722790
            @test vnl.proj[1].eref == -1.19528760596
            @test vnl.proj[1].type == "oncv"
            @test vnl.proj[1].radfunc.data.npts == 692
            @test length(vnl.proj[1].radfunc.data._) == 692

            @test length(file.pseudo_wave_functions) == 1
            wfn = first(file.pseudo_wave_functions)
            @test wfn.set == "non_relativistic"
            @test length(wfn.pswf) == 3
            @test wfn.pswf[1].n == 5
            @test wfn.pswf[1].l == "s"
            @test isnothing(wfn.pswf[1].j)
            @test wfn.pswf[1].energy_level == -1.19528760596
        end
    end
end
