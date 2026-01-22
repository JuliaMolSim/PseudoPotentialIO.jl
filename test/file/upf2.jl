@testset "UPF v2.0.1" begin

    @testset "Internal data consistency" begin
        # Test after a circle, not data is lost
        for filepath in values(UPF2_CASE_FILEPATHS)
            psp = load_psp_file(filepath)
            
            # create a tmp file to save for test purpose
            tmpdir = mktempdir(; cleanup = true)
            outpath = joinpath(tmpdir, "tmp.upf")
            save_psp_file(outpath, psp, format(psp))

            newpsp = load_psp_file(outpath)

            # test no fields are lost (the fields that can be nothing)
            for n in fieldnames(UpfFile)
                if !isnothing(getfield(psp, n))
                    @test !isnothing(getfield(newpsp, n)) || "$(filepath): $(n) of $(outpath) file is nothing"
                end
            end

            @test psp.header == newpsp.header
            @test psp.info == newpsp.info
            @test psp.mesh.r ≈ newpsp.mesh.r
            @test psp.mesh.rab ≈ newpsp.mesh.rab
            @test psp.local_ ≈ newpsp.local_
            if length(psp.nonlocal.betas) > 0
                @test psp.nonlocal.betas[1].beta ≈ newpsp.nonlocal.betas[1].beta
            end
            @test psp.nonlocal.dij ≈ newpsp.nonlocal.dij
            @test psp.rhoatom ≈ newpsp.rhoatom

            # check io and recur_io are the same

            @test format(psp) == "UPF v2.0.1"

            # UPF v2.0.1 has a different augmentation data format from UPF v1.old
            if psp.header.is_ultrasoft | psp.header.is_paw
                augmentation = psp.nonlocal.augmentation

                @test isnothing(augmentation.rinner)

                if psp.header.is_paw
                    @test !isnothing(augmentation.multipoles)
                    #TODO @test length(augmentation.multipoles) ==
                end
            end
        end
    end


    # TODO: either qijs or qijls should be present
    # But UPF 2.0.1 parser not yet have qijs supported

    #@testset "UPF1 -> UPF2" begin
    #    for filepath in values(UPF1_CASE_FILEPATHS)
    #        psp = load_psp_file(filepath)

    #        @test format(psp) == "UPF v1.old"

    #        # create a tmp file to save for test purpose
    #        tmpdir = mktempdir(; cleanup = true)
    #        outpath = joinpath(tmpdir, "tmp.upf")
    #        save_psp_file(outpath, psp)

    #        println(outpath)

    #        newpsp = load_psp_file(outpath)

    #        @test format(newpsp) == "UPF v2.0.1"

    #        # test no fields are lost (the fields that can be nothing)
    #        for n in fieldnames(UpfFile)
    #            if !isnothing(getfield(psp, n))
    #                @test !isnothing(getfield(newpsp, n)) || "$(filepath): $(n) of $(outpath) file is nothing"
    #            end
    #        end

    #        @test psp.mesh.r ≈ newpsp.mesh.r
    #        @test psp.mesh.rab ≈ newpsp.mesh.rab
    #        @test psp.local_ ≈ newpsp.local_
    #        if length(psp.nonlocal.betas) > 0
    #            @test psp.nonlocal.betas[1].beta ≈ newpsp.nonlocal.betas[1].beta
    #        end
    #        @test psp.nonlocal.dij ≈ newpsp.nonlocal.dij
    #        @test psp.rhoatom ≈ newpsp.rhoatom
    #    end
    #end

    @testset "Mg.upf" begin
        filename = "Mg.upf"
        file = load_psp_file(UPF2_CASE_FILEPATHS[filename])

        header = file.header
        @test header.generated == "Generated using ONCVPSP code by D. R. Hamann"
        @test header.author == "anonymous"
        @test header.date == "180509"
        @test header.comment == ""
        @test header.element == "Mg"
        @test header.pseudo_type == "NC"
        @test header.relativistic == "full"
        @test !header.is_ultrasoft
        @test !header.is_paw
        @test !header.is_coulomb
        @test header.has_so
        @test !header.has_wfc
        @test !header.has_gipaw
        @test isnothing(header.paw_as_gipaw)
        @test !header.core_correction
        @test header.functional == "PBESOL"
        @test header.z_valence == 10.00
        @test header.total_psenergy == -1.18479378822E+02
        @test isnothing(header.wfc_cutoff)
        @test header.rho_cutoff == 1.50900000000E+01
        @test header.l_max == 1
        @test isnothing(header.l_max_rho)
        @test header.l_local == -1
        @test header.mesh_size == 1510
        @test header.number_of_wfc == 4
        @test header.number_of_proj == 6

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 0.0400
        @test mesh.rab[5] == 0.0100
        @test mesh.mesh == 1510
        @test isnothing(mesh.dx)
        @test isnothing(mesh.xmin)
        @test isnothing(mesh.rmax)
        @test isnothing(mesh.zmesh)

        @test isnothing(file.nlcc)

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -3.9295603163E+01

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj
        for (i, beta) in enumerate(betas)
            @test beta.index == i
            @test (length(beta.beta) == beta.cutoff_radius_index) | (length(beta.beta) == file.header.mesh_size)
            @test isnothing(beta.norm_conserving_radius)
            @test isnothing(beta.ultrasoft_cutoff_radius)
            @test isnothing(beta.label)
        end

        @test betas[1].angular_momentum == 0
        @test betas[1].cutoff_radius_index == 160
        @test betas[1].cutoff_radius == 1.5900000000E+00
        @test betas[1].beta[5] == 1.2872222840E-01

        @test betas[2].angular_momentum == 0
        @test betas[2].cutoff_radius_index == 160
        @test betas[2].cutoff_radius == 1.5900000000E+00
        @test betas[2].beta[5] == 4.8344886498E-01

        @test betas[3].angular_momentum == 1
        @test betas[3].cutoff_radius_index == 160
        @test betas[3].cutoff_radius == 1.5900000000E+00
        @test betas[3].beta[5] == 3.1550982941E-02

        @test betas[4].angular_momentum == 1
        @test betas[4].cutoff_radius_index == 160
        @test betas[4].cutoff_radius == 1.5900000000E+00
        @test betas[4].beta[5] == 3.1541055365E-02

        @test betas[5].angular_momentum == 1
        @test betas[5].cutoff_radius_index == 160
        @test betas[5].cutoff_radius == 1.5900000000E+00
        @test betas[5].beta[5] == -1.3905607575E-02

        @test betas[6].angular_momentum == 1
        @test betas[6].cutoff_radius_index == 160
        @test betas[6].cutoff_radius == 1.5900000000E+00
        @test betas[6].beta[5] == -1.3902597252E-02

        @test size(nonlocal.dij) == (6, 6)
        @test nonlocal.dij[1, 1] == -1.6525827200E+00
        @test nonlocal.dij[2, 2] == 4.4618201576E+00
        @test nonlocal.dij[3, 3] == -1.3035237938E+01
        @test nonlocal.dij[4, 4] == -1.2990770678E+01
        @test nonlocal.dij[5, 5] == -5.0193584007E+00
        @test nonlocal.dij[6, 6] == -4.9940955764E+00
        @test count(d -> !iszero(d), nonlocal.dij) == 6

        @test isnothing(nonlocal.augmentation)

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
            @test isnothing(chi.n)
            @test isnothing(chi.ultrasoft_cutoff_radius)
        end

        @test pswfc[1].label == "2S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.00
        @test pswfc[1].pseudo_energy == -0.5849619791E+01
        @test pswfc[1].chi[5] == 1.3903075777E-01

        @test pswfc[2].label == "2P"
        @test pswfc[2].l == 1
        @test pswfc[2].occupation == 4.00
        @test pswfc[2].pseudo_energy == -0.3413995091E+01
        @test pswfc[2].chi[5] == 1.2433653829E-02

        @test pswfc[4].label == "3S"
        @test pswfc[4].l == 0
        @test pswfc[4].occupation == 2.00
        @test pswfc[4].pseudo_energy == -0.3439210732E+00
        @test pswfc[4].chi[5] == 3.1825890104E-02

        @test isnothing(file.full_wfc)

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 4.1614005220E-02

        spin_orb = file.spin_orb

        @test length(spin_orb.relbetas) == 6
        for (i, relbeta) in enumerate(spin_orb.relbetas)
            @test relbeta.index == i
        end

        @test spin_orb.relbetas[1].lll == 0
        @test spin_orb.relbetas[1].jjj == 0.5

        @test spin_orb.relbetas[6].lll == 1
        @test spin_orb.relbetas[6].jjj == 1.5

        @test length(spin_orb.relwfcs) == 4
        for (i, relwfc) in enumerate(spin_orb.relwfcs)
            @test relwfc.index == i
        end

        @test spin_orb.relwfcs[1].lchi == 0
        @test spin_orb.relwfcs[1].jchi == 0.5
        @test spin_orb.relwfcs[1].nn == 1

        @test spin_orb.relwfcs[4].lchi == 0
        @test spin_orb.relwfcs[4].jchi == 0.5
        @test spin_orb.relwfcs[4].nn == 3

        @test isnothing(file.paw)
        @test isnothing(file.gipaw)
    end

    @testset "Si.pbe-n-rrkjus_psl.1.0.0.UPF" begin
        filename = "Si.pbe-n-rrkjus_psl.1.0.0.UPF"
        file = load_psp_file(UPF2_CASE_FILEPATHS[filename])

        header = file.header
        @test header.generated == "Generated using \"atomic\" code by A. Dal Corso  v.5.1"
        @test header.author == "ADC"
        @test header.date == "10Oct2014"
        @test header.comment == ""
        @test header.element == "Si"
        @test header.pseudo_type == "USPP"
        @test header.relativistic == "scalar"
        @test header.is_ultrasoft
        @test !header.is_paw
        @test !header.is_coulomb
        @test !header.has_so
        @test !header.has_wfc
        @test !header.has_gipaw
        @test !header.paw_as_gipaw
        @test header.core_correction
        @test header.functional == "PBE"
        @test header.z_valence == 4.000000000000000E+000
        @test header.total_psenergy == -1.102230803678973E+001
        @test header.wfc_cutoff == 4.374353160781668E+001
        @test header.rho_cutoff == 1.749741264312667E+002
        @test header.l_max == 2
        @test header.l_max_rho == 4
        @test header.l_local == -1
        @test header.mesh_size == 1141
        @test header.number_of_wfc == 2
        @test header.number_of_proj == 6

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 6.847393954957285E-005
        @test mesh.rab[5] == 8.559242443696606E-007
        @test mesh.mesh == 1141
        @test mesh.dx == 1.250000000000000E-002
        @test mesh.xmin == -7.000000000000000E+000
        @test mesh.rmax == 1.000000000000000E+002
        @test mesh.zmesh == 1.400000000000000E+001

        @test length(file.nlcc) == file.header.mesh_size
        @test file.nlcc[5] == 1.538616848020000E+000

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -8.438188853780485E+000

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj
        for (i, beta) in enumerate(betas)
            @test beta.index == i
            @test (length(beta.beta) == beta.cutoff_radius_index) | (length(beta.beta) == file.header.mesh_size)
        end

        @test betas[1].label == "3S"
        @test betas[1].angular_momentum == 0
        @test betas[1].cutoff_radius_index == 829
        @test betas[1].cutoff_radius == 1.600000000000000E+000
        @test betas[1].ultrasoft_cutoff_radius == 1.800000000000000E+000
        @test betas[1].beta[5] == -1.233266788373589E-003

        @test betas[2].label == "3S"
        @test betas[2].angular_momentum == 0
        @test betas[2].cutoff_radius_index == 829
        @test betas[2].cutoff_radius == 1.600000000000000E+000
        @test betas[2].ultrasoft_cutoff_radius == 1.800000000000000E+000
        @test betas[2].beta[5] == 1.539247050443425E-003

        @test betas[6].label == "3D"
        @test betas[6].angular_momentum == 2
        @test betas[6].cutoff_radius_index == 829
        @test betas[6].cutoff_radius == 1.600000000000000E+000
        @test betas[6].ultrasoft_cutoff_radius == 1.800000000000000E+000
        @test betas[6].beta[5] == 2.641408572045704E-010

        @test size(nonlocal.dij) == (6, 6)
        @test nonlocal.dij[1, 1] == 4.959233384252660E-001
        @test nonlocal.dij[2, 1] == 4.451348006463387E-001
        @test nonlocal.dij[1, 2] == 4.451348006463387E-001
        @test nonlocal.dij[2, 2] == 1.035512382098274E+000
        @test nonlocal.dij[3, 3] == 1.663108600181353E-001
        @test nonlocal.dij[3, 4] == 1.141776366606293E-001
        @test nonlocal.dij[4, 3] == 1.141776366606293E-001
        @test nonlocal.dij[4, 4] == 2.057883145155203E-001
        @test nonlocal.dij[5, 5] == -4.479438914118474E-002
        @test nonlocal.dij[5, 6] == -4.428967719198919E-002
        @test nonlocal.dij[6, 5] == -4.428967719198919E-002
        @test nonlocal.dij[6, 6] == -4.357014627662627E-002
        @test count(d -> !iszero(d), nonlocal.dij) == 12

        augmentation = nonlocal.augmentation
        @test augmentation.q_with_l
        @test augmentation.nqf == 0
        @test augmentation.nqlc == 5

        @test augmentation.q[1, 1] == -7.715557295409126E-002
        @test augmentation.q[2, 1] == -6.659970266131190E-002
        @test augmentation.q[1, 2] == -6.659970266131190E-002
        @test augmentation.q[2, 2] == -4.191304554821577E-002
        @test augmentation.q[3, 3] == -9.532064592377364E-003
        @test augmentation.q[3, 4] == -2.494352700293546E-002
        @test augmentation.q[4, 3] == -2.494352700293546E-002
        @test augmentation.q[4, 4] == -5.614360081862774E-002
        @test augmentation.q[5, 5] == 5.892629052420917E-003
        @test augmentation.q[5, 6] == 6.456888217977823E-003
        @test augmentation.q[6, 5] == 6.456888217977823E-003
        @test augmentation.q[6, 6] == 7.046059441924070E-003

        @test isnothing(augmentation.qijs)

        qijls = augmentation.qijls

        for (i, qijl) in enumerate(qijls)
            @test length(qijl.qijl) == file.header.mesh_size
        end

        @test qijls[1].first_index == 1
        @test qijls[1].second_index == 1
        @test qijls[1].composite_index == 1
        @test qijls[1].angular_momentum == 0
        @test qijls[1].qijl[5] == 7.047762479021314E-010

        @test qijls[2].first_index == 1
        @test qijls[2].second_index == 2
        @test qijls[2].composite_index == 2
        @test qijls[2].angular_momentum == 0
        @test qijls[2].qijl[5] == -3.457052002690957E-009

        @test qijls[34].first_index == 6
        @test qijls[34].second_index == 6
        @test qijls[34].composite_index == 21
        @test qijls[34].angular_momentum == 4
        @test qijls[34].qijl[5] == 7.271904089448543E-026

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
        end

        @test pswfc[1].label == "3S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.000000000000000E+000
        @test pswfc[1].n == 1
        @test pswfc[1].pseudo_energy == -7.947277707288829E-001
        @test pswfc[1].cutoff_radius == 1.600000000000000E+000
        @test pswfc[1].ultrasoft_cutoff_radius == 1.800000000000000E+000
        @test pswfc[1].chi[5] == 3.542904978471305E-005

        @test pswfc[2].label == "3P"
        @test pswfc[2].l == 1
        @test pswfc[2].occupation == 2.000000000000000E+000
        @test pswfc[2].n == 2
        @test pswfc[2].pseudo_energy == -2.999646802634566E-001
        @test pswfc[2].cutoff_radius == 1.600000000000000E+000
        @test pswfc[2].ultrasoft_cutoff_radius == 1.800000000000000E+000
        @test pswfc[2].chi[5] == 2.067383051207364E-009

        @test isnothing(file.full_wfc)

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 6.986543942450972E-009

        @test isnothing(file.spin_orb)
        @test isnothing(file.paw)
        @test isnothing(file.gipaw)
    end

    @testset "Al.pbe-n-kjpaw_psl.1.0.0.UPF" begin
        filename = "Al.pbe-n-kjpaw_psl.1.0.0.UPF"
        file = load_psp_file(UPF2_CASE_FILEPATHS[filename])

        header = file.header
        @test header.generated == "Generated using \"atomic\" code by A. Dal Corso  v.5.1"
        @test header.author == "ADC"
        @test header.date == "10Oct2014"
        @test header.comment == ""
        @test header.element == "Al"
        @test header.pseudo_type == "PAW"
        @test header.relativistic == "scalar"
        @test header.is_ultrasoft
        @test header.is_paw
        @test !header.is_coulomb
        @test !header.has_so
        @test header.has_wfc
        @test header.has_gipaw
        @test header.paw_as_gipaw
        @test header.core_correction
        @test header.functional == "SLA PW PBX PBC"
        @test header.z_valence == 3.000000000000000E+000
        @test header.total_psenergy == -3.922759814471728E+001
        @test header.wfc_cutoff == 2.949430082155506E+001
        @test header.rho_cutoff == 1.434836740307051E+002
        @test header.l_max == 2
        @test header.l_max_rho == 4
        @test header.l_local == -1
        @test header.mesh_size == 1135
        @test header.number_of_wfc == 2
        @test header.number_of_proj == 6

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 7.374116566877076E-005
        @test mesh.rab[5] == 9.217645708596345E-007
        @test mesh.mesh == 1135
        @test mesh.dx == 1.250000000000000E-002
        @test mesh.xmin == -7.000000000000000E+000
        @test mesh.rmax == 1.000000000000000E+002
        @test mesh.zmesh == 1.300000000000000E+001

        @test length(file.nlcc) == file.header.mesh_size
        @test file.nlcc[5] == 2.656001002192721E-001

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -6.296496837141098E+000

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj
        for (i, beta) in enumerate(betas)
            @test beta.index == i
            @test (length(beta.beta) == beta.cutoff_radius_index) | (length(beta.beta) == file.header.mesh_size)
        end

        @test betas[1].label == "3S"
        @test betas[1].angular_momentum == 0
        @test betas[1].cutoff_radius_index == 827
        @test betas[1].cutoff_radius == 1.700000000000000E+000
        @test betas[1].ultrasoft_cutoff_radius == 1.900000000000000E+000
        @test betas[1].beta[5] == -1.792962988082909E-003

        @test betas[2].label == "3S"
        @test betas[2].angular_momentum == 0
        @test betas[2].cutoff_radius_index == 827
        @test betas[2].cutoff_radius == 1.700000000000000E+000
        @test betas[2].ultrasoft_cutoff_radius == 1.900000000000000E+000
        @test betas[2].beta[5] == 1.737827866806716E-003

        @test betas[6].label == "3D"
        @test betas[6].angular_momentum == 2
        @test betas[6].cutoff_radius_index == 827
        @test betas[6].cutoff_radius == 1.700000000000000E+000
        @test betas[6].ultrasoft_cutoff_radius == 1.900000000000000E+000
        @test betas[6].beta[5] == 2.699493457967658E-010

        @test size(nonlocal.dij) == (6, 6)
        @test nonlocal.dij[1, 1] == 4.162665025130763E-001
        @test nonlocal.dij[2, 1] == 5.099080954670963E-001
        @test nonlocal.dij[1, 2] == 5.099080954670963E-001
        @test nonlocal.dij[2, 2] == 1.134914583776183E+000
        @test nonlocal.dij[3, 3] == 1.526405720963211E-001
        @test nonlocal.dij[3, 4] == 2.252999736729244E-001
        @test nonlocal.dij[4, 3] == 2.252999736729244E-001
        @test nonlocal.dij[4, 4] == 3.368020928299262E-001
        @test nonlocal.dij[5, 5] == -4.599057251334548E-002
        @test nonlocal.dij[5, 6] == -4.575373726218579E-002
        @test nonlocal.dij[6, 5] == -4.575373726218579E-002
        @test nonlocal.dij[6, 6] == -4.534796520906204E-002
        @test count(d -> !iszero(d), nonlocal.dij) == 12

        augmentation = nonlocal.augmentation
        @test augmentation.q_with_l
        @test augmentation.nqf == 0
        @test augmentation.nqlc == 5
        @test augmentation.shape == "PSQ"
        @test augmentation.cutoff_r == -1.000000000000000E+000
        @test augmentation.cutoff_r_index == 839
        @test augmentation.augmentation_epsilon == 1.000000000000000E-012
        @test augmentation.l_max_aug == 4

        @test augmentation.q[1, 1] == -3.070968103524502E-002
        @test augmentation.q[2, 1] == -3.491452538270029E-002
        @test augmentation.q[1, 2] == -3.491452538270029E-002
        @test augmentation.q[2, 2] == -2.379025736223055E-002
        @test augmentation.q[3, 3] == -2.007950908549670E-003
        @test augmentation.q[3, 4] == -1.350416114335737E-002
        @test augmentation.q[4, 3] == -1.350416114335737E-002
        @test augmentation.q[4, 4] == -4.115951195614498E-002
        @test augmentation.q[5, 5] == 5.002993776692934E-003
        @test augmentation.q[5, 6] == 5.470781396762834E-003
        @test augmentation.q[6, 5] == 5.470781396762834E-003
        @test augmentation.q[6, 6] == 5.960693599449139E-003

        @test augmentation.multipoles[7] == -3.491452538270029E-002

        @test isnothing(augmentation.qijs)

        qijls = augmentation.qijls

        for (i, qijl) in enumerate(qijls)
            @test length(qijl.qijl) == file.header.mesh_size
        end

        @test qijls[1].first_index == 1
        @test qijls[1].second_index == 1
        @test qijls[1].composite_index == 1
        @test qijls[1].angular_momentum == 0
        @test qijls[1].qijl[5] == 1.301443338375301E-009

        @test qijls[2].first_index == 1
        @test qijls[2].second_index == 2
        @test qijls[2].composite_index == 2
        @test qijls[2].angular_momentum == 0
        @test qijls[2].qijl[5] == -1.862965770722507E-009

        @test qijls[34].first_index == 6
        @test qijls[34].second_index == 6
        @test qijls[34].composite_index == 21
        @test qijls[34].angular_momentum == 4
        @test qijls[34].qijl[5] == 6.608621564156752E-026

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
        end

        @test pswfc[1].label == "3S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.000000000000000E+000
        @test pswfc[1].n == 1
        @test pswfc[1].pseudo_energy == -5.698258955729895E-001
        @test pswfc[1].cutoff_radius == 1.700000000000000E+000
        @test pswfc[1].ultrasoft_cutoff_radius == 1.900000000000000E+000
        @test pswfc[1].chi[5] == 2.099526933006344E-005

        @test pswfc[2].label == "3P"
        @test pswfc[2].l == 1
        @test pswfc[2].occupation == 1.000000000000000E+000
        @test pswfc[2].n == 2
        @test pswfc[2].pseudo_energy == -1.993375352101455E-001
        @test pswfc[2].cutoff_radius == 1.700000000000000E+000
        @test pswfc[2].ultrasoft_cutoff_radius == 2.000000000000000E+000
        @test pswfc[2].chi[5] == 1.330682509144020E-009

        full_wfc = file.full_wfc

        @assert length(full_wfc.aewfcs) == 6
        for (i, aewfc) in enumerate(full_wfc.aewfcs)
            @test aewfc.index == i
        end

        @test full_wfc.aewfcs[1].label == "3S"
        @test full_wfc.aewfcs[1].l == 0
        @test full_wfc.aewfcs[1].wfc[5] == 4.609881930168180E-004

        @test full_wfc.aewfcs[2].label == "3S"
        @test full_wfc.aewfcs[2].l == 0
        @test full_wfc.aewfcs[2].wfc[5] == 7.371129588893765E-004

        @test full_wfc.aewfcs[6].label == "3D"
        @test full_wfc.aewfcs[6].l == 2
        @test full_wfc.aewfcs[6].wfc[5] == 2.534444746026031E-012

        @assert length(full_wfc.pswfcs) == 6
        for (i, pswfc) in enumerate(full_wfc.pswfcs)
            @test pswfc.index == i
        end

        @test full_wfc.pswfcs[1].label == "3S"
        @test full_wfc.pswfcs[1].l == 0
        @test full_wfc.pswfcs[1].wfc[5] == 2.099528103374044E-005

        @test full_wfc.pswfcs[2].label == "3S"
        @test full_wfc.pswfcs[2].l == 0
        @test full_wfc.pswfcs[2].wfc[5] == 1.540244976023415E-005

        @test full_wfc.pswfcs[6].label == "3D"
        @test full_wfc.pswfcs[6].l == 2
        @test full_wfc.pswfcs[6].wfc[5] == 1.794811387770046E-013

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 4.550846706828015E-009

        @test isnothing(file.spin_orb)

        paw = file.paw
        @test paw.paw_data_format == 2
        @test paw.core_energy == -4.461427815823291E+002

        @test length(paw.occupations) == 6
        @test paw.occupations[3] == 1.000000000000000E+000

        @test length(paw.ae_nlcc) == file.header.mesh_size
        @test paw.ae_nlcc[5] == 1.490091617040069E+003

        @test length(paw.ae_vloc) == file.header.mesh_size
        @test paw.ae_vloc[5] == -3.524982145519007E+005

        gipaw = file.gipaw
        @test gipaw.gipaw_data_format == 2

        @test length(gipaw.core_orbitals) == 3
        for (i, core_orbital) in enumerate(gipaw.core_orbitals)
            @test length(core_orbital.core_orbital) == file.header.mesh_size
            @test core_orbital.index == i
        end

        @test gipaw.core_orbitals[1].label == "1S"
        @test gipaw.core_orbitals[1].l == 0
        @test gipaw.core_orbitals[1].n == 1
        @test gipaw.core_orbitals[1].core_orbital[5] == 6.905763114231386E-003

        @test gipaw.core_orbitals[2].label == "2S"
        @test gipaw.core_orbitals[2].l == 0
        @test gipaw.core_orbitals[2].n == 2
        @test gipaw.core_orbitals[2].core_orbital[5] == 1.794874938233875E-003

        @test gipaw.core_orbitals[3].label == "2P"
        @test gipaw.core_orbitals[3].l == 1
        @test gipaw.core_orbitals[3].n == 2
        @test gipaw.core_orbitals[3].core_orbital[5] == 7.103427686008862E-007
    end

    @testset "He.pbe-hgh.UPF" begin
        filename = "He.pbe-hgh.UPF"
        file = load_psp_file(UPF2_CASE_FILEPATHS[filename])

        header = file.header
        @test header.generated == "Generated in analytical, separable form"
        @test header.author == "Goedecker/Hartwigsen/Hutter/Teter"
        @test header.date == "Phys.Rev.B58, 3641 (1998); B54, 1703 (1996)"
        @test header.comment ==
              "Contains atomic orbitals generated by ld1.x - use with care"
        @test header.element == "He"
        @test header.pseudo_type == "NC"
        @test header.relativistic == "scalar"
        @test !header.is_ultrasoft
        @test !header.is_paw
        @test !header.is_coulomb
        @test !header.has_so
        @test !header.has_wfc
        @test !header.has_gipaw
        @test !header.paw_as_gipaw
        @test !header.core_correction
        @test header.functional == "SLA PW PBX PBC"
        @test header.z_valence == 2.000000000000000E+000
        @test header.total_psenergy == 0.000000000000000E+000
        @test header.wfc_cutoff == 0.000000000000000E+000
        @test header.rho_cutoff == 0.000000000000000E+000
        @test header.l_max == -1
        @test header.l_max_rho == 0
        @test header.l_local == -3
        @test header.mesh_size == 985
        @test header.number_of_wfc == 1
        @test header.number_of_proj == 0

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 4.793175768470099E-004
        @test mesh.rab[5] == 5.991469710587625E-006
        @test mesh.mesh == 985
        @test mesh.dx == 1.250000000000000E-002
        @test mesh.xmin == -7.000000000000000E+000
        @test mesh.rmax == 1.001684049873959E+002
        @test mesh.zmesh == 2.000000000000000E+000

        @test isnothing(file.nlcc)

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -3.420189165663822E+001

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj

        @test size(nonlocal.dij) == (0, 0)

        @test isnothing(nonlocal.augmentation)

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
        end

        @test pswfc[1].label == "1S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.000000000000000E+000
        @test pswfc[1].n == 1
        @test pswfc[1].pseudo_energy == -1.158234898214897E+000
        @test pswfc[1].cutoff_radius == 6.154897413020295E-001
        @test pswfc[1].ultrasoft_cutoff_radius == 6.154897413020295E-001
        @test pswfc[1].chi[5] == 1.831401620367775E-003

        @test isnothing(file.full_wfc)

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 6.708063790171425E-006

        @test isnothing(file.spin_orb)
        @test isnothing(file.paw)
        @test isnothing(file.gipaw)
    end
end
