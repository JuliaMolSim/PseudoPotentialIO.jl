@testset "UPF v1.old" begin
    @testset "Internal data consistency" begin
        for filepath in values(UPF1_CASE_FILEPATHS)
            file = load_psp_file(filepath)

            @test format(file) == "UPF v1.old"

            # UPF v1.old has a different augmentation data format from UPF v2.0.1
            if file.header.is_ultrasoft | file.header.is_paw
                augmentation = file.nonlocal.augmentation

                @test isnothing(augmentation.multipoles)

                if augmentation.nqf > 0
                    @test !isnothing(augmentation.qfcoefs)
                    @test length(augmentation.qfcoefs) ==
                          (file.header.number_of_proj * (file.header.number_of_proj + 1) /
                           2)
                    for qfcoef in augmentation.qfcoefs
                        @test length(qfcoef.qfcoef) ==
                              augmentation.nqf * (2file.header.l_max + 1)
                    end
                else
                    @test isnothing(augmentation.qfcoefs)
                end

                @test !isnothing(augmentation.rinner)
                @test length(augmentation.rinner) == 2file.header.l_max + 1
            end
        end
    end

    @testset "ag_lda_v1.4.uspp.F.UPF" begin
        filename = "ag_lda_v1.4.uspp.F.UPF"
        file = load_psp_file(UPF1_CASE_FILEPATHS[filename])

        header = file.header
        @test isnothing(header.generated)
        @test isnothing(header.author)
        @test isnothing(header.date)
        @test isnothing(header.comment)
        @test header.element == "Ag"
        @test header.pseudo_type == "US"
        @test header.relativistic == "scalar"
        @test header.is_ultrasoft
        @test !header.is_paw
        @test !header.is_coulomb
        @test !header.has_so
        @test !header.has_wfc
        @test !header.has_gipaw
        @test isnothing(header.paw_as_gipaw)
        @test header.core_correction
        @test header.functional == "SLA PZ NOGX NOGC PZ"
        @test header.z_valence == 19.00000000000
        @test header.total_psenergy == -295.20667856400
        @test header.wfc_cutoff == 0.00000
        @test header.rho_cutoff == 0.00000
        @test header.l_max == 2
        @test isnothing(header.l_max_rho)
        @test isnothing(header.l_local)
        @test header.mesh_size == 887
        @test header.number_of_wfc == 5
        @test header.number_of_proj == 6

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 5.05789189118E-09
        @test mesh.rab[5] == 1.32875017333E-09
        @test isnothing(mesh.mesh)
        @test isnothing(mesh.dx)
        @test isnothing(mesh.xmin)
        @test isnothing(mesh.rmax)
        @test isnothing(mesh.zmesh)

        @test length(file.nlcc) == file.header.mesh_size
        @test file.nlcc[5] == 4.66361704568E-01

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -5.08978714579E+01

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj
        for (i, beta) in enumerate(betas)
            @test beta.index == i
            @test isnothing(beta.cutoff_radius)
            @test (length(beta.beta) == beta.cutoff_radius_index) | (length(beta.beta) == file.header.mesh_size)
            @test isnothing(beta.norm_conserving_radius)
            @test isnothing(beta.ultrasoft_cutoff_radius)
            @test isnothing(beta.label)
        end

        @test betas[1].angular_momentum == 0
        @test betas[1].cutoff_radius_index == 717
        @test betas[1].beta[5] == -5.63186386318E-05

        @test betas[2].angular_momentum == 0
        @test betas[2].cutoff_radius_index == 717
        @test betas[2].beta[5] == -1.63507161561E-04

        @test betas[3].angular_momentum == 1
        @test betas[3].cutoff_radius_index == 717
        @test betas[3].beta[5] == -2.79143566540E-08

        @test betas[4].angular_momentum == 1
        @test betas[4].cutoff_radius_index == 717
        @test betas[4].beta[5] == -4.68355152378E-09

        @test betas[5].angular_momentum == 2
        @test betas[5].cutoff_radius_index == 717
        @test betas[5].beta[5] == -1.89086716506E-14

        @test betas[6].angular_momentum == 2
        @test betas[6].cutoff_radius_index == 717
        @test betas[6].beta[5] == -3.04887580518E-13

        @test size(nonlocal.dij) == (6, 6)
        @test nonlocal.dij[1, 1] == -1.76074960113E+00
        @test nonlocal.dij[1, 2] == 1.92857178205E-01
        @test nonlocal.dij[2, 2] == -1.22330466146E+00
        @test nonlocal.dij[3, 3] == 3.13785928053E+00
        @test nonlocal.dij[3, 4] == -9.84909312467E+00
        @test nonlocal.dij[4, 4] == 1.83276570475E+01
        @test nonlocal.dij[5, 5] == 4.83747588237E+00
        @test nonlocal.dij[5, 6] == 1.17462592708E+01
        @test nonlocal.dij[6, 6] == 2.56599483603E+01
        @test count(d -> !iszero(d), nonlocal.dij) == 12

        augmentation = nonlocal.augmentation
        @test isnothing(augmentation.multipoles)
        @test !augmentation.q_with_l
        @test augmentation.nqf == 6

        @test size(augmentation.q) == (6, 6)
        @test augmentation.q[1, 1] == -3.56687909692E-01
        @test augmentation.q[1, 2] == augmentation.q[2, 1] == 2.26500592008E-01
        @test augmentation.q[6, 6] == 1.84734674885E+00

        @test length(augmentation.qfcoefs) == 21
        @test augmentation.qfcoefs[1].qfcoef[5] == -3.73585933729E+00
        @test augmentation.qfcoefs[2].qfcoef[5] == 1.37859649853E+00

        @test all(r -> r == 1.30000000000E+00, augmentation.rinner)

        @test length(augmentation.qijs) == 21
        for qij in augmentation.qijs
            @test length(qij.qij) == file.header.mesh_size
        end
        @test augmentation.qijs[1].first_index == 1
        @test augmentation.qijs[1].second_index == 1
        @test augmentation.qijs[1].qij[5] == -1.55354181933E-16
        @test augmentation.qijs[2].first_index == 1
        @test augmentation.qijs[2].second_index == 2
        @test augmentation.qijs[2].qij[5] == 7.90456435011E-17

        @test isnothing(augmentation.qijls)
        @test !augmentation.q_with_l
        @test augmentation.nqf == 6
        @test isnothing(augmentation.nqlc)
        @test isnothing(augmentation.shape)
        @test isnothing(augmentation.iraug)
        @test isnothing(augmentation.raug)
        @test isnothing(augmentation.l_max_aug)
        @test isnothing(augmentation.augmentation_epsilon)
        @test isnothing(augmentation.cutoff_r)
        @test isnothing(augmentation.cutoff_r_index)

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
            @test isnothing(chi.n)
            @test isnothing(chi.pseudo_energy)
            @test isnothing(chi.cutoff_radius)
            @test isnothing(chi.ultrasoft_cutoff_radius)
        end

        @test pswfc[1].label == "4S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.00
        @test pswfc[1].chi[5] == 1.53041946142E-08

        @test pswfc[2].label == "4P"
        @test pswfc[2].l == 1
        @test pswfc[2].occupation == 6.00
        @test pswfc[2].chi[5] == 7.18927166524E-17

        @test pswfc[5].label == "5P"
        @test pswfc[5].l == 1
        @test pswfc[5].occupation == 0.00
        @test pswfc[5].chi[5] == 1.48358454137E-17

        @test isnothing(file.full_wfc)

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 2.89294930407E-16

        @test isnothing(file.spin_orb)
        @test isnothing(file.paw)
        @test isnothing(file.gipaw)
    end

    @testset "B_pbe_v1.01.uspp.F.UPF" begin
        filename = "B_pbe_v1.01.uspp.F.UPF"
        file = load_psp_file(UPF1_CASE_FILEPATHS[filename])

        header = file.header
        @test isnothing(header.generated)
        @test isnothing(header.author)
        @test isnothing(header.date)
        @test isnothing(header.comment)
        @test header.element == "B"
        @test header.pseudo_type == "US"
        @test header.relativistic == "scalar"
        @test header.is_ultrasoft
        @test !header.is_paw
        @test !header.is_coulomb
        @test !header.has_so
        @test !header.has_wfc
        @test !header.has_gipaw
        @test isnothing(header.paw_as_gipaw)
        @test header.core_correction
        @test header.functional == "SLA PW PBE PBE PBE"
        @test header.z_valence == 3.00000000000
        @test header.total_psenergy == -6.00744913082
        @test header.wfc_cutoff == 0.00000
        @test header.rho_cutoff == 0.00000
        @test header.l_max == 1
        @test isnothing(header.l_max_rho)
        @test isnothing(header.l_local)
        @test header.mesh_size == 781
        @test header.number_of_wfc == 2
        @test header.number_of_proj == 4

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 1.25728654505E-05
        @test mesh.rab[5] == 3.24915430936E-06
        @test isnothing(mesh.mesh)
        @test isnothing(mesh.dx)
        @test isnothing(mesh.xmin)
        @test isnothing(mesh.rmax)
        @test isnothing(mesh.zmesh)

        @test length(file.nlcc) == file.header.mesh_size
        @test file.nlcc[5] == 1.20782515483E+00

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -1.00137766482E+01

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj
        for (i, beta) in enumerate(betas)
            @test beta.index == i
            @test isnothing(beta.cutoff_radius)
            @test (length(beta.beta) == beta.cutoff_radius_index) | (length(beta.beta) == file.header.mesh_size)
            @test isnothing(beta.norm_conserving_radius)
            @test isnothing(beta.ultrasoft_cutoff_radius)
            @test isnothing(beta.label)
        end

        @test betas[1].angular_momentum == 0
        @test betas[1].cutoff_radius_index == 559
        @test betas[1].beta[5] == -3.82483956865E-05

        @test betas[2].angular_momentum == 0
        @test betas[2].cutoff_radius_index == 559
        @test betas[2].beta[5] == -8.67580981562E-05

        @test betas[3].angular_momentum == 1
        @test betas[3].cutoff_radius_index == 559
        @test betas[3].beta[5] == 1.99207032008E-07

        @test betas[4].angular_momentum == 1
        @test betas[4].cutoff_radius_index == 559
        @test betas[4].beta[5] == 1.69150912265E-07

        @test size(nonlocal.dij) == (4, 4)
        @test nonlocal.dij[1, 1] == -1.22873815109E-01
        @test nonlocal.dij[1, 2] == -5.61912458144E+00
        @test nonlocal.dij[2, 2] == -5.87789890874E+00
        @test nonlocal.dij[3, 3] == 9.71393906034E+00
        @test nonlocal.dij[3, 4] == 2.46444027912E+01
        @test nonlocal.dij[4, 4] == 6.13899360806E+01
        @test count(d -> !iszero(d), nonlocal.dij) == 8

        augmentation = nonlocal.augmentation
        @test isnothing(augmentation.multipoles)
        @test !augmentation.q_with_l
        @test augmentation.nqf == 8

        @test size(augmentation.q) == (4, 4)
        @test augmentation.q[1, 1] == -5.06674546795E-01
        @test augmentation.q[1, 2] == augmentation.q[2, 1] == -3.23289503777E-01
        @test augmentation.q[1, 3] == augmentation.q[3, 1] == 0.00000000000E+00
        @test augmentation.q[1, 4] == augmentation.q[4, 1] == 0.00000000000E+00
        @test augmentation.q[2, 2] == 1.89884871724E-01
        @test augmentation.q[2, 3] == augmentation.q[2, 3] == 0.00000000000E+00
        @test augmentation.q[2, 4] == augmentation.q[4, 2] == 0.00000000000E+00
        @test augmentation.q[3, 3] == 9.88682989485E-01
        @test augmentation.q[4, 3] == augmentation.q[3, 4] == 1.83449812909E+00
        @test augmentation.q[4, 4] == 3.27504485114E+00

        @test length(augmentation.rinner) == 3
        @test all(r -> r == 1.10000000000E+00, augmentation.rinner)

        @test length(augmentation.qijs) == 10
        for qij in augmentation.qijs
            @test length(qij.qij) == file.header.mesh_size
        end
        @test length(augmentation.qfcoefs) == 10

        @test augmentation.qijs[1].first_index == 1
        @test augmentation.qijs[1].second_index == 1
        @test augmentation.qijs[1].qij[5] == -1.94859515081E-09
        @test augmentation.qfcoefs[1].qfcoef[5] == -7.21237350451E+01

        @test augmentation.qijs[2].first_index == 1
        @test augmentation.qijs[2].second_index == 2
        @test augmentation.qijs[2].qij[5] == -1.53409489684E-09
        @test augmentation.qfcoefs[2].qfcoef[5] == -5.60931726792E+01

        @test augmentation.qijs[10].first_index == 4
        @test augmentation.qijs[10].second_index == 4
        @test augmentation.qijs[10].qij[5] == 1.57362794452E-08
        @test augmentation.qfcoefs[10].qfcoef[5] == 7.28024914160E+02

        @test isnothing(augmentation.qijls)
        @test isnothing(augmentation.nqlc)
        @test isnothing(augmentation.shape)
        @test isnothing(augmentation.iraug)
        @test isnothing(augmentation.raug)
        @test isnothing(augmentation.l_max_aug)
        @test isnothing(augmentation.augmentation_epsilon)
        @test isnothing(augmentation.cutoff_r)
        @test isnothing(augmentation.cutoff_r_index)

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
            @test isnothing(chi.n)
            @test isnothing(chi.pseudo_energy)
            @test isnothing(chi.cutoff_radius)
            @test isnothing(chi.ultrasoft_cutoff_radius)
        end

        @test pswfc[1].label == "2S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.00
        @test pswfc[1].chi[5] == 5.96294641979E-06

        @test pswfc[2].label == "2P"
        @test pswfc[2].l == 1
        @test pswfc[2].occupation == 1.00
        @test pswfc[2].chi[5] == 1.64367609052E-10

        @test isnothing(file.full_wfc)

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 1.04916913718E-10

        @test isnothing(file.spin_orb)
        @test isnothing(file.paw)
        @test isnothing(file.gipaw)
    end

    @testset "si_pbesol_v1.uspp.F.UPF" begin
        filename = "si_pbesol_v1.uspp.F.UPF"
        file = load_psp_file(UPF1_CASE_FILEPATHS[filename])

        header = file.header
        @test isnothing(header.generated)
        @test isnothing(header.author)
        @test isnothing(header.date)
        @test isnothing(header.comment)
        @test header.element == "Si"
        @test header.pseudo_type == "US"
        @test header.relativistic == "scalar"
        @test header.is_ultrasoft
        @test !header.is_paw
        @test !header.is_coulomb
        @test !header.has_so
        @test !header.has_wfc
        @test !header.has_gipaw
        @test isnothing(header.paw_as_gipaw)
        @test header.core_correction
        @test header.functional == "SLA PW PSX PSC PBEsol"
        @test header.z_valence == 4.00000000000
        @test header.total_psenergy == -9.13712037618
        @test header.wfc_cutoff == 0.00000
        @test header.rho_cutoff == 0.00000
        @test header.l_max == 2
        @test isnothing(header.l_max_rho)
        @test isnothing(header.l_local)
        @test header.mesh_size == 899
        @test header.number_of_wfc == 2
        @test header.number_of_proj == 6

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 4.49030908945E-06
        @test mesh.rab[5] == 1.16041225334E-06
        @test isnothing(mesh.mesh)
        @test isnothing(mesh.dx)
        @test isnothing(mesh.xmin)
        @test isnothing(mesh.rmax)
        @test isnothing(mesh.zmesh)

        @test length(file.nlcc) == file.header.mesh_size
        @test file.nlcc[5] == 6.79372951207E-01

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -1.02787851219E+01

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj
        for (i, beta) in enumerate(betas)
            @test beta.index == i
            @test isnothing(beta.cutoff_radius)
            @test (length(beta.beta) == beta.cutoff_radius_index) | (length(beta.beta) == file.header.mesh_size)
            @test isnothing(beta.norm_conserving_radius)
            @test isnothing(beta.ultrasoft_cutoff_radius)
            @test isnothing(beta.label)
        end

        @test betas[1].angular_momentum == 0
        @test betas[1].cutoff_radius_index == 627
        @test betas[1].beta[5] == -8.50729690678E-06

        @test betas[2].angular_momentum == 0
        @test betas[2].cutoff_radius_index == 627
        @test betas[2].beta[5] == 2.10924706615E-05

        @test betas[3].angular_momentum == 1
        @test betas[3].cutoff_radius_index == 627
        @test betas[3].beta[5] == 3.11188017130E-09

        @test betas[4].angular_momentum == 1
        @test betas[4].cutoff_radius_index == 627
        @test betas[4].beta[5] == -1.50732607827E-08

        @test betas[5].angular_momentum == 2
        @test betas[5].cutoff_radius_index == 627
        @test betas[5].beta[5] == 1.13855763125E-11

        @test betas[6].angular_momentum == 2
        @test betas[6].cutoff_radius_index == 627
        @test betas[6].beta[5] == -4.89183257166E-11

        @test size(nonlocal.dij) == (6, 6)
        @test nonlocal.dij[1, 1] == 4.17504606973E-01
        @test nonlocal.dij[1, 2] == 4.05461929704E+00
        @test nonlocal.dij[2, 2] == -3.34190790394E+00
        @test nonlocal.dij[3, 3] == -1.85637785661E+00
        @test nonlocal.dij[3, 4] == 4.31799349543E+00
        @test nonlocal.dij[4, 4] == -1.54087360331E+00
        @test nonlocal.dij[5, 5] == 3.24729498129E+00
        @test nonlocal.dij[5, 6] == -5.46143064005E+00
        @test nonlocal.dij[6, 6] == 9.04320814597E+00
        @test count(d -> !iszero(d), nonlocal.dij) == 12

        augmentation = nonlocal.augmentation
        @test isnothing(augmentation.multipoles)
        @test !augmentation.q_with_l
        @test augmentation.nqf == 8

        @test size(augmentation.q) == (6, 6)
        @test augmentation.q[1, 1] == -5.99456496341E-01
        @test augmentation.q[1, 2] == augmentation.q[2, 1] == 4.20874806979E-01
        @test augmentation.q[6, 6] == 6.19942163555E-01

        @test length(augmentation.rinner) == 5
        @test all(r -> r == 9.00000000000E-01, augmentation.rinner)

        @test length(augmentation.qijs) == 21
        for qij in augmentation.qijs
            @test length(qij.qij) == file.header.mesh_size
        end
        @test length(augmentation.qfcoefs) == 21
        
        @test augmentation.qijs[1].first_index == 1
        @test augmentation.qijs[1].second_index == 1
        @test augmentation.qijs[1].qij[5] == 2.30032305636E-12
        @test augmentation.qfcoefs[1].qfcoef[5] == -2.28802388683E+02

        @test augmentation.qijs[2].first_index == 1
        @test augmentation.qijs[2].second_index == 2
        @test augmentation.qijs[2].qij[5] == 1.69389615913E-10
        @test augmentation.qfcoefs[2].qfcoef[5] == -7.09913657883E+01

        @test augmentation.qijs[21].first_index == 6
        @test augmentation.qijs[21].second_index == 6
        @test augmentation.qijs[21].qij[5] == 5.14863218160E-12
        @test augmentation.qfcoefs[21].qfcoef[5] == 1.02352656213E+02

        @test isnothing(augmentation.qijls)
        @test isnothing(augmentation.nqlc)
        @test isnothing(augmentation.shape)
        @test isnothing(augmentation.iraug)
        @test isnothing(augmentation.raug)
        @test isnothing(augmentation.l_max_aug)
        @test isnothing(augmentation.augmentation_epsilon)
        @test isnothing(augmentation.cutoff_r)
        @test isnothing(augmentation.cutoff_r_index)

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
            @test isnothing(chi.n)
            @test isnothing(chi.pseudo_energy)
            @test isnothing(chi.cutoff_radius)
            @test isnothing(chi.ultrasoft_cutoff_radius)
        end

        @test pswfc[1].label == "3S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.00
        @test pswfc[1].chi[5] == 6.44327877084E-07

        @test pswfc[2].label == "3P"
        @test pswfc[2].l == 1
        @test pswfc[2].occupation == 2.00
        @test pswfc[2].chi[5] == 6.73449169908E-12

        @test isnothing(file.full_wfc)

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 9.84444587528E-11

        @test isnothing(file.spin_orb)
        @test isnothing(file.paw)
        @test isnothing(file.gipaw)
    end

    @testset "mg_pbe_v1.4.uspp.F.UPF" begin
        filename = "mg_pbe_v1.4.uspp.F.UPF"
        file = load_psp_file(UPF1_CASE_FILEPATHS[filename])

        header = file.header
        @test isnothing(header.generated)
        @test isnothing(header.author)
        @test isnothing(header.date)
        @test isnothing(header.comment)
        @test header.element == "Mg"
        @test header.pseudo_type == "US"
        @test header.relativistic == "scalar"
        @test header.is_ultrasoft
        @test !header.is_paw
        @test !header.is_coulomb
        @test !header.has_so
        @test !header.has_wfc
        @test !header.has_gipaw
        @test isnothing(header.paw_as_gipaw)
        @test !header.core_correction
        @test header.functional == "SLA PW PBX PBC PBE"
        @test header.z_valence == 10.00000000000
        @test header.total_psenergy == -125.08981684000
        @test header.wfc_cutoff == 0.00000
        @test header.rho_cutoff == 0.00000
        @test header.l_max == 2
        @test isnothing(header.l_max_rho)
        @test isnothing(header.l_local)
        @test header.mesh_size == 821
        @test header.number_of_wfc == 3
        @test header.number_of_proj == 7

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 1.98100765738E-08
        @test mesh.rab[5] == 5.20427151222E-09
        @test isnothing(mesh.mesh)
        @test isnothing(mesh.dx)
        @test isnothing(mesh.xmin)
        @test isnothing(mesh.rmax)
        @test isnothing(mesh.zmesh)

        @test isnothing(file.nlcc)

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -3.96180531727E+01

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj
        for (i, beta) in enumerate(betas)
            @test beta.index == i
            @test isnothing(beta.cutoff_radius)
            @test (length(beta.beta) == beta.cutoff_radius_index) | (length(beta.beta) == file.header.mesh_size)
            @test isnothing(beta.norm_conserving_radius)
            @test isnothing(beta.ultrasoft_cutoff_radius)
            @test isnothing(beta.label)
        end

        @test betas[1].angular_momentum == 0
        @test betas[1].cutoff_radius_index == 661
        @test betas[1].beta[5] == -8.06372437750E-05

        @test betas[2].angular_momentum == 0
        @test betas[2].cutoff_radius_index == 661
        @test betas[2].beta[5] == 1.45732756410E-05

        @test betas[3].angular_momentum == 0
        @test betas[3].cutoff_radius_index == 661
        @test betas[3].beta[5] == 6.78200641134E-04

        @test betas[4].angular_momentum == 1
        @test betas[4].cutoff_radius_index == 661
        @test betas[4].beta[5] == -4.68700375660E-09

        @test betas[5].angular_momentum == 1
        @test betas[5].cutoff_radius_index == 661
        @test betas[5].beta[5] == 1.69113364129E-07

        @test betas[6].angular_momentum == 2
        @test betas[6].cutoff_radius_index == 661
        @test betas[6].beta[5] == -6.58196838008E-14

        @test betas[7].angular_momentum == 2
        @test betas[7].cutoff_radius_index == 661
        @test betas[7].beta[5] == -3.36344971742E-13

        @test size(nonlocal.dij) == (7, 7)
        @test nonlocal.dij[1, 1] == 2.22757816326E+01
        @test nonlocal.dij[1, 2] == -5.12800750666E+00
        @test nonlocal.dij[1, 3] == 1.22489807665E+01
        @test nonlocal.dij[2, 2] == -5.11535718643E-01
        @test nonlocal.dij[2, 3] == -4.01674388929E+00
        @test nonlocal.dij[3, 3] == -1.72050231320E+00
        @test nonlocal.dij[4, 4] == 2.31491498927E+01
        @test nonlocal.dij[4, 5] == -3.23205785707E+01
        @test nonlocal.dij[5, 5] == 4.32591292604E+01
        @test nonlocal.dij[6, 6] == 5.06898497757E+00
        @test nonlocal.dij[6, 7] == -3.89652282790E+00
        @test nonlocal.dij[7, 7] == 3.00051462460E+00
        @test count(d -> !iszero(d), nonlocal.dij) == 17

        augmentation = nonlocal.augmentation
        @test isnothing(augmentation.multipoles)
        @test !augmentation.q_with_l
        @test augmentation.nqf == 8

        @test size(augmentation.q) == (7, 7)
        @test augmentation.q[1, 1] == 1.79874106443E+00
        @test augmentation.q[1, 2] == augmentation.q[2, 1] == -1.29615056877E+00
        @test augmentation.q[7, 7] == 9.57081971594E-02

        @test length(augmentation.rinner) == 5
        @test all(r -> r == 1.20000000000E+00, augmentation.rinner)

        @test length(augmentation.qijs) == 28
        for qij in augmentation.qijs
            @test length(qij.qij) == file.header.mesh_size
        end
        @test length(augmentation.qfcoefs) == 28

        @test augmentation.qijs[1].first_index == 1
        @test augmentation.qijs[1].second_index == 1
        @test augmentation.qijs[1].qij[5] == 2.32704360102E-14
        @test augmentation.qfcoefs[1].qfcoef[5] == 3.43676888927E+02

        @test augmentation.qijs[2].first_index == 1
        @test augmentation.qijs[2].second_index == 2
        @test augmentation.qijs[2].qij[5] == -1.08277556385E-14
        @test augmentation.qfcoefs[2].qfcoef[5] == -1.10516062036E+02

        @test augmentation.qijs[28].first_index == 7
        @test augmentation.qijs[28].second_index == 7
        @test augmentation.qijs[28].qij[5] == -4.33337772184E-16
        @test augmentation.qfcoefs[28].qfcoef[5] == -5.52962873366E+01

        @test isnothing(augmentation.qijls)
        @test isnothing(augmentation.nqlc)
        @test isnothing(augmentation.shape)
        @test isnothing(augmentation.iraug)
        @test isnothing(augmentation.raug)
        @test isnothing(augmentation.l_max_aug)
        @test isnothing(augmentation.augmentation_epsilon)
        @test isnothing(augmentation.cutoff_r)
        @test isnothing(augmentation.cutoff_r_index)

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
            @test isnothing(chi.n)
            @test isnothing(chi.pseudo_energy)
            @test isnothing(chi.cutoff_radius)
            @test isnothing(chi.ultrasoft_cutoff_radius)
        end

        @test pswfc[1].label == "2S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.00
        @test pswfc[1].chi[5] == 4.74019527544E-08

        @test pswfc[3].label == "3S"
        @test pswfc[3].l == 0
        @test pswfc[3].occupation == 1.70
        @test pswfc[3].chi[5] == 1.02365544333E-08

        @test isnothing(file.full_wfc)

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 3.70544359504E-14

        @test isnothing(file.spin_orb)
        @test isnothing(file.paw)
        @test isnothing(file.gipaw)
    end
end
