@testset "upf_augmentation" begin
    # Most fixtures come from the existing test artifacts: pslibrary ultrasoft (Si) and
    # PAW (Al), plus a PAW file from a different generator (Dy) for l channels up to f.
    # The fully-relativistic case (PP_AEWFC_REL appears only in fully-relativistic PAW
    # files) is fetched on demand -- see the FR_PAW_FIXTURE machinery in fixtures.jl.

    @testset "ultrasoft (Si rrkjus)" begin
        file = load_psp_file(UPF2_CASE_FILEPATHS["Si.pbe-n-rrkjus_psl.1.0.0.UPF"])
        aug = upf_augmentation(file)
        @test aug isa UpfAugmentationData
        @test !aug.is_paw
        @test aug.aewfcs === nothing
        @test aug.pswfcs === nothing
        @test aug.ae_rel_wfcs === nothing
        @test aug.ae_core_density === nothing
        @test aug.ae_vloc === nothing
        # No augmentation-sphere truncation without PAW partial waves to compare.
        @test aug.i_cut == length(file.mesh.r)

        nproj = file.header.number_of_proj
        @test nproj == 6
        @test aug.proj_l == [0, 0, 1, 1, 2, 2]
        @test size(aug.q) == (nproj, nproj)
        @test size(aug.qq_at) == (nproj, nproj)
        @test issymmetric(aug.qq_at)
        # Upper triangle of a 6x6 index pairing.
        @test length(aug.qfuncl) == 21

        # Literal reference values pinning the reshaping and integration conventions.
        @test aug.qq_at[1, 1] ≈ -0.07715561497721883 rtol=1e-12
        @test aug.qq_at[1, 2] ≈ -0.06659968553281527 rtol=1e-12
        @test aug.qfuncl[(1, 1)][0][500] ≈ 0.0001652228511050223 rtol=1e-12

        # qq_at recomputed from qfuncl agrees with the file's own PP_Q at l=0.
        for i = 1:nproj, j = 1:nproj
            aug.proj_l[i] == aug.proj_l[j] || continue
            haskey(aug.qfuncl, (min(i, j), max(i, j))) || continue
            @test isapprox(aug.qq_at[i, j], aug.q[i, j]; atol=1e-4)
        end

        # An s-p pair carries only the l=1 channel by the triangle/parity rule.
        @test sort(collect(keys(aug.qfuncl[(1, 3)]))) == [1]

        # qfuncl keys respect the triangle/parity selection rule throughout.
        for ((i, j), channels) in aug.qfuncl
            li, lj = aug.proj_l[i], aug.proj_l[j]
            for l in keys(channels)
                @test abs(li - lj) <= l <= li + lj
                @test iseven(l - li - lj)
                @test length(channels[l]) == length(file.mesh.r)
            end
        end
    end

    @testset "PAW (Al kjpaw)" begin
        file = load_psp_file(UPF2_CASE_FILEPATHS["Al.pbe-n-kjpaw_psl.1.0.0.UPF"])
        aug = upf_augmentation(file)
        @test aug.is_paw
        nproj = file.header.number_of_proj
        @test nproj == 6
        @test aug.aewfcs !== nothing && length(aug.aewfcs) == nproj
        @test aug.pswfcs !== nothing && length(aug.pswfcs) == nproj
        @test aug.ae_rel_wfcs === nothing  # scalar-relativistic: no small components
        @test all(length(w) == length(file.mesh.r) for w in aug.aewfcs)
        @test all(length(w) == length(file.mesh.r) for w in aug.pswfcs)

        # Taken from the file's own PP_AUGMENTATION/cutoff_r_index.
        @test aug.i_cut == 839
        @test aug.qq_at[1, 1] ≈ -0.030709687725573528 rtol=1e-12
        @test aug.qfuncl[(1, 1)][0][500] ≈ 0.0003065559814337421 rtol=1e-12

        # AE and PS partial waves are identical beyond i_cut: the defining property of
        # the augmentation sphere, and the reason one-center sums truncate there.
        for beta in eachindex(aug.aewfcs)
            @test aug.aewfcs[beta][aug.i_cut:end] == aug.pswfcs[beta][aug.i_cut:end]
        end

        @test aug.ae_core_density !== nothing
        @test length(aug.ae_core_density) == length(file.mesh.r)
        @test aug.ae_core_density[1] ≈ 1490.764267098482 rtol=1e-12
        # ae_vloc keeps the file's raw Rydberg values, like PP_LOCAL.
        @test aug.ae_vloc !== nothing
        @test length(aug.ae_vloc) == length(file.mesh.r)
        @test aug.ae_vloc == file.paw.ae_vloc
        @test aug.ae_vloc[1] ≈ -370575.6139724471 rtol=1e-12
    end

    @testset "PAW from a different generator (Dy)" begin
        file = load_psp_file(UPF2_CASE_FILEPATHS["Dy.GGA-PBE-paw-v1.0.UPF"])
        aug = upf_augmentation(file)
        @test aug.is_paw
        nproj = file.header.number_of_proj
        @test nproj == 8
        @test aug.proj_l == [0, 0, 1, 1, 2, 2, 3, 3]
        @test length(aug.qfuncl) == 36
        @test aug.i_cut == 969
        @test aug.qq_at[1, 1] ≈ -0.1369412313387547 rtol=1e-12
        @test issymmetric(aug.qq_at)

        # f channels: the l range must reach 2*l_max = 6.
        @test maximum(l for channels in values(aug.qfuncl) for l in keys(channels)) == 6
        for ((i, j), channels) in aug.qfuncl
            li, lj = aug.proj_l[i], aug.proj_l[j]
            for l in keys(channels)
                @test abs(li - lj) <= l <= li + lj
                @test iseven(l - li - lj)
            end
        end
    end

    if !haskey(UPF2_CASE_FILEPATHS, FR_PAW_FIXTURE)
        @test_skip "fully-relativistic PAW (C rel-kjpaw): fixture download unavailable"
    else
        @testset "fully-relativistic PAW (C rel-kjpaw)" begin
            file = load_psp_file(UPF2_CASE_FILEPATHS[FR_PAW_FIXTURE])
            aug = upf_augmentation(file)
            @test aug.is_paw
            @test file.header.has_so
            nproj = file.header.number_of_proj
            @test nproj == 6
            @test aug.proj_l == [0, 0, 1, 1, 1, 1]

            # Fully-relativistic files concatenate PP_AEWFC (Dirac large components) and
            # PP_AEWFC_REL (small components) into a single 2*nproj block, which the
            # accessor splits back apart.
            @test length(file.full_wfc.aewfcs) == 2nproj
            @test length(aug.aewfcs) == nproj
            @test length(aug.ae_rel_wfcs) == nproj
            @test all(aug.aewfcs[i] != aug.ae_rel_wfcs[i] for i = 1:nproj)
            @test all(length(w) == length(file.mesh.r) for w in aug.ae_rel_wfcs)

            @test aug.i_cut == 749
            @test aug.qq_at[1, 1] ≈ -0.06118443141872321 rtol=1e-12
            @test aug.aewfcs[1][500] ≈ 0.2988187274128993 rtol=1e-12
            @test aug.ae_rel_wfcs[1][500] ≈ -0.007146824131769734 rtol=1e-12
        end
    end

    @testset "errors" begin
        # A norm-conserving file has no PP_AUGMENTATION block.
        ncfile = load_psp_file(UPF2_CASE_FILEPATHS["56_Ba_m.upf"])
        @test_throws ErrorException upf_augmentation(ncfile)

        # Legacy q_with_l = false ultrasoft files (UPF v1 qfcoef/rinner pseudization,
        # e.g. the GBRV family) are deliberately not decoded rather than mis-decoded.
        gbrv = load_psp_file(UPF1_CASE_FILEPATHS["ag_lda_v1.4.uspp.F.UPF"])
        @test PseudoPotentialIO.is_ultrasoft(gbrv)
        @test_throws ErrorException upf_augmentation(gbrv)
    end
end
