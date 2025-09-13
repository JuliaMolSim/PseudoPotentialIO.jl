@testset "UPF" begin
    @testset "Internal data consistency" begin
        @testset "$filepath" for filepath in UPF_FILEPATHS
            file = load_psp_file(filepath)

            @test haskey(PeriodicTable.elements, Symbol(file.header.element))
            z_atom = PeriodicTable.elements[Symbol(titlecase(file.header.element))].number
            @test file.header.z_valence <= z_atom

            if file.header.pseudo_type == "PAW"
                @test file.header.is_paw
                @test file.header.is_ultrasoft
                @test !file.header.is_coulomb
            elseif file.header.pseudo_type in ("US", "USPP")
                @test file.header.is_ultrasoft
                @test !any((file.header.is_paw, file.header.is_coulomb))
            elseif file.header.pseudo_type == "1/r"
                @test file.header.is_coulomb
                @test !any((file.header.is_paw, file.header.is_ultrasoft))
            else
                @test file.header.pseudo_type in ("NC", "SL")
            end

            if file.header.is_ultrasoft | file.header.is_paw
                @test !isnothing(file.nonlocal.augmentation)
                augmentation = file.nonlocal.augmentation

                @test size(augmentation.q) == (file.header.number_of_proj,
                                               file.header.number_of_proj)

                if augmentation.q_with_l
                    @test isnothing(augmentation.qijs)
                    @test !isnothing(augmentation.qijls)
                    nqijl = 0
                    for i in 1:file.header.number_of_proj, j in i:file.header.number_of_proj
                        li = file.nonlocal.betas[i].angular_momentum
                        lj = file.nonlocal.betas[j].angular_momentum
                        for l in abs(li - lj):2:(li + lj)
                            nqijl += 1
                        end
                    end
                    @test length(augmentation.qijls) == nqijl
                else
                    @test !isnothing(augmentation.qijs)
                    @test isnothing(augmentation.qijls)
                    @test (length(augmentation.qijs) ==
                           file.header.number_of_proj * (file.header.number_of_proj + 1) / 2)
                end
            else
                @test isnothing(file.nonlocal.augmentation)
            end

            @test !isnothing(file.local_)
            @test length(file.local_) == file.header.mesh_size

            if file.header.core_correction
                @test !isnothing(file.nlcc)
                @test length(file.nlcc) == file.header.mesh_size
            end

            if file.header.has_so
                @test !isnothing(file.spin_orb)
                @test length(file.spin_orb.relbetas) == file.header.number_of_proj
                @test length(file.spin_orb.relwfcs) == file.header.number_of_wfc
            end

            if file.header.has_gipaw
                @test !isnothing(file.gipaw)
                #TODO check the contents
            else
                @test isnothing(file.gipaw)
            end

            if file.header.number_of_wfc > 0
                @test length(file.pswfc) == file.header.number_of_wfc
                for chi in file.pswfc
                    @test length(chi.chi) == file.header.mesh_size
                end
            else
                @test isnothing(file.pswfc)
            end

            if file.header.has_wfc
                @test !isnothing(file.full_wfc)
                #TODO check contents
            else
                @test isnothing(file.full_wfc)
            end

            @test (length(file.mesh.r) == length(file.mesh.rab) == file.header.mesh_size)
            if !isnothing(file.mesh.mesh)
                @test file.mesh.mesh == file.header.mesh_size
            end

            @test length(file.nonlocal.betas) == file.header.number_of_proj
            for beta in file.nonlocal.betas
                @test (length(beta.beta) == file.header.mesh_size) | (length(beta.beta) == beta.cutoff_radius_index)
            end
        end
    end
end
