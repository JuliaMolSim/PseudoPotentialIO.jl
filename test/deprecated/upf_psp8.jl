@testset "PSP8--UPF" begin
    pseudo_pairs = [(upf=load_upf(pair[1]), psp=load_psp8(pair[2])) for pair in UPF2_PSP8_FILEPATHS]

    @testset "header" begin
        for pair in pseudo_pairs
            psphead = pair.psp["header"]
            upfhead = pair.upf["header"]
            @test psphead["z_atom"] == elements[Symbol(upfhead["element"])].number
            @test psphead["z_valence"] == upfhead["z_valence"]
            @test psphead["l_max"] == upfhead["l_max"]
            @test psphead["core_correction"] == upfhead["core_correction"]
            @test psphead["has_so"] == upfhead["has_so"]
            if !psphead["has_so"]
                @test sum(psphead["number_of_proj"]) == upfhead["number_of_proj"]
            end
        end
    end

    @testset "radial_grid" begin
        for pair in pseudo_pairs
            mesh_size = pair.psp["header"]["mesh_size"]
            @test isapprox(
                pair.psp["radial_grid"],
                pair.upf["radial_grid"][begin:mesh_size]
            )
        end
    end

    @testset "local_potential" begin
        for pair in pseudo_pairs
            mesh_size = pair.psp["header"]["mesh_size"]
            @test isapprox(
                pair.psp["local_potential"],
                pair.upf["local_potential"][begin:mesh_size] ./ 2,
            )
        end
    end

    @testset "nlcc" begin
        for pair in pseudo_pairs
            if pair.psp["header"]["core_correction"]
                mesh_size = pair.psp["header"]["mesh_size"]
                @test isapprox(
                    pair.psp["nlcc"]["core_charge_density"] ./ 4π,
                    pair.upf["core_charge_density"][begin:mesh_size],
                )
            end
        end
    end

    @testset "nonlocal_potential" begin
        non_so_pairs = filter(pair -> !pair.psp["header"]["has_so"], pseudo_pairs)
        # @testset "Vnl $(pair.psp["header"]["filename"])" for pair in non_so_pairs
        for pair in non_so_pairs
            mesh_size = pair.psp["header"]["mesh_size"]
            idx_upf = 1
            for l = 0:pair.psp["header"]["l_max"]
                for i = 1:length(pair.psp["beta_projectors"][l+1])
                    psp_β = pair.psp["beta_projectors"][l+1][i]
                    upf_β = pair.upf["beta_projectors"][idx_upf]["radial_function"]
                    ir_max = min(length(psp_β), length(upf_β))
                    
                    psp_ekb = pair.psp["ekb"][l+1][i,i]
                    upf_Dij = pair.upf["D_ion"][idx_upf,idx_upf]

                    psp_βekbβ = psp_β[1:ir_max] * psp_ekb * psp_β[1:ir_max]'
                    upf_βekbβ = upf_β[1:ir_max] * upf_Dij * upf_β[1:ir_max]'

                    @test psp_βekbβ ≈ upf_βekbβ ./ 2 atol=1e-4
                    idx_upf += 1
                end
            end
        end
    end
end
