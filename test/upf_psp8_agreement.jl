function test_psp8_upf_agreement(psp8::Psp8File, upf::UpfFile)
    psp8head = psp8.header
    upfhead = upf.header

    @testset "header" begin
        @test psp8head.zatom == elements[Symbol(upfhead.element)].number
        @test psp8head.zion == upfhead.z_valence
        @test psp8head.lmax == upfhead.l_max
        @test (psp8head.fchrg > 0) == upfhead.core_correction
        @test !isnothing(psp8head.nprojso) == upfhead.has_so
        if !isnothing(psp8head.nprojso)
            @test sum(psp8head.nprojso) == upfhead.number_of_proj
        end
    end

    @testset "radial_grid" begin
        mesh_size = psp8head.mmax
        @test isapprox(
            psp8.rgrid,
            upf.mesh.r[begin:mesh_size]
        )
    end

    @testset "local_potential" begin
        mesh_size = psp8head.mmax
        @test isapprox(
            psp8.v_local,
            upf.local_[begin:mesh_size] ./ 2,
        )
    end

    @testset "nlcc" begin
        if upfhead.core_correction
            mesh_size = psp8head.mmax
            @test isapprox(
                psp8.rhoc ./ 4π,
                upf.nlcc[begin:mesh_size],
            )
        end
    end

    @testset "nonlocal_potential" begin
        if !upfhead.has_so
            mesh_size = psp8head.mmax
            idx_upf = 1
            for l = 0:psp8head.lmax
                for i = 1:length(psp8.projectors[l+1])
                    psp_β = psp8.projectors[l+1][i]
                    upf_β_data = upf.nonlocal.betas[idx_upf]
                    upf_β = upf_β_data.beta[1:upf_β_data.cutoff_radius_index]
                    ir_max = min(length(psp_β), length(upf_β))
                    
                    psp_ekb = psp8.ekb[l+1][i]
                    upf_Dij = upf.nonlocal.dij[idx_upf,idx_upf]

                    psp_βekbβ = psp_β[1:ir_max] * psp_ekb * psp_β[1:ir_max]'
                    upf_βekbβ = upf_β[1:ir_max] * upf_Dij * upf_β[1:ir_max]'

                    @test psp_βekbβ ≈ upf_βekbβ ./ 2 atol=1e-4
                    idx_upf += 1
                end
            end
        end
    end
end

@testset "PSP8--UPF agreement" begin
    for (element, upf_path, psp8_path) in UPF2_PSP8_FILEPATHS
        @testset "$element" begin
            upf = load_psp_file(upf_path)
            psp8 = load_psp_file(psp8_path)

            test_psp8_upf_agreement(psp8, upf)
        end
    end
end

@testset "PSP8 -> UPF conversion agreement" begin
    psp8 = load_psp_file(
        joinpath(_resolve_family("pd_nc_sr_pbe_standard_0.4.1_psp8"), "Pb.psp8"),
    )
    upf = UpfFile(psp8)

    test_psp8_upf_agreement(psp8, upf)
end
