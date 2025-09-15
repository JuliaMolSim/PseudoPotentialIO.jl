@testset "AbstractPsPFile" begin
    for filepath in TEST_FILEPATHS
        file = load_psp_file(filepath)
        @test isa(identifier(file), AbstractString)
        @test isa(format(file), AbstractString)
        @test isa(element(file), Element)
        @test all(isa.(functional(file), Functional))
        @test isa(valence_charge(file), Int)
        @test 0 <= valence_charge(file)
        @test isa(is_norm_conserving(file), Bool)
        @test isa(is_ultrasoft(file), Bool)
        @test isa(is_paw(file), Bool)
        @test (is_norm_conserving(file) || is_ultrasoft(file) || is_paw(file))  # at least one
        @test isa(has_spin_orbit(file), Bool)
        @test isa(has_model_core_charge_density(file), Bool)
        @test isa(formalism(file), AbstractString)
        @test formalism(file) in ("norm-conserving", "ultrasoft", "projector-augmented wave")
    end
end
