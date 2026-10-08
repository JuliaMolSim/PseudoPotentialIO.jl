# Write every UPF v2 norm-conserving and every psp8 fixture, read it back, and compare every
# field of the two structs (floating-point data to the precision of the writers).

# Paths of the fields that differ between `a` and `b`, recursively through structs and arrays
function field_mismatches(a, b; rtol=1e-12)
    mismatches = String[]
    _collect_mismatches!(mismatches, a, b, "file", rtol)
    return mismatches
end

function _collect_mismatches!(mismatches, a, b, path, rtol)
    if a isa AbstractFloat && b isa AbstractFloat
        isapprox(a, b; rtol) || push!(mismatches, path)
    elseif a isa AbstractArray{<:AbstractFloat} && b isa AbstractArray{<:AbstractFloat}
        size(a) == size(b) && isapprox(a, b; rtol) || push!(mismatches, path)
    elseif a isa AbstractArray && b isa AbstractArray
        if size(a) != size(b)
            push!(mismatches, path)
        else
            for i in eachindex(a, b)
                _collect_mismatches!(mismatches, a[i], b[i], "$path[$i]", rtol)
            end
        end
    elseif isstructtype(typeof(a)) && typeof(a) == typeof(b) && !(a isa AbstractString) &&
           fieldcount(typeof(a)) > 0
        for n in fieldnames(typeof(a))
            _collect_mismatches!(mismatches, getfield(a, n), getfield(b, n), "$path.$n", rtol)
        end
    else
        isequal(a, b) || push!(mismatches, path)
    end
    return mismatches
end

function roundtrip(file, ext)
    outpath = joinpath(mktempdir(; cleanup=true), "roundtrip" * ext)
    save_psp_file(outpath, file)
    return load_psp_file(outpath)
end

# the identifier is the file name
without_identifier(file::UpfFile) = [getfield(file, n) for n in fieldnames(UpfFile)
                                     if n != :identifier]
without_identifier(file::Psp8File) = [getfield(file, n) for n in fieldnames(Psp8File)
                                      if n != :identifier]

@testset "Round trip" begin
    @testset "UPF v2.0.1 NC" begin
        for filepath in values(UPF2_CASE_FILEPATHS)
            file = load_psp_file(filepath)
            is_norm_conserving(file) || continue
            newfile = roundtrip(file, ".upf")
            mismatches = field_mismatches(without_identifier(file),
                                          without_identifier(newfile))
            @test isempty(mismatches) || (filepath, mismatches)
        end
    end

    @testset "PSP8" begin
        for filepath in values(PSP8_FILEPATHS)
            file = load_psp_file(filepath)
            newfile = roundtrip(file, ".psp8")
            mismatches = field_mismatches(without_identifier(file),
                                          without_identifier(newfile))
            @test isempty(mismatches) || (filepath, mismatches)
        end
    end
end

@testset "Libxc functional strings" begin
    @test PseudoPotentialIO.libxc_to_qe_libxc("LDA_X LDA_C_PW") == "XC-001L-012L"
    @test PseudoPotentialIO.libxc_to_qe_libxc("GGA_X_PBE GGA_C_PBE") ==
          "XC-000I-000I-101L-130L"
    @test PseudoPotentialIO.libxc_to_qe_libxc("GGA_C_PBE GGA_X_PBE") ==
          "XC-000I-000I-101L-130L"
    @test PseudoPotentialIO.libxc_to_qe_libxc("MGGA_X_R2SCAN MGGA_C_R2SCAN") ==
          "XC-000I-000I-000I-000I-497L-498L"
end
