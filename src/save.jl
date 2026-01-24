"""
Save a [`PsPFile`](@ref) struct into a file. The file type depends on the
type of the file to save (e.g. a `Psp8File` file is saved to a `.psp8` or
a `UpfFile` to a `.upf`). An important keyword argument is `format`, which
enables to request a particular version of the file format to be employed.
"""
function save_psp_file(path::AbstractString, psp::PsPFile, args...; kwargs...)
    open(path, "w") do io
        save_psp(io, psp, args...; kwargs...)
    end
end
