_FILE_EXT_LOADERS = Dict(".upf" => UpfFile,
                         ".psp8" => Psp8File,
                         ".hgh" => HghFile)

"""
Parse a pseudopotential file into a [`PsPFile`](@ref) struct.
"""
function load_psp_file(path::AbstractString)
    _, ext = splitext(path)
    ext = lowercase(ext)
    ext in keys(_FILE_EXT_LOADERS) && return _FILE_EXT_LOADERS[ext](path)
    return error("Unsupported PsP file extension $(ext)")
end
