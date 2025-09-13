_FILE_EXT_LOADERS = Dict(".upf" => UpfFile,
                         ".psp8" => Psp8File,
                         ".hgh" => HghFile)

"""
Parse a pseudopotential file into a `PsPFile` struct.
"""
function load_psp_file(path::AbstractString)
    _, ext = splitext(path)
    ext = lowercase(ext)
    ext in keys(_FILE_EXT_LOADERS) && return _FILE_EXT_LOADERS[ext](path)
    return error("Unsupported PsP file extension $(ext)")
end

function load_psp_file(family_name_or_dir::AbstractString, filename::AbstractString)
    dir = _resolve_family(family_name_or_dir)
    path = joinpath(dir, filename)
    isfile(path) && return load_psp_file(path)
    return error("PsP $filename does not exist")
end

"""
List all known pseudopotential families with additional information about the number and
format(s) of the pseudopotentials within if requested (slow!).
"""
function list_families(;with_info=false)
    artifacts_toml_path = Artifacts.find_artifacts_toml(pathof(PseudoPotentialIO))
    artifacts = Artifacts.load_artifacts_toml(artifacts_toml_path)
    if with_info
        families = Dict("name" => [], "n_psp" => [], "format" => [])
        for (name, artifact) in artifacts
            if artifact_exists(Base.SHA1(artifact["git-tree-sha1"]))
                family = load_family.(name)
                n_psp = length(family["Filename"])
                format = unique(family["Format"])
                push!(families["name"], name)
                push!(families["n_psp"], n_psp)
                push!(families["format"], format)
            else
                push!(families["name"], name)
                push!(families["n_psp"], missing)
                push!(families["format"], missing)
            end
        end
        header = ["Name", "No. PsP", "Format(s)"]
        data = [families["name"];; families["n_psp"];; families["format"]]
        @ptconf crop=:none
        @pt :header=header data
        @ptconfclean
    else
        families = collect(keys(artifacts))
        @ptconf crop=:none
        @pt :header=["Name"] families
        @ptconfclean
    end
    return families
end

"""
Show the elements contained in a pseudopotential family as a periodic table.
"""
function show_family_periodic_table(family_name_or_dir::AbstractString)
    family = load_family(family_name_or_dir)
    _show_family_header(family, family_name_or_dir)
    show_family_periodic_table(family)
    return family
end
function show_family_periodic_table(family)
    elements = [PeriodicTable.elements[Symbol(element)] for element in family["Element"]]
    elements_struct = PeriodicTable.Elements(elements)

    println()
    show(Core.stdout, MIME("text/plain"), elements_struct)
    return nothing
end

"""
List the pseudopotentials in a pseudopotential family in a pretty table.
The elements for which pseudos are shown can be restricted by passing a list of strings,
e.g. ["Ag"].
"""
function show_family_list(family_name_or_dir::AbstractString; elements=[])
    family = load_family(family_name_or_dir)
    _show_family_header(family, family_name_or_dir)
    show_family_list(family; elements)
    return family
end
function show_family_list(family; elements=[])
    _keys = ["Element", "Filename", "Valence Charge", "NLCC", "Spin Orbit", "Format",
             "Formalism", "File"]
    if !isempty(elements)
        indices = filter(i -> family["Element"][i] in elements, eachindex(family["Element"]))
        family = Dict(k => v[indices] for (k, v) in family)
    end
    hlines = [0, 1]
    for i in 2:length(family["Element"])
        if family["Element"][i - 1] != family["Element"][i]
            push!(hlines, i)
        end
    end
    push!(hlines, length(family["Element"]) + 1)
    data = [family["Element"];; family["Filename"];; family["Valence Charge"];;
            family["NLCC"];; family["Spin Orbit"];; family["Format"];; family["Formalism"];;
            family["File"]]
    @ptconf vlines=:all hlines=hlines
    @pt :header=_keys data
    @ptconfclean
    return nothing
end

function _show_family_header(family, family_name_or_dir)
        @printf "%10s: %s\n" "Family" family_name_or_dir
        @printf "%10s: %s\n" "Directory" _resolve_family(family_name_or_dir)
        @printf "%10s: %s\n" "Formalism" join(string.(unique(family["Formalism"])), ", ")
        @printf "%10s: %d\n" "Pseudos" length(family["Filename"])
        return nothing
end

"""
Load all pseudopotentials from a given family.
"""
function load_family(family_name_or_dir::AbstractString)
    dir = _resolve_family(family_name_or_dir)
    (_, _, filenames) = first(walkdir(dir))
    psp_filenames = filter(f -> !startswith(f, "."), filenames)
    psp_filenames = filter(f -> lowercase(splitext(f)[2]) in keys(_FILE_EXT_LOADERS),
                           psp_filenames)
    psp_paths = map(f -> joinpath(dir, f), psp_filenames)
    psp_files = load_psp_file.(psp_paths)
    data = Dict(
        "Element" => element.(psp_files),
        "Filename" => psp_filenames,
        "Valence Charge" => valence_charge.(psp_files),
        "NLCC" => has_model_core_charge_density.(psp_files),
        "Spin Orbit" => has_spin_orbit.(psp_files),
        "Format" => typeof.(psp_files),
        "Formalism" => formalism.(psp_files),
        "File" => psp_files
    )
    return data
end

# Resolve the name of a PsP family artifact into its directory
# or return the input if it's already a directory
function _resolve_family(family_name_or_dir::AbstractString)
    isdir(family_name_or_dir) && return family_name_or_dir
    try
        return @artifact_str(family_name_or_dir)
    catch
        error("PsP family $family does not exist")
    end
end
