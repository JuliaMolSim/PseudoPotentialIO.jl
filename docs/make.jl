using PseudoPotentialIO
using Documenter
using Literate
import LibGit2
import Pkg

DEBUG = false  # Set to true to disable some checks and cleanup
SRCPATH   = joinpath(@__DIR__, "src")
BUILDPATH = joinpath(@__DIR__, "build")
ROOTPATH  = joinpath(@__DIR__, "..")
CONTINUOUS_INTEGRATION = get(ENV, "CI", nothing) == "true"
PPIOREV    = LibGit2.head(ROOTPATH)
PPIOBRANCH = try LibGit2.branch(LibGit2.GitRepo(ROOTPATH)) catch end
PPIOGH     = "github.com/azadoks/PseudoPotentialIO.jl"
PPIOREPO   = PPIOGH * ".git"
PAGES = ["Home" => "index.md",
         "Formats" => "formats.md",
         "api.md"]

# Setup julia dependencies for docs generation if not yet done
Pkg.activate(@__DIR__)
if !isfile(joinpath(@__DIR__, "Manifest.toml"))
    Pkg.develop(Pkg.PackageSpec(path=ROOTPATH))
    Pkg.instantiate()
end

# Get list of files from PAGES
extract_paths(pages::AbstractArray) = collect(Iterators.flatten(extract_paths.(pages)))
extract_paths(file::AbstractString) = [file]
extract_paths(pair::Pair) = extract_paths(pair.second)

# Transform files to *.md
transform_to_md(pages::AbstractArray) = transform_to_md.(pages)
transform_to_md(file::AbstractString) = first(splitext(file)) * ".md"
transform_to_md(pair::Pair) = (pair.first => transform_to_md(pair.second))

# Collect files to treat with Literate (i.e. the examples and the .jl files in the docs)
# The examples go to docs/literate_build/examples, the .jl files stay where they are
literate_files = map(filter!(endswith(".jl"), extract_paths(PAGES))) do file
       if startswith(file, "examples/")
           (src=joinpath(ROOTPATH, file), dest=joinpath(SRCPATH, "examples"), example=true)
       else
           (src=joinpath(SRCPATH, file), dest=joinpath(SRCPATH, dirname(file)), example=false)
       end
   end

# Run Literate on them all
for file in literate_files
       preprocess = file.example ? add_badges : identity
       Literate.markdown(file.src, file.dest;
                         flavor=Literate.DocumenterFlavor(),
                         credit=false, preprocess)
       Literate.notebook(file.src, file.dest; credit=false,
                         execute=CONTINUOUS_INTEGRATION || DEBUG)
   end

DocMeta.setdocmeta!(PseudoPotentialIO, :DocTestSetup, :(using PseudoPotentialIO);
                    recursive=true)

makedocs(;
         modules=[PseudoPotentialIO],
         authors="Austin Zadoks",
         repo="https://github.com/azadoks/PseudoPotentialIO.jl/blob/{commit}{path}#{line}",
         sitename="PseudoPotentialIO.jl",
         format=Documenter.HTML(;
                                prettyurls=CONTINUOUS_INTEGRATION,
                                canonical="https://azadoks.github.io/PseudoPotentialIO.jl",
                                edit_link="main",
                                assets=String[]),
         pages=transform_to_md(PAGES),
         checkdocs=:exports)

# Deploy docs to gh-pages branch
deploydocs(; repo=PPIOREPO, devbranch="main")

# Remove generated example files
if !DEBUG
    for file in literate_files
        base = splitext(basename(file.src))[1]
        for ext in [".ipynb", ".md"]
            rm(joinpath(file.dest, base * ext), force=true)
        end
    end
end

if !CONTINUOUS_INTEGRATION
    println("\nDocs generated, try $(joinpath(BUILDPATH, "index.html"))")
end
