using PseudoPotentialIO
using Documenter

DocMeta.setdocmeta!(DftFunctionals, :DocTestSetup, :(using PseudoPotentialIO); recursive=true)

makedocs(;
    modules=[PseudoPotentialIO],
    authors="Austin Zadoks, Bruno Ploumhans, Michael F. Herbst and contributors",
    sitename="PseudoPotentialIO.jl",
    format=Documenter.HTML(;
        canonical="https://juliamolsim.github.io/PseudoPotentialIO.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "formats.md",
        "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/JuliaMolSim/PseudoPotentialIO.jl",
    devbranch="main",
)
