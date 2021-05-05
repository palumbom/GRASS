using GRASS
using Documenter

DocMeta.setdocmeta!(GRASS, :DocTestSetup, :(using GRASS); recursive=true)

makedocs(;
    modules=[GRASS],
    authors="Michael Palumbo",
    repo="https://github.com/palumbom/GRASS.jl/blob/{commit}{path}#{line}",
    sitename="GRASS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://palumbom.github.io/GRASS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/palumbom/GRASS.jl",
)
