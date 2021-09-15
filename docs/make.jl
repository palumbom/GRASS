using GRASS
using Documenter

DocMeta.setdocmeta!(GRASS, :DocTestSetup, :(using GRASS); recursive=true)

Home = "Home" => ["index.md", "installation.md"]

Examples = "Examples" => ["examples/spectra.md"]

License = "License" => "license.md"

Index = "Index" => longlist.md

pages = [Home, Examples, License, Index]

makedocs(;
    modules=[GRASS],
    authors="Michael Palumbo",
    repo="https://github.com/palumbom/GRASS/blob/{commit}{path}#{line}",
    sitename="GRASS",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://palumbom.github.io/GRASS",
        assets=String[],
    ),
    pages=pages,
)

deploydocs(;
    repo="github.com/palumbom/GRASS",
    devbranch = "main",
)
