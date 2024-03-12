using GRASS
using Documenter

DocMeta.setdocmeta!(GRASS, :DocTestSetup, :(using GRASS); recursive=true)

Introduction = "Introduction" => "index.md"

Getting_Started = "Getting Started" => ["installation.md", "examples/basic.md", "examples/advanced.md"]

Input_Data = "Input Data" => ["input_data/description.md", "input_data/manipulation.md"]

License = "License" => "license.md"

Index = "Index" => "longlist.md"

pages = [Introduction, Getting_Started, License, Index]

makedocs(;
    modules=[GRASS],
    authors="Michael Palumbo",
    # repo="https://github.com/palumbom/GRASS/blob/{commit}{path}#{line}",
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
)
