using OceanGreensFunctionMethods
using Documenter

DocMeta.setdocmeta!(OceanGreensFunctionMethods, :DocTestSetup, :(using OceanGreensFunctionMethods); recursive=true)

makedocs(;
    modules=[OceanGreensFunctionMethods],
    authors="G Jake Gebbie <ggebbie@whoi.edu>",
    sitename="OceanGreensFunctionMethods.jl",
    format=Documenter.HTML(;
        canonical="https://ggebbie.github.io/OceanGreensFunctionMethods.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/OceanGreensFunctionMethods.jl",
    devbranch="main",
)
