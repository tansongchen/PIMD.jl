using PIMD
using Documenter

makedocs(;
    modules=[PIMD],
    authors="Songchen Tan",
    repo="https://github.com/tansongchen/PIMD.jl/blob/{commit}{path}#L{line}",
    sitename="PIMD.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tansongchen.github.io/PIMD.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tansongchen/PIMD.jl",
)
