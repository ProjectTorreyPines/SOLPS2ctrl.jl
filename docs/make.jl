using Documenter
using SOLPS2ctrl

makedocs(;
    modules=[SOLPS2ctrl],
    format=Documenter.HTML(),
    sitename="SOLPS2ctrl",
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/ProjectTorreyPines/SOLPS2ctrl.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    versions=["stable" => "v^", "v#.#"],
)
