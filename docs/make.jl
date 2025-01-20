push!(LOAD_PATH, "../src/")

using Documenter, ReefBiodiversityAccountSetup

makedocs(;
    sitename="ReefBiodiversityAccountSetup Documentation", pages=[
        "index.md", "API.md"]
)

deploydocs(;
    repo="github.com/open-AIMS/ReefBiodiversityAccountSetup.jl.git")
