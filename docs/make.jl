push!(LOAD_PATH, "../src/")

using Documenter, JefimenkoModels

makedocs(
    sitename="JefimenkoModels.jl",
    authors = "Michael Ingold",
    format = Documenter.HTML(analytics = "G-6SFXR709Z9",
			     canonical = "https://mikeingold.github.io/JefimenkoModels.jl/stable/"),
    pages = [
            "Home" => "index.md",
            "Tutorial" => "tutorial.md",
            "Reference" => "reference.md"
	    ]
	)

deploydocs(repo="github.com/mikeingold/JefimenkoModels.jl.git")

