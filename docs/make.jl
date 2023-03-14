push!(LOAD_PATH, "../src/")

using Documenter, JefimenkoModels

makedocs(sitename="JefimenkoModels.jl",
	 authors = "Michael Ingold",
	 pages = [
		"Home" => "index.md",
		"Tutorial" => "tutorial.md",
		"Reference" => "reference.md"
	 ]
	)

deploydocs(repo="github.com/mikeingold/JefimenkoModels.jl")

