push!(LOAD_PATH, "../src/")

using Documenter, JefimenkoModels

makedocs(sitename="JefimenkoModels.jl",
	 authors = "Michael Ingold",
	 pages = [
		"Home" => "index.md",
		"Examples" => "examples.md",
		"Reference" => "reference.md"
	 ]
	)

deploydocs(repo="github.com/mikeingold/JefimenkoModels.jl")

