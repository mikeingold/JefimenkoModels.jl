push!(LOAD_PATH, "../src/")

using Documenter
#using JefimenkoModels

makedocs(
	 sitename="JefimenkoModels.jl",
	 pages = [
		"Home" => "index.md",
		"Examples" => "examples.md",
		"Reference" => "reference.md"
	 ]
	)

deploydocs(repo="github.com/mikeingold/JefimenkoModels.jl")

