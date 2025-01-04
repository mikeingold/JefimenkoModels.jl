module JefimenkoModels
    import LinearAlgebra
    import Meshes
    import MeshIntegrals
    import PhysicalConstants
    
    import UnitfulChainRules
    import Zygote

    using StaticArrays
    using Unitful
    using Unitful.DefaultSymbols: A, C, V, W, Wb, m, s, rad

    include("structs.jl")
    export AbstractPropagationMedia, SimpleMedia, CLASSICAL_VACUUM
    export RadiationSource, JefimenkoModel

    include("fields.jl")
    export E, H, P
end
