module JefimenkoModels
    import ForwardDiff
    import LinearAlgebra
    import Meshes
    import MeshIntegrals

    using PhysicalConstants.CODATA2018: c_0, ε_0, μ_0
    using StaticArrays
    using Unitful
    using Unitful.DefaultSymbols: W, A, V, C, m, s, rad

    include("structs.jl")
    export AbstractPropagationMedia, SimpleMedia
    export RadiationSource, JefimenkoModel

    include("utils.jl")
    export CLASSICAL_VACUUM, NULL_CHARGE, NULL_CURRENT, t′

    include("fields.jl")
    export E, H, P
end
