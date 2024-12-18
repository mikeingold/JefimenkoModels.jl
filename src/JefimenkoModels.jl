module JefimenkoModels
    import Meshes

    using Unitful, UnitfulCoordinateSystems
    using Unitful.DefaultSymbols: W, A, V, C, m, s, rad
    using PhysicalConstants.CODATA2018: c_0, ε_0, μ_0
    using ForwardDiff, Integrals, LinearAlgebra, StaticArrays

    include("structs.jl")
    export RadiationSource, JefimenkoModel, PropagationMedia_Simple, PropagationMedia_DiagonallyAnisotropic

    include("accessors.jl")

    include("utils.jl")
    export CLASSICAL_VACUUM, NULL_CHARGE, NULL_CURRENT, t′

    include("integrands.jl")

    include("fields.jl")
    export E, H, P
end
