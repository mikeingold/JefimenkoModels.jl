################################################################################
#                        PROPAGATION MEDIA
################################################################################

abstract type AbstractPropagationMedia end

struct SimpleMedia{E <: Quantity, U <: Quantity, C <: Quantity} <: AbstractPropagationMedia
    permittivity::E
    permeability::U
    c::C
end

function SimpleMedia(permittivity::E, permeability::U) where {E <: Quantity, U <: Quantity}
    ε = Unitful.uconvert((A * s) / (V * m), permittivity)
    μ = Unitful.uconvert((V * s) / (A * m), permeability)
    c = Unitful.uconvert(m / s, sqrt(1 / (ε * μ)))
    return SimpleMedia(ε, μ, c)
end

CLASSICAL_VACUUM = let
    ε₀ = Unitful.uconvert((A * s) / (V * m), float(PhysicalConstants.CODATA2018.ε_0))
    μ₀ = Unitful.uconvert((V * s) / (A * m), float(PhysicalConstants.CODATA2018.μ_0))
    c₀ = Unitful.uconvert(m / s, float(PhysicalConstants.CODATA2018.c_0))
    SimpleMedia(ε₀, μ₀, c₀)
end

#=
struct PropagationMedia_DiagonallyAnisotropic{E <: Quantity, U <: Quantity, C <: Quantity} <: AbstractPropagationMedia
    epsilon::Diagonal{E}
    mu::Diagonal{U}
    c::Diagonal{C}
end
=#

################################################################################
#                        RADIATION SOURCES
################################################################################

struct RadiationSource{G <: Meshes.Geometry, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}
    geometry::G
    rho_e::F1
    rho_h::F2
    J_e::F3
    J_h::F4
end

# Allow keyword construction where unspecified==NULL
function RadiationSource(
    geometry::Meshes.Geometry;
    rho_e = missing,
    rho_h = missing,
    J_e = missing,
    J_h = missing
)
    # Number of parametric dimensions in geometry is used to determine proper density units
    N = Meshes.paramdim(geometry)

    # Replace missing sources with null functions
    if ismissing(rho_e)
        rho_e = (r̄, t) -> 0.0 * (u"C" / u"m"^N)
    end
    if ismissing(rho_h)
        rho_h = (r̄, t) -> 0.0 * (u"Wb" / u"m"^N)
    end
    if ismissing(J_e)
        J_e = (r̄, t) -> 0.0 * (u"A" / u"m"^(N-1))
    end
    if ismissing(J_h)
        J_h = (r̄, t) -> 0.0 * (u"V" / u"m"^(N-1))
    end

    return RadiationSource(geometry, rho_e, rho_h, J_e, J_h)
end

################################################################################
#                        JEFIMENKO MODELS
################################################################################

struct JefimenkoModel{M <: AbstractPropagationMedia, S <: RadiationSource}
    media::M
    sources::Vector{S}
    metadata::Dict{Symbol, Any}

    # Metadata is optional
    function JefimenkoModel{M, S}(
        media::M,
        sources::Vector{S},
        metadata::Dict{Symbol, Any} = Dict{Symbol, Any}()
    ) where {M <: AbstractPropagationMedia, S <: RadiationSource}
        return new(media, sources, metadata)
    end
end

# Allow construction with single source
function JefimenkoModel(
    media::M,
    source::S,
    metadata::Dict{Symbol, Any} = Dict{Symbol, Any}()
) where {M <: AbstractPropagationMedia, S <: RadiationSource}
    return JefimenkoModel{M, S}(media, [source], metadata)
end
