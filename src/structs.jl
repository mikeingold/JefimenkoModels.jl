################################################################################
#                        JEFIMENKOMODEL SOURCES
################################################################################

# Type T defines the data type used for calculation, typically <: AbstractFloat
struct RadiationSource{G, T} where {G <: Meshes.Geometry, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}
    geometry::G
    rho_e::F1
    rho_h::F2
    J_e::F3
    J_h::F4
end

################################################################################
#                        PROPAGATION MEDIA
################################################################################

abstract type AbstractPropagationMedia end

    struct PropagationMedia_Simple <: AbstractPropagationMedia
        epsilon::Quantity
        mu::Quantity
        c::Quantity
    end

    struct PropagationMedia_DiagonallyAnisotropic <: AbstractPropagationMedia
        epsilon::Diagonal{Quantity}
        mu::Diagonal{Quantity}
        c::Diagonal{Quantity}
    end

export PropagationMedia_Simple, PropagationMedia_DiagonallyAnisotropic

################################################################################
#                        JEFIMENKO MODELS
################################################################################

struct JefimenkoModel{T}
    media::AbstractPropagationMedia
    sources::Vector{AbstractJefimenkoSource{T}}
    metadata::Dict{Symbol,Any}
end

export JefimenkoModel