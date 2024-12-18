################################################################################
#                        RADIATION SOURCES
################################################################################

struct RadiationSource{G, F1, F2, F3, F4} where {G <: Meshes.Geometry, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}
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

################################################################################
#                        JEFIMENKO MODELS
################################################################################

struct JefimenkoModel{M, S} where {M <: AbstractPropagationMedia, S <: RadiationSource}
    media::M
    sources::Vector{S}
    metadata::Dict{Symbol,Any}
end
