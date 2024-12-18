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

function Base.getproperty(s::RadiationSource, sym::Symbol)
    if sym in (:geometry, :rho_e, :rho_h, :J_e, :J_h)  # included
        return getfield(s, sym)
    elseif sym in (:ρₑ, :ρe)  # aliases
        return getfield(s, :rho_e)
    elseif sym in (:ρₕ, :ρh)  # aliases
        return getfield(s, :rho_h)
    elseif sym in (:Jₑ, :Je)  # aliases
        return getfield(s, :J_e)
    elseif sym in (:Jₕ, :Jh)  # aliases
        return getfield(s, :J_h)
    else  # fallback
        return getfield(s, sym)
    end
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

function Base.getproperty(m::AbstractPropagationMedia, sym::Symbol)
    if sym in (:epsilon, :mu, :c)  # included
        return getfield(s, sym)
    elseif sym in (:ε, :ϵ)  # aliases
        return getfield(s, :epsilon)
    elseif sym in (:μ)  # aliases
        return getfield(s, :mu)
    else  # fallback
        return getfield(s, sym)
    end
end

################################################################################
#                        JEFIMENKO MODELS
################################################################################

struct JefimenkoModel{M, S} where {M <: AbstractPropagationMedia, S <: RadiationSource}
    media::M
    sources::Vector{S}
    metadata::Dict{Symbol,Any}
end
