################################################################################
#                        PROPAGATION MEDIA
################################################################################

abstract type AbstractPropagationMedia end

struct SimpleMedia{E <: Quantity, U <: Quantity, C <: Quantity} <: AbstractPropagationMedia
    epsilon::E
    mu::U
    c::C
end

#=
struct PropagationMedia_DiagonallyAnisotropic{E <: Quantity, U <: Quantity, C <: Quantity} <: AbstractPropagationMedia
    epsilon::Diagonal{E}
    mu::Diagonal{U}
    c::Diagonal{C}
end
=#

function Base.getproperty(media::AbstractPropagationMedia, sym::Symbol)
    if sym in (:epsilon, :mu, :c)  # included
        return getfield(media, sym)
    elseif sym in (:ε, :ϵ)  # aliases
        return getfield(media, :epsilon)
    elseif sym in (:μ)  # aliases
        return getfield(media, :mu)
    else  # fallback
        return getfield(media, sym)
    end
end

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
    N = Meshes.paramdim(geometry)

    if ismissing(rho_e)
        rho_e = (r̄, t) -> 0.0 * u"C" / u"m"^N
    end

    if ismissing(rho_h)
        rho_h = (r̄, t) -> 0.0 * u"Wb" / u"m"^N
    end

    if ismissing(J_e)
        J_e = (r̄, t) -> 0.0 * u"A" / u"m"^(N-1)
    end

    if ismissing(J_h)
        J_h = (r̄, t) -> 0.0 * u"V" / u"m"^(N-1)
    end

    return RadiationSource(geometry, rho_e, rho_h, J_e, J_h)
end

function Base.getproperty(source::RadiationSource, sym::Symbol)
    if sym in (:geometry, :rho_e, :rho_h, :J_e, :J_h)  # included
        return getfield(source, sym)
    elseif sym in (:ρₑ, :ρe)  # aliases
        return getfield(source, :rho_e)
    elseif sym in (:ρₕ, :ρh)  # aliases
        return getfield(source, :rho_h)
    elseif sym in (:Jₑ, :Je)  # aliases
        return getfield(source, :J_e)
    elseif sym in (:Jₕ, :Jh)  # aliases
        return getfield(source, :J_h)
    else  # fallback
        return getfield(source, sym)
    end
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
