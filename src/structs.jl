################################################################################
#                        JEFIMENKOMODEL SOURCES
################################################################################

# Type T defines the data type used for calculation, typically <: AbstractFloat
abstract type AbstractJefimenkoSource{T} end

    ############################################################################
    #                        VOLUME SOURCES
    ############################################################################

    abstract type AbstractVolumeSource{T} <: AbstractJefimenkoSource{T} end

        struct VolumeSource_Rectangular{T,pe,ph,Je,Jh} <: AbstractVolumeSource{T}
            xlims::Tuple{Unitful.Length, Unitful.Length}
            ylims::Tuple{Unitful.Length, Unitful.Length}
            zlims::Tuple{Unitful.Length, Unitful.Length}
            ρ::pe
            ρₕ::ph
            Jₑ::Je
            Jₕ::Jh
        end

        struct VolumeSource_Cylinder{T,pe,ph,Je,Jh} <: AbstractVolumeSource{T}
            r::Tuple{Unitful.Length, Unitful.Length}
            ϕlims::Tuple{Unitful.Length, Unitful.Length}
            zlims::Tuple{Unitful.Length, Unitful.Length}
            ρ::pe
            ρₕ::ph
            Jₑ::Je
            Jₕ::Jh
        end

        struct VolumeSource_Sphere{T,pe,ph,Je,Jh} <: AbstractVolumeSource{T}
            r::Tuple{Unitful.Length, Unitful.Length}
            θlims::Tuple{Unitful.Length, Unitful.Length}
            ϕlims::Tuple{Unitful.Length, Unitful.Length}
            ρ::pe
            ρₕ::ph
            Jₑ::Je
            Jₕ::Jh
        end
    
    export VolumeSource_Rectangular, VolumeSource_Cylinder, VolumeSource_Sphere

    ############################################################################
    #                        SURFACE SOURCES
    ############################################################################

    abstract type AbstractSurfaceSource{T} <: AbstractJefimenkoSource{T} end

        struct SurfaceSource_Rectangle{T,pe,ph,Je,Jh} <: AbstractSurfaceSource{T}
            xlims::Tuple{Unitful.Length, Unitful.Length}
            ylims::Tuple{Unitful.Length, Unitful.Length}
            ρ::pe
            ρₕ::ph
            Jₑ::Je
            Jₕ::Jh
        end

        struct SurfaceSource_Disk{T,pe,ph,Je,Jh} <: AbstractSurfaceSource{T}
            ρ₀::Unitful.Length
            ρ::pe
            ρₕ::ph
            Jₑ::Je
            Jₕ::Jh
        end
        
    export SurfaceSource_Rectangle, SurfaceSource_Disk

    ############################################################################
    #                        LINE SOURCES
    ############################################################################

    abstract type AbstractLineSource{T} <: AbstractJefimenkoSource{T} end

        struct LineSource_Straight{T,pe,ph,Je,Jh} <: AbstractLineSource{T}
            ā::AbstractCoordinate
            b̄::AbstractCoordinate
            ρ::pe
            ρₕ::ph
            Jₑ::Je
            Jₕ::Jh
        end

    export LineSource_Straight

################################################################################
#                        PROPAGATION MEDIA
################################################################################

abstract type AbstractPropagationMedia end

    struct PropagationMedia_Simple <: AbstractPropagationMedia
        ε::Quantity
        μ::Quantity
        c::Quantity
    end

    struct PropagationMedia_DiagonallyAnisotropic <: AbstractPropagationMedia
        ε::Diagonal{Quantity}
        μ::Diagonal{Quantity}
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
