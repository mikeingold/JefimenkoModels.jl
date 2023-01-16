################################################################################
#                        JEFIMENKOMODEL SOURCES
################################################################################

abstract type AbstractJefimenkoSource{T} end

    ############################################################################
    #                        VOLUME SOURCES
    ############################################################################

    abstract type AbstractVolumeSource{T} <: AbstractJefimenkoSource{T} end

        struct VolumeSource_Rectangular{T} <: AbstractVolumeSource{T}
            xlims::Tuple{Unitful.Length, Unitful.Length}
            ylims::Tuple{Unitful.Length, Unitful.Length}
            zlims::Tuple{Unitful.Length, Unitful.Length}
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
        end

        struct VolumeSource_Cylinder{T} <: AbstractVolumeSource{T}
            r::Tuple{Unitful.Length, Unitful.Length}
            ϕlims::Tuple{Unitful.Length, Unitful.Length}
            zlims::Tuple{Unitful.Length, Unitful.Length}
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
        end

        struct VolumeSource_Sphere{T} <: AbstractVolumeSource{T}
            r::Tuple{Unitful.Length, Unitful.Length}
            θlims::Tuple{Unitful.Length, Unitful.Length}
            ϕlims::Tuple{Unitful.Length, Unitful.Length}
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
        end
    
    export VolumeSource_Rectangular, VolumeSource_Cylinder, VolumeSource_Sphere

    ############################################################################
    #                        SURFACE SOURCES
    ############################################################################

    abstract type AbstractSurfaceSource{T} <: AbstractJefimenkoSource{T} end

        struct SurfaceSource_Rectangle{T} <: AbstractSurfaceSource{T}
            xlims::Tuple{Unitful.Length, Unitful.Length}
            ylims::Tuple{Unitful.Length, Unitful.Length}
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
        end

        struct SurfaceSource_Disk{T} <: AbstractSurfaceSource{T}
            ρ₀::Unitful.Length
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
        end
        
    export SurfaceSource_Rectangle, SurfaceSource_Disk

    ############################################################################
    #                        LINE SOURCES
    ############################################################################

    abstract type AbstractLineSource{T} <: AbstractJefimenkoSource{T} end

        struct LineSource_Straight{T} <: AbstractLineSource{T}
            ā::CoordinateCartesian
            b̄::CoordinateCartesian
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
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