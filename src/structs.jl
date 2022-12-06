###########################################################################
#                        JEFIMENKOMODEL SOURCES
###########################################################################

abstract type JefimenkoSource{T} end

    #######################################################################
    #                        VOLUME SOURCES
    #######################################################################

    abstract type VolumeSource{T} <: JefimenkoSource{T} end

        struct VolumeSource_Rectangular{T} <: VolumeSource{T}
            xlims::Tuple{Unitful.Length, Unitful.Length}
            ylims::Tuple{Unitful.Length, Unitful.Length}
            zlims::Tuple{Unitful.Length, Unitful.Length}
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
        end

        struct VolumeSource_Cylinder{T} <: VolumeSource{T}
            r::Tuple{Unitful.Length, Unitful.Length}
            ϕlims::Tuple{Unitful.Length, Unitful.Length}
            zlims::Tuple{Unitful.Length, Unitful.Length}
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
        end

        struct VolumeSource_Sphere{T} <: VolumeSource{T}
            r::Tuple{Unitful.Length, Unitful.Length}
            θlims::Tuple{Unitful.Length, Unitful.Length}
            ϕlims::Tuple{Unitful.Length, Unitful.Length}
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
        end
    
    export VolumeSource_Rectangular, VolumeSource_Cylinder, VolumeSource_Sphere

    #######################################################################
    #                        SURFACE SOURCES
    #######################################################################

    abstract type SurfaceSource{T} <: JefimenkoSource{T} end

        struct SurfaceSource_Rectangle{T} <: SurfaceSource{T}
            xlims::Tuple{Unitful.Length, Unitful.Length}
            ylims::Tuple{Unitful.Length, Unitful.Length}
            ρₑ::Function
            ρₕ::Function
            Jₑ::Function
            Jₕ::Function
        end

        abstract type SurfaceSource_Disk{T} <: SurfaceSource{T} end

            struct SurfaceSource_Disk_General{T} <: SurfaceSource_Disk{T}
                ρ₀::Unitful.Length
                ρₑ::Function
                ρₕ::Function
                Jₑ::Function
                Jₕ::Function
            end

            struct SurfaceSource_Disk_CurrentsOnly{T} <: SurfaceSource_Disk{T}
                ρ₀::Unitful.Length
                Jₑ::Function
                Jₕ::Function
            end

            struct SurfaceSource_Disk_ElectricOnly{T} <: SurfaceSource_Disk{T}
                ρ₀::Unitful.Length
                ρₑ::Function
                Jₑ::Function
            end
        
    export SurfaceSource_Rectangle
    export SurfaceSource_Disk_General,  SurfaceSource_Disk_CurrentsOnly,  SurfaceSource_Disk_ElectricOnly

    #######################################################################
    #                        LINE SOURCES
    #######################################################################

    abstract type LineSource{T} <: JefimenkoSource{T} end

        abstract type LineSource_Straight{T} <: LineSource{T} end

            struct LineSource_Straight_General{T} <: LineSource_Straight{T}
                ā::CoordinateCartesian
                b̄::CoordinateCartesian
                ρₑ::Function
                ρₕ::Function
                Jₑ::Function
                Jₕ::Function
            end

            struct LineSource_Straight_ElectricOnly{T} <: LineSource_Straight{T}
                ā::CoordinateCartesian
                b̄::CoordinateCartesian
                ρₑ::Function
                Jₑ::Function
            end

            struct LineSource_Straight_CurrentsOnly{T} <: LineSource_Straight{T}
                ā::CoordinateCartesian
                b̄::CoordinateCartesian
                Jₑ::Function
                Jₕ::Function
            end

    export LineSource_Straight_General, LineSource_Straight_CurrentsOnly, LineSource_Straight_ElectricOnly

###########################################################################
#                        PROPAGATION MEDIA
###########################################################################

abstract type PropagationMedia end

    struct PropagationMedia_Simple <: PropagationMedia
        ε::Quantity
        μ::Quantity
        c::Quantity
    end

    struct PropagationMedia_DiagonallyAnisotropic <: PropagationMedia
        ε::Diagonal{Quantity}
        μ::Diagonal{Quantity}
        c::Diagonal{Quantity}
    end

    export PropagationMedia_Simple, PropagationMedia_DiagonallyAnisotropic

###########################################################################
#                        JEFIMENKO MODELS
###########################################################################

struct JefimenkoModel{T}
    media::PropagationMedia
    sources::Vector{JefimenkoSource{T}}
    metadata::Dict{Symbol,Any}
end

export JefimenkoModel