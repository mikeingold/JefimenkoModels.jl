module JefimenkoModels
    using ForwardDiff
    using HCubature
    using LinearAlgebra
    using StaticArrays
    using Unitful
    using Unitful.DefaultSymbols: A, V, m, s, rad
    using UnitfulCoordinateSystems

    ###########################################################################
    #                        DATA STRUCTURES
    ###########################################################################

    abstract type JefimenkoSource end

        abstract type VolumeSource  <: JefimenkoSource end

            struct VolumeSource_Rectangular <: VolumeSource
                xlims::Tuple{Unitful.Length, Unitful.Length}
                ylims::Tuple{Unitful.Length, Unitful.Length}
                zlims::Tuple{Unitful.Length, Unitful.Length}
                ρₑ::Function
                ρₕ::Function
                Jₑ::Function
                Jₕ::Function
            end

            struct VolumeSource_Cylinder <: VolumeSource
                r::Tuple{Unitful.Length, Unitful.Length}
                ϕlims::Tuple{Unitful.Length, Unitful.Length}
                zlims::Tuple{Unitful.Length, Unitful.Length}
                ρₑ::Function
                ρₕ::Function
                Jₑ::Function
                Jₕ::Function
            end

            struct VolumeSource_Sphere <: VolumeSource
                r::Tuple{Unitful.Length, Unitful.Length}
                θlims::Tuple{Unitful.Length, Unitful.Length}
                ϕlims::Tuple{Unitful.Length, Unitful.Length}
                ρₑ::Function
                ρₕ::Function
                Jₑ::Function
                Jₕ::Function
            end

        abstract type SurfaceSource <: JefimenkoSource end

            struct SurfaceSource_Rectangle <: SurfaceSource
                xlims::Tuple{Unitful.Length, Unitful.Length}
                ylims::Tuple{Unitful.Length, Unitful.Length}
                ρₑ::Function
                ρₕ::Function
                Jₑ::Function
                Jₕ::Function
            end

            struct SurfaceSource_Disk <: SurfaceSource
                ρ₀::Unitful.Length
                ρₑ::Function
                ρₕ::Function
                Jₑ::Function
                Jₕ::Function
            end

        abstract type LineSource <: JefimenkoSource end

            struct LineSource_Straight <: LineSource
                ā::CoordinateCartesian
                b̄::CoordinateCartesian
                ρₑ::Function
                ρₕ::Function
                Jₑ::Function
                Jₕ::Function
            end

    export VolumeSource_Rectangular, VolumeSource_Cylinder, VolumeSource_Sphere
    export SurfaceSource_Rectangle, SurfaceSource_Disk
    export LineSource_Straight

    struct PropagationMedia
        ε::Quantity
        μ::Quantity
        c::Quantity
    end

    struct JefimenkoModel
        media::PropagationMedia
        sources::Vector{JefimenkoSource}
        metadata::Dict{Symbol,Any}
    end

    export PropagationMedia, JefimenkoModel

    ###########################################################################
    #                            FUNCTIONS
    ###########################################################################

    tᵣ(r̄::Coordinate, t::Unitful.Time, r̄′::Coordinate, c::Quantity)::Unitful.Time = t - (norm(r̄-r̄′)/c)

    ###########################################################################
    #                       𝐄-FIELD FUNCTIONS
    ###########################################################################

    function 𝐄(model::JefimenkoModel, r̄::Coordinate, t::Unitful.Time; rtol=sqrt(eps()))
        # Sum the contributions of the 𝐄(r̄,t) produced by each source in model
        E_contrib(source) = 𝐄(source, model.media, r̄, t; rtol=rtol)
        return mapreduce(E_contrib, +, model.sources) 
    end

    function 𝐄(source::SurfaceSource_Disk, media::PropagationMedia, r̄::Coordinate, t::Unitful.Time; rtol=sqrt(eps()))
        # Define the integrand function, returns Cartesian vector with implied units of A/(s*m)
        function integrand_u(r̄′::CoordinatePolar)::SVector{3,Float64}
            r̄′_cart = CoordinateCartesian(r̄′)
            tr_s::Float64 = ustrip(Float64, s, tᵣ(r̄,t,r̄′,media.c))   # in s
            Jₑ(t::Real) = source.Jₑ(r̄′_cart, t)                       # in A/m
            ∂Jₑ_∂t = ForwardDiff.derivative(Jₑ, tr_s)                # in A/(s*m)
            r_m = ustrip(m, norm(r̄ - r̄′_cart))                       # in m
            ρ_m = ustrip(m, r̄′.r)                                   # in m
            return ∂Jₑ_∂t .* (ρ_m / r_m)                       # in A/(s*m^2) times integration factor ρ in m
        end
        # Define an shim function since HCubature doesn't currently support Unitful
        integrand(coord) = integrand_u(CoordinatePolar(coord[1]*m, coord[2]*rad))

        coeff::Quantity = -1 / (4π * media.ε * media.c^2)

        # Integrate over circular aperture
        ρ₀_m = ustrip(m, source.ρ₀)
        iint = hcubature(integrand, [0.0, 0.0], [ρ₀_m, 2.0π], rtol=rtol)
        iint_u = iint[1] .* (A/s)  # dimensions after 2D integration are A/s
        return ( coeff .* iint_u ) .|> V/m  # convert to V/m
    end

    function 𝐄(model::JefimenkoModel, r̄::Coordinate, t::Unitful.Time; rtol=sqrt(eps()))
        mapreduce(source -> 𝐄(source, model.media, r̄, t), +, model.sources) 
    end

    export 𝐄

    ###########################################################################
    #                       𝐇-FIELD FUNCTIONS
    ###########################################################################

    function 𝐇(model::JefimenkoModel, r̄::Coordinate, t::Unitful.Time; rtol=sqrt(eps()))
        # Sum the contributions of the 𝐇(r̄,t) produced by each source in model
        return mapreduce(source -> 𝐇(source, model.media, r̄, t; rtol=rtol), +, model.sources) 
    end

    function 𝐇(source::SurfaceSource_Disk, media::PropagationMedia, r̄::Coordinate, t::Unitful.Time; rtol=sqrt(eps()))
        error("Not implemented yet.")
    end

    export 𝐇
end
