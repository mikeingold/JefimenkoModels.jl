module JefimenkoModels
    using ForwardDiff
    using HCubature
    using LinearAlgebra
    using StaticArrays
    using Unitful
    using Unitful.DefaultSymbols: A, V, C, m, s, rad
    using UnitfulCoordinateSystems

    ###########################################################################
    #                        DATA STRUCTURES
    ###########################################################################

    abstract type JefimenkoSource{T} end

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

    export VolumeSource_Rectangular, VolumeSource_Cylinder, VolumeSource_Sphere
    export SurfaceSource_Rectangle
    export SurfaceSource_Disk_General,  SurfaceSource_Disk_CurrentsOnly,  SurfaceSource_Disk_ElectricOnly
    export LineSource_Straight_General, LineSource_Straight_CurrentsOnly, LineSource_Straight_ElectricOnly

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

    struct JefimenkoModel{T}
        media::PropagationMedia
        sources::Vector{JefimenkoSource{T}}
        metadata::Dict{Symbol,Any}
    end

    export PropagationMedia, JefimenkoModel

    ###########################################################################
    #                            INTERNAL FUNCTIONS
    ###########################################################################

    tᵣ(r̄::Coordinate, t::Unitful.Time, r̄′::Coordinate, c::Quantity)::Unitful.Time = t - (norm(r̄-r̄′)/c)

    function 𝐈e(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time, media::PropagationMedia_Simple, source::SurfaceSource_Disk_General{T}) where {T<:AbstractFloat}
        r̄′_cart = CoordinateCartesian(r̄′)
        Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′_cart))
        r_m = norm(Δr̄_m)
        ρ_m = ustrip(T, m, UnitfulCoordinateSystems.ρ(r̄′))
        c = ustrip(T, m/s, media.c)
        ε = ustrip(T, A*s/(V*m), media.ε)

        tr_s::T = ustrip(T, s, tᵣ(r̄,t,r̄′,media.c))        # retarded time in s

        # Source functions
        ρₑ(t::Real) = source.ρₑ(r̄′_cart, t)               # in C m^-2
        ∂ρₑ_∂t(t::Real) = ForwardDiff.derivative(ρₑ, t)   # in C m^-2 s^-1
        Jₑ(t::Real) = source.Jₑ(r̄′_cart, t)               # in A m^-1
        ∂Jₑ_∂t(t::Real) = ForwardDiff.derivative(Jₑ, t)   # in A m^-1 s^-1
        Jₕ(t::Real) = source.Jₕ(r̄′_cart, t)                # in V m^-1
        ∂Jₕ_∂t(t::Real) = ForwardDiff.derivative(Jₕ, t)    # in V m^-1 s^-1

        # Calculate first term
        term1a = ( (Δr̄_m ./ r_m^3) .* ρₑ(tr_s) )                # [m/m^3 * C/m^2]         -> [A*s/m^4]
        term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₑ_∂t(tr_s) )  # [m/m^2 * s/m * C/sm^-2] -> [A*s/m^4]
        term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₑ_∂t(tr_s) )        # [1/m * s^2/m^2 * A/sm]  -> [A*s/m^4]
        term1  = ( (ε^-1) .* (term1a + term1b - term1c) )       # [Vm/As * As/m^4] -> [V/m^3]
        
        # Calculate second term
        term2a = ( Jₕ(tr_s) ./ r_m^3 )                          # [V/m / m^3] -> [V/m^4]
        term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₕ_∂t(tr_s) )      # [1/m^2 * s/m * V/sm] -> [V/m^4]
        term2  = cross((term2a + term2b), Δr̄_m)                 # [V/m^4 * m] -> [V/m^3]

        # Combine terms and apply integration factor
        return ( (term1 - term2) * ρ_m )  # [V/m^3 * m] -> [V/m^2]
    end

    function 𝐈e(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time, media::PropagationMedia_Simple, source::SurfaceSource_Disk_ElectricOnly{T}) where {T<:AbstractFloat}
        r̄′_cart = CoordinateCartesian(r̄′)
        Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′_cart))
        r_m = norm(Δr̄_m)
        ρ_m = ustrip(T, m, UnitfulCoordinateSystems.ρ(r̄′))
        c = ustrip(T, m/s, media.c)
        ε = ustrip(T, A*s/(V*m), media.ε)

        tr_s::T = ustrip(T, s, tᵣ(r̄,t,r̄′,media.c))        # retarded time in s

        # Source functions
        ρₑ(t::Real) = source.ρₑ(r̄′_cart, t)               # in C m^-2
        ∂ρₑ_∂t(t::Real) = ForwardDiff.derivative(ρₑ, t)   # in C m^-2 s^-1
        Jₑ(t::Real) = source.Jₑ(r̄′_cart, t)               # in A m^-1
        ∂Jₑ_∂t(t::Real) = ForwardDiff.derivative(Jₑ, t)   # in A m^-1 s^-1

        # Calculate first term
        term1a = ( (Δr̄_m ./ r_m^3) .* ρₑ(tr_s) )                # [m/m^3 * C/m^2]         -> [A*s/m^4]
        term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₑ_∂t(tr_s) )  # [m/m^2 * s/m * C/sm^-2] -> [A*s/m^4]
        term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₑ_∂t(tr_s) )        # [1/m * s^2/m^2 * A/sm]  -> [A*s/m^4]
        term1  = ( (ε^-1) .* (term1a + term1b - term1c) )       # [Vm/As * As/m^4] -> [V/m^3]

        # Apply integration factor
        return ( term1 * ρ_m )  # [V/m^3 * m] -> [V/m^2]
    end

    function 𝐈e(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time, media::PropagationMedia_Simple, source::SurfaceSource_Disk_CurrentsOnly{T}) where {T<:AbstractFloat}
        r̄′_cart = CoordinateCartesian(r̄′)
        Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′_cart))
        r_m = norm(Δr̄_m)
        ρ_m = ustrip(T, m, UnitfulCoordinateSystems.ρ(r̄′))
        c = ustrip(T, m/s, media.c)
        ε = ustrip(T, A*s/(V*m), media.ε)

        tr_s::T = ustrip(T, s, tᵣ(r̄,t,r̄′,media.c))        # retarded time in s

        # Source functions
        Jₑ(t::Real) = source.Jₑ(r̄′_cart, t)               # in A m^-1
        ∂Jₑ_∂t(t::Real) = ForwardDiff.derivative(Jₑ, t)   # in A m^-1 s^-1
        Jₕ(t::Real) = source.Jₕ(r̄′_cart, t)                # in V m^-1
        ∂Jₕ_∂t(t::Real) = ForwardDiff.derivative(Jₕ, t)    # in V m^-1 s^-1

        # Calculate first term
        term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₑ_∂t(tr_s) )    # [1/m * s^2/m^2 * A/sm]  -> [A*s/m^4]
        term1  = ( (ε^-1) .* (-term1c) )                    # [Vm/As * As/m^4] -> [V/m^3]
        
        # Calculate second term
        term2a = ( Jₕ(tr_s) ./ r_m^3 )                          # [V/m / m^3] -> [V/m^4]
        term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₕ_∂t(tr_s) )      # [1/m^2 * s/m * V/sm] -> [V/m^4]
        term2  = cross((term2a + term2b), Δr̄_m)                 # [V/m^4 * m] -> [V/m^3]

        # Combine terms and apply integration factor
        return ( (term1 - term2) * ρ_m )  # [V/m^3 * m] -> [V/m^2]
    end

    ###########################################################################
    #                       𝐄-FIELD FUNCTIONS
    ###########################################################################

    function 𝐄(r̄::Coordinate, t::Unitful.Time, model::JefimenkoModel; rtol=sqrt(eps()))
        # Sum the contributions of the 𝐄(r̄,t) produced by each source in model
        E_contrib(source) = 𝐄(r̄, t, source; media=model.media, rtol=rtol)
        return mapreduce(E_contrib, +, model.sources) 
    end

    function 𝐄(r̄::Coordinate, t::Unitful.Time, source::SurfaceSource_Disk{T}; media::PropagationMedia, rtol=sqrt(eps())) where {T<:AbstractFloat}
        # Define an shim function since HCubature doesn't currently support Unitful integration
        function integrand(coord)
            # coord -> [ρ in m, ϕ in rad]
            r̄′ = CoordinatePolar(coord[1]*m, coord[2]*rad)
            return 𝐈e(r̄′; r̄=r̄, t=t, media=media, source=source)
        end

        # Integrate over circular aperture.   [V/m^2 * m * []] -> [V/m]
        ρ₀_m = ustrip(T, m, source.ρ₀)
        iint = hcubature(integrand, [zero(T), zero(T)], [ρ₀_m, T(2π)], rtol=rtol)
        return ( (1/4π) .* iint[1] .* (V/m) )
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
