module JefimenkoModels
    using ForwardDiff
    using HCubature
    using LinearAlgebra
    using StaticArrays
    using Unitful
    using Unitful.DefaultSymbols: W, A, V, C, m, s, rad
    using UnitfulCoordinateSystems

    __DEFAULT_RTOL = sqrt(eps())

    include("structs.jl")

    ###########################################################################
    #                      RETARDED TIME CALCULATIONS
    ###########################################################################

    tᵣ(r̄::Coordinate, t::Unitful.Time, r̄′::Coordinate, c::Quantity)::Unitful.Time = t - (norm(r̄-r̄′)/c)

    tᵣ(r̄::Coordinate, t::Unitful.Time, r̄′::Coordinate, media::PropagationMedia_Simple)::Unitful.Time = t - (norm(r̄-r̄′)/media.c)

    function tᵣ(r̄::Coordinate, t::Unitful.Time, r̄′::Coordinate, media::PropagationMedia_DiagonallyAnisotropic)::Unitful.Time
        Δr̄ = SVector(r̄ - r̄′)
        Δt = norm(media.c^-1 * Δr̄) |> unit(t)
        return (t - Δt)
    end

    ###########################################################################
    #                       𝐄-FIELD FUNCTIONS
    ###########################################################################

    include("integrands_E.jl")

    function 𝐄(r̄::Coordinate, t::Unitful.Time, model::JefimenkoModel; rtol=__DEFAULT_RTOL)
        # Sum the contributions of the 𝐄(r̄,t) produced by each source in model
        E_contrib(source) = 𝐄(r̄, t, source; media=model.media, rtol=rtol)
        return mapreduce(E_contrib, +, model.sources) 
    end

    function 𝐄(r̄::Coordinate, t::Unitful.Time, source::SurfaceSource_Disk{T}; media::PropagationMedia, rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
        # Define an shim function since HCubature doesn't currently support Unitful integration
        function integrand(coord)
            # coord -> [ρ in m, ϕ in rad]
            r̄′ = CoordinatePolar(coord[1]*m, coord[2]*rad)
            return 𝐈e(r̄′, source; r̄=r̄, t=t, media=media)
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

    include("integrands_H.jl")

    function 𝐇(r̄::Coordinate, t::Unitful.Time, model::JefimenkoModel; rtol=__DEFAULT_RTOL)
        # Sum the contributions of the 𝐇(r̄,t) produced by each source in model
        H_contrib(source) = 𝐇(r̄, t, source; media=model.media, rtol=rtol)
        return mapreduce(H_contrib, +, model.sources) 
    end

    function 𝐇(r̄::Coordinate, t::Unitful.Time, source::SurfaceSource_Disk{T}; media::PropagationMedia, rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
        # Define an shim function since HCubature doesn't currently support Unitful integration
        function integrand(coord)
            # coord -> [ρ in m, ϕ in rad]
            r̄′ = CoordinatePolar(coord[1]*m, coord[2]*rad)
            return 𝐈h(r̄′, source; r̄=r̄, t=t, media=media)
        end

        # Integrate over circular aperture.   [A/m^2 * m * []] -> [A/m]
        ρ₀_m = ustrip(T, m, source.ρ₀)
        iint = hcubature(integrand, [zero(T), zero(T)], [ρ₀_m, T(2π)], rtol=rtol)
        return ( (1/4π) .* iint[1] .* (A/m) )
    end

    export 𝐇

    ###########################################################################
    #                       POYNTING VECTOR FUNCTIONS
    ###########################################################################

    function 𝐏(r̄::Coordinate, t::Unitful.Time, model::JefimenkoModel; rtol=__DEFAULT_RTOL)
        E = 𝐄(r̄,t,model; rtol=rtol)
        H = 𝐇(r̄,t,model; rtol=rtol)
        return cross(E,H) .|> W/m^2
    end

    function 𝐏(r̄::Coordinate, t::Unitful.Time, source::SurfaceSource_Disk{T}; media::PropagationMedia, rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
        E = 𝐄(r̄,t,source; media=media, rtol=rtol)
        H = 𝐇(r̄,t,source; media=media, rtol=rtol)
        return cross(E,H) .|> W/m^2
    end

    export 𝐏
end
