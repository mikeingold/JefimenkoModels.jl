function _𝐇(r̄::Coordinate, t::Unitful.Time, source::SurfaceSource_Disk{T},
    media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
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