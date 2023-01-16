"""
    _𝐇(r̄::Coordinate, t::Time, source::JefimenkoSource, media::PropagationMedia; rtol)

Calculate the magnetic field at (`r̄`,`t`) using the electric Jefimenko equation due to a
particular `source`, transmitted through a particular homogeneous `propagation media`.
Calculate the integral using a specified `relative tolerance`.

# Arguments
- `r̄::UnitfulCoordinateSystems.Coordinate`: spatial location of the observation point
- `t::Unitful.Time`: time at which the magnetic field is observed
- `source::JefimenkoSource`: source of the magnetic field
- `media::PropagationMedia`: properties of the propagation media

# Keywords
- `rtol::Real`: relative tolerance at which to solve the integral (optional)
"""
function _𝐇(r̄::Coordinate, t::Unitful.Time, source::LineSource_Straight_General{T},
            media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    # Calculate the length of the line source from starting point ā to ending point b̄
    dmax::Unitful.Length = norm(source.b̄ - source.ā)

    function integrand(u,p)
        d::Unitful.Length = u[1] * m
        # Parameterize a straight line from ā to b̄ according to the distance traveled
        # Get a unit vector pointing from ā -> b̄
        û = (source.b̄ - source.ā) ./ dmax

        #= TODO Replace the following hack
        #       Can't add a Coordinate to an SVector
        #       Implement SVector(Coordinate) and Coordinate(SVector)
        #       Implement +(::CoordinateCartesian, ::SVector{T,3}) where T<:Real
        # Start at ā, progress the specified distance in direction û
        source.ā + ( d .* û )
        =#
        trek = d .* û
        r̄′ = CoordinateCartesian(trek[1], trek[2], trek[3])

        val = _integrand_H_R1(r̄′; source=source, media=media, r̄=r̄, t=t)
        ustrip.(T, A/m^2, val)
    end

    # Define the integrand as a f(d) traveled along line source, solve it
    prob = IntegralProblem(integrand, zeros(T,1), [ustrip(T,m,dmax)])
    sol = solve(prob, HCubatureJL(), reltol=rtol)
    return ( (1/4π) .* (sol.u) .* (A/m) )             # in [A/m² * m] -> [A/m]
end

function _𝐇(r̄::Coordinate, t::Unitful.Time, source::SurfaceSource_Disk_General{T},
    media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    function disk_integrand(ū,p)
        # Assign aliases to ū values and convert to a Coordinate
        (ρ_m, ϕ_rad) = ū
        r̄′ = CoordinatePolar(ρ_m*m, ϕ_rad*rad)
        # Return integrand scaled by the radial integration factor,
        #   in implied units [V/m³ * m] -> [V/m²]
        return _integrand_E_R2(r̄′; r̄=r̄, t=t, source=source, media=media) * ρ_m
    end

    # Define and solve the integral problem over a circular aperture,
    #   in implied units [V/m² * m] -> [V/m]
    ρ₀_m = ustrip(T, m, source.ρ₀)
    lb = [zero(T), zero(T)]
    ub = [ρ₀_m, T(2π)]
    prob = IntegralProblem(disk_integrand, lb, ub)
    sol = solve(prob, HCubatureJL(), reltol=rtol)
    return ( (1/4π) .* (sol.u) .* (A/m) )
end
