###########################################################################
#                        LINEAR SOURCES
###########################################################################

"""
    _𝐄(r̄::Coordinate, t::Time, source::JefimenkoSource, media::PropagationMedia; rtol)

Calculate the electric field at (`r̄`,`t`) using the electric Jefimenko equation due to a
particular `source`, transmitted through a particular homogeneous `propagation media`.
Calculate the integral using a specified `relative tolerance`.

# Arguments
- `r̄::UnitfulCoordinateSystems.Coordinate`: spatial location of the observation point
- `t::Unitful.Time`: time at which the electric field is observed
- `source::JefimenkoSource`: source of the electric field
- `media::PropagationMedia`: properties of the propagation media

# Keywords
- `rtol::Real`: relative tolerance at which to solve the integral (optional)
"""
function _𝐄(r̄::Coordinate, t::Unitful.Time, source::LineSource_Straight{T},
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

        val = _integrand_E_R1(r̄′; source=source, media=media, r̄=r̄, t=t)
        return ustrip.(T, V/m^2, val)
    end

    # Define the integrand as a f(d) traveled along line source, solve it
    prob = IntegralProblem(integrand, zeros(T,1), [ustrip(T,m,dmax)])
    sol = solve(prob, HCubatureJL(), reltol=rtol)
    return ( (1/4π) .* (sol.u) .* (V/m) )             # in [V/m² * m] -> [V/m]
end

###########################################################################
#                        SURFACE SOURCES
###########################################################################

function _𝐄(r̄::Coordinate, t::Unitful.Time, source::SurfaceSource_Disk{T},
            media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    function disk_integrand(ū,p)
        # Assign aliases to ū values and convert to a Coordinate
        (ρ_m, ϕ_rad) = ū
        r̄′ = CoordinatePolar(ρ_m*m, ϕ_rad*rad)
        # Return integrand scaled by the radial integration factor,
        #   in implied units [V/m³ * m] -> [V/m²]
        return _integrand_E_R2(r̄′; r̄=r̄, t=t, source=source, media=media) * ρ_m
    end

    # Get integration limits: ρ ∈ [0,ρ₀], ϕ ∈ [0,2π]
    ρ₀_m = ustrip(T, m, source.ρ₀)
    lb = [zero(T), zero(T)]
    ub = [ρ₀_m, T(2π)]

    # Define and solve the integral problem over a circular aperture,
    #   in implied units [V/m² * m] -> [V/m]
    prob = IntegralProblem(disk_integrand, lb, ub)
    sol = solve(prob, HCubatureJL(), reltol=rtol)
    return ( (1/4π) .* (sol.u) .* (V/m) )
end

function _𝐄(r̄::Coordinate, t::Unitful.Time, source::SurfaceSource_Rectangle{T},
    media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
function integrand(u,p)
(x_m, y_m) = u
r̄′ = CoordinateCartesian(x_m*m, y_m*m, 0.0m)
return _integrand_E_R2(r̄′; r̄=r̄, t=t, source=source, media=media)  # implied [V/m³]
end

# Get integration limits
(lim_min_x_m, lim_max_x_m) = ustrip.(T, m, source.xlims)
(lim_min_y_m, lim_max_y_m) = ustrip.(T, m, source.ylims)
lb = [lim_min_x_m, lim_min_y_m]
ub = [lim_max_x_m, lim_max_y_m]

# Define and solve the integral problem over rectangular aperture
prob = IntegralProblem(integrand, lb, ub)
sol = solve(prob, HCubatureJL(), reltol=rtol)     # implied [V/m³ * m²] -> [V/m]
return ( (1/4π) .* (sol.u) .* (V/m) )
end

###########################################################################
#                        VOLUME SOURCES
###########################################################################

function _𝐄(r̄::Coordinate, t::Unitful.Time, source::VolumeSource_Cylinder{T},
    media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    error("Solver not yet implemented.")
end

function _𝐄(r̄::Coordinate, t::Unitful.Time, source::VolumeSource_Rectangular{T},
    media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    function integrand(u,p)
        (x_m, y_m, z_m) = u
        r̄′ = CoordinateCartesian(x_m*m, y_m*m, z_m)
        return _integrand_E_R3(r̄′; r̄=r̄, t=t, source=source, media=media)  # implied [V/m⁴]
    end

    # Get integration limits
    (lim_min_x_m, lim_max_x_m) = ustrip.(T, m, source.xlims)
    (lim_min_y_m, lim_max_y_m) = ustrip.(T, m, source.ylims)
    (lim_min_z_m, lim_max_z_m) = ustrip.(T, m, source.zlims)
    lb = [lim_min_x_m, lim_min_y_m, lim_min_z_m]
    ub = [lim_max_x_m, lim_max_y_m, lim_max_z_m]

    # Define and solve the integral problem over rectangular aperture
    prob = IntegralProblem(integrand, lb, ub)
    sol = solve(prob, HCubatureJL(), reltol=rtol)     # implied [V/m⁴ * m³] -> [V/m]
    return ( (1/4π) .* (sol.u) .* (V/m) )
end

function _𝐄(r̄::Coordinate, t::Unitful.Time, source::VolumeSource_Sphere{T},
    media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    error("Solver not yet implemented.")
end
