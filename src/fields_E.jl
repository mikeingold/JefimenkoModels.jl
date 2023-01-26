###########################################################################
#                        LINEAR SOURCES
###########################################################################

"""
    __E(r̄::AbstractCoordinate, t::Time, source::AbstractJefimenkoSource,
        media::PropagationMedia; rtol=sqrt(eps))

Calculate the electric field at (`r̄`,`t`) using the electric Jefimenko equation due to a
particular `source`, transmitted through a particular homogeneous `propagation media`.
Calculate the integral using a specified `relative tolerance`.

# Arguments
- `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the observation point
- `t::Unitful.Time`: time at which the electric field is observed
- `source::AbstractJefimenkoSource{T}`: source of the electric field
- `media::PropagationMedia`: properties of the propagation media

# Keywords
- `rtol::Real`: relative tolerance at which to solve the integral (optional)
"""
function __E(r̄::AbstractCoordinate, t::Unitful.Time, source::LineSource_Straight{T},
             media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    # Calculate the length of the line source from starting point ā to ending point b̄
    dmax::Unitful.Length = norm(source.b̄ - source.ā)

    # Calculate the integrand E-field vector in implied units [V/m²]
    function integrand_Vm2(u, p)::SVector{3,T}
        d::Unitful.Length = u * m
        
        # Parameterize a straight line from ā to b̄ according to the distance traveled
        # Start at ā, progress the specified distance in direction û
        û = (source.b̄ - source.ā) ./ dmax
        r̄′::CoordinateCartesian = source.ā + (d .* û)

        return __integrand_E(r̄′; source=source, media=media,
                             r̄=CoordinateCartesian(r̄), t=t)::SVector{3,T}
    end

    # Define the integrand as a f(d) traveled along line source, solve it
    prob = IntegralProblem(integrand_Vm2, zero(T), ustrip(T,m,dmax))
    sol = solve(prob, QuadGKJL(), reltol=rtol)
    return ( (1/4π) .* (sol.u) .* (V/m) )             # in [V/m² * m] -> [V/m]
end

###########################################################################
#                        SURFACE SOURCES
###########################################################################

function __E(r̄::AbstractCoordinate, t::Unitful.Time, source::SurfaceSource_Disk{T},
             media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    function disk_integrand_Vm2(u, p)::SVector{3,T}
        # Convert given (ρ[m],φ[rad]) to a Coordinate
        r̄′ = CoordinateCartesian(CoordinatePolar(u[1]*m, u[2]*rad))

        # Return integrand scaled by the radial integration factor
        return (__integrand_E(r̄′; source=source, media=media,
                             r̄=CoordinateCartesian(r̄), t=t) .* u[1])::SVector{3,T}
    end

    # Get integration limits: ρ ∈ [0,ρ₀], ϕ ∈ [0,2π]
    ρ₀_m = ustrip(T, m, source.ρ₀)
    lb = [zero(T), zero(T)]
    ub = [T(ρ₀_m), T(2π)]

    # Define and solve the integral problem over a circular aperture
    prob = IntegralProblem(disk_integrand_Vm2, lb, ub)
    sol = solve(prob, HCubatureJL(), reltol=rtol)      # implied units [V/m² * m] -> [V/m]
    return ( (1/4π) .* (sol.u) .* (V/m) )
end

function __E(r̄::AbstractCoordinate, t::Unitful.Time, source::SurfaceSource_Rectangle{T},
             media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    function integrand_Vm3(u, p)::SVector{3,T}
        r̄′ = CoordinateCartesian(u[1]*m, u[2]*m, 0.0m)

        return __integrand_E(r̄′; source=source, media=media,
                             r̄=CoordinateCartesian(r̄), t=t)::SVector{3,T}
    end

    # Get integration limits
    (lim_min_x_m, lim_max_x_m) = ustrip.(T, m, source.xlims)
    (lim_min_y_m, lim_max_y_m) = ustrip.(T, m, source.ylims)
    lb = [lim_min_x_m, lim_min_y_m]
    ub = [lim_max_x_m, lim_max_y_m]

    # Define and solve the integral problem over rectangular aperture
    prob = IntegralProblem(integrand_Vm3, lb, ub)
    sol = solve(prob, HCubatureJL(), reltol=rtol)     # implied units [V/m³ * m²] -> [V/m]
    return ( (1/4π) .* (sol.u) .* (V/m) )
end

###########################################################################
#                        VOLUME SOURCES
###########################################################################

function __E(r̄::AbstractCoordinate, t::Unitful.Time, source::VolumeSource_Cylinder{T},
             media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    error("Solver not yet implemented.")
end

function __E(r̄::AbstractCoordinate, t::Unitful.Time, source::VolumeSource_Rectangular{T},
             media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    function integrand_Vm4(u, p)::SVector{3,T}
        r̄′ = CoordinateCartesian(u[1]*m, u[2]*m, u[3]*m)

        return __integrand_E(r̄′; source=source, media=media,
                             r̄=CoordinateCartesian(r̄), t=t)::SVector{3,T}
    end

    # Get integration limits
    (lim_min_x_m, lim_max_x_m) = ustrip.(T, m, source.xlims)
    (lim_min_y_m, lim_max_y_m) = ustrip.(T, m, source.ylims)
    (lim_min_z_m, lim_max_z_m) = ustrip.(T, m, source.zlims)
    lb = [lim_min_x_m, lim_min_y_m, lim_min_z_m]
    ub = [lim_max_x_m, lim_max_y_m, lim_max_z_m]

    # Define and solve the integral problem over rectangular aperture
    prob = IntegralProblem(integrand_Vm4, lb, ub)
    sol = solve(prob, HCubatureJL(), reltol=rtol)     # implied units [V/m⁴ * m³] -> [V/m]
    return ( (1/4π) .* (sol.u) .* (V/m) )
end

function __E(r̄::AbstractCoordinate, t::Unitful.Time, source::VolumeSource_Sphere{T},
             media::PropagationMedia_Simple; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
    error("Solver not yet implemented.")
end
