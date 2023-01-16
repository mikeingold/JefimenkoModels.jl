"""
    __integrand_H_R1(r̄′::Coordinate; source::LineSource{T}, media::PropagationMedia
                    r̄::Coordinate, t::Time)  where {T<:AbstractFloat}

Calculate the integrand function for the magnetic Jefimenko equation of a one-dimensional
`line source` at the `source point r̄′`. Parameterize the integrand function for a particular
`field source`, `propagation media`, and an observer positioned at space-time point (`r̄`,`t`).

# Arguments
- `r̄′::Coordinate`: spatial coordinate of the source point

# Parameters
- `r̄::UnitfulCoordinateSystems.Coordinate`: spatial coordinate at the observation point
- `t::Unitful.Time`: time at the observation point
- `source::JefimenkoSource`: the source model generating the magnetic field
- `media::PropagationMedia`: properties of the propagation media

# Returns
- `Vector{Quantity}`: the predicted vector-valued integrand value with appropriate units
"""
function __integrand_H_R1(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time,
                         source::AbstractLineSource{T}, media::PropagationMedia_Simple
                         ) where {T<:AbstractFloat}
    # Get spatial properties, in implicit units of meters
    r̄′ = CoordinateCartesian(r̄′)                 # ensure r̄′ is in Cartesian format
    r̄ = CoordinateCartesian(r̄)                   # ensure r̄′ is in Cartesian format
    Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′))        #  vector r̄-r̄′
    r_m = norm(Δr̄_m)                             #  magnitude |r̄-r̄′|

    # Get media properties, in implicit units as specified
    c = ustrip(T, m/s, media.c)               #  speed of light in [m/s]
    μ = ustrip(T, (V*s)/(A*m), media.μ)       #  permeability in [Vs/Am]

    # Calculate source-observer retarded time, in implicit units of seconds
    t′_s::T = ustrip(T, s, t′(r̄,t,r̄′,media))

    # Source functions
    ρₕ = source.ρₕ(r̄′, t′_s)                                           # in [Wb m^-1]
    ∂ρₕ_∂t = ForwardDiff.derivative(t_s -> source.ρₕ(r̄′,t_s), t′_s)    # in [Wb m^-1 s^-1]
    Jₑ = source.Jₑ(r̄′, t′_s)                                          # in [A]
    ∂Jₑ_∂t = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), t′_s)   # in [A s^-1]
    ∂Jₕ_∂t = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), t′_s)    # in [V s^-1]

    # Calculate first term
    term1a = ( (Δr̄_m ./ r_m^3) .* ρₕ )                   # [m/m³ * Wb/m] -> [Vs/m³]
    term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₕ_∂t )     # [m/m² * s/m * Wb/ms] -> [Vs/m³]
    term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₕ_∂t )           # [1/m * s²/m² * V/s]  -> [Vs/m³]
    term1  = ( (μ^-1) .* (term1a + term1b - term1c) )    # [Am/Vs * Vs/m³] -> [A/m²]
    
    # Calculate second term
    term2a = ( Jₑ ./ r_m^3 )                       # [A / m³] -> [A/m³]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₑ_∂t )   # [1/m² * s/m * A/s] -> [A/m³]
    term2  = cross((term2a + term2b), Δr̄_m)        # [A/m³ * m] -> [A/m²]

    # Combine terms, apply appropriate units
    return (term1 + term2) * (A/m^2)
end

"""
    __integrand_H_R2(r̄′::Coordinate; source::SurfaceSource{T}, media::PropagationMedia
                    r̄::Coordinate, t::Time) where {T<:AbstractFloat}

Calculate the integrand function for the magnetic Jefimenko equation of a two-dimensional
`surface source` at the `source point r̄′`. Parameterize the integrand function for a particular
`field source`, `propagation media`, and an observer positioned at space-time point (`r̄`,`t`).

# Arguments
- `r̄′::Coordinate`: spatial coordinate of the source point

# Parameters
- `r̄::UnitfulCoordinateSystems.Coordinate`: spatial coordinate at the observation point
- `t::Unitful.Time`: time at the observation point
- `source::JefimenkoSource`: the source model generating the magnetic field
- `media::PropagationMedia`: properties of the propagation media

# Returns
- `Vector{T}`: the predicted vector-valued integrand value, in implicit units [V/m³]. (The
`HCubature` solver does not yet support `Unitful` types.)
"""
function __integrand_H_R2(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time,
                         source::AbstractSurfaceSource{T}, media::PropagationMedia_Simple
                         ) where {T<:AbstractFloat}
    # Get spatial properties, in implicit units of meters
    r̄′ = CoordinateCartesian(r̄′)                 # ensure r̄′ is in Cartesian format
    r̄ = CoordinateCartesian(r̄)                   # ensure r̄′ is in Cartesian format
    Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′))        #  vector r̄-r̄′
    r_m = norm(Δr̄_m)                             #  magnitude |r̄-r̄′|

    # Get media properties, in implicit units as specified
    c = ustrip(T, m/s, media.c)               #  speed of light in [m/s]
    μ = ustrip(T, (V*s)/(A*m), media.μ)       #  permeability in [Vs/Am]

    # Calculate source-observer retarded time, in implicit units of seconds
    t′_s::T = ustrip(T, s, t′(r̄,t,r̄′,media))

    # Source functions
    ρₕ = source.ρₕ(r̄′, t′_s)                                           # in [Wb m^-2]
    ∂ρₕ_∂t = ForwardDiff.derivative(t_s -> source.ρₕ(r̄′,t_s), t′_s)    # in [Wb m^-2 s^-1]
    Jₑ = source.Jₑ(r̄′, t′_s)                                          # in [A m^-1]
    ∂Jₑ_∂t = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), t′_s)   # in [A m^-1 s^-1]
    ∂Jₕ_∂t = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), t′_s)    # in [V m^-1 s^-1]

    # Calculate first term
    term1a = ( (Δr̄_m ./ r_m^3) .* ρₕ )                   # [m/m³ * Wb/m²] -> [Vs/m⁴]
    term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₕ_∂t )     # [m/m² * s/m * Wb/m²s] -> [Vs/m⁴]
    term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₕ_∂t )           # [1/m * s²/m² * V/ms]  -> [Vs/m⁴]
    term1  = ( (μ^-1) .* (term1a + term1b - term1c) )    # [Am/Vs * Vs/m⁴] -> [A/m³]
    
    # Calculate second term
    term2a = ( Jₑ ./ r_m^3 )                       # [A/m / m³] -> [A/m⁴]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₑ_∂t )   # [1/m² * s/m * A/ms] -> [A/m⁴]
    term2  = cross((term2a + term2b), Δr̄_m)        # [A/m⁴ * m] -> [A/m³]

    # Combine terms and return, in implicit units of [A/m³]
    return (term1 + term2)
end

"""
    __integrand_H_R3(r̄′::Coordinate; source::VolumeSource{T}, media::PropagationMedia
                    r̄::Coordinate, t::Time) where {T<:AbstractFloat}

Calculate the integrand function for the magnetic Jefimenko equation of a two-dimensional
`surface source` at the `source point r̄′`. Parameterize the integrand function for a particular
`field source`, `propagation media`, and an observer positioned at space-time point (`r̄`,`t`).

# Arguments
- `r̄′::Coordinate`: spatial coordinate of the source point

# Parameters
- `r̄::UnitfulCoordinateSystems.Coordinate`: spatial coordinate at the observation point
- `t::Unitful.Time`: time at the observation point
- `source::JefimenkoSource`: the source model generating the magnetic field
- `media::PropagationMedia`: properties of the propagation media

# Returns
- `Vector{T}`: the predicted vector-valued integrand value, in implicit units [V/m³]. (The
`HCubature` solver does not yet support `Unitful` types.)
"""
function __integrand_H_R3(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time,
                         source::AbstractVolumeSource{T}, media::PropagationMedia_Simple
                         ) where {T<:AbstractFloat}
    # Get spatial properties, in implicit units of meters
    r̄′ = CoordinateCartesian(r̄′)                 # ensure r̄′ is in Cartesian format
    r̄ = CoordinateCartesian(r̄)                   # ensure r̄′ is in Cartesian format
    Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′))        #  vector r̄-r̄′
    r_m = norm(Δr̄_m)                             #  magnitude |r̄-r̄′|

    # Get media properties, in implicit units as specified
    c = ustrip(T, m/s, media.c)               #  speed of light in [m/s]
    μ = ustrip(T, (V*s)/(A*m), media.μ)       #  permeability in [Vs/Am]

    # Calculate source-observer retarded time, in implicit units of seconds
    t′_s::T = ustrip(T, s, t′(r̄,t,r̄′,media))

    # Source functions
    ρₕ = source.ρₕ(r̄′, t′_s)                                           # in [Wb m^-3]
    ∂ρₕ_∂t = ForwardDiff.derivative(t_s -> source.ρₕ(r̄′,t_s), t′_s)    # in [Wb m^-3 s^-1]
    Jₑ = source.Jₑ(r̄′, t′_s)                                          # in [A m^-2]
    ∂Jₑ_∂t = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), t′_s)   # in [A m^-2 s^-1]
    ∂Jₕ_∂t = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), t′_s)    # in [V m^-2 s^-1]

    # Calculate first term
    term1a = ( (Δr̄_m ./ r_m^3) .* ρₕ )                   # [m/m³ * Wb/m³] -> [Vs/m⁵]
    term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₕ_∂t )     # [m/m² * s/m * Wb/m³s] -> [Vs/m⁵]
    term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₕ_∂t )           # [1/m * s²/m² * V/m²s]  -> [Vs/m⁵]
    term1  = ( (μ^-1) .* (term1a + term1b - term1c) )    # [Am/Vs * Vs/m⁵] -> [A/m⁴]
    
    # Calculate second term
    term2a = ( Jₑ ./ r_m^3 )                       # [A/m² / m³] -> [A/m⁵]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₑ_∂t )   # [1/m² * s/m * A/m²s] -> [A/m⁵]
    term2  = cross((term2a + term2b), Δr̄_m)        # [A/m⁵ * m] -> [A/m⁴]

    # Combine terms and return, in implicit units of [A/m⁴]
    return (term1 + term2)
end
