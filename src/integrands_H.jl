"""
    _integrand_H_R1(r̄′::Coordinate; source::LineSource{T}, media::PropagationMedia
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
function _integrand_H_R1(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time,
                         source::LineSource{T}, media::PropagationMedia_Simple
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
    tr::Unitful.Time = tᵣ(r̄,t,r̄′,media)
    tr_s::T = ustrip(T, s, tr)

    # Source functions
    ρₕ = source.ρₕ(r̄′, tr_s)                                           # in [Wb m^-1]
    ∂ρₕ_∂t = ForwardDiff.derivative(t_s -> source.ρₕ(r̄′,t_s), tr_s)    # in [Wb m^-1 s^-1]
    Jₑ = source.Jₑ(r̄′, tr_s)                                          # in [A]
    ∂Jₑ_∂t = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), tr_s)   # in [A s^-1]
    #Jₕ = source.Jₕ(r̄′, tr_s)                                          # in [V]
    ∂Jₕ_∂t = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), tr_s)    # in [V s^-1]

    # Calculate first term
    term1a = ( (Δr̄_m ./ r_m^3) .* ρₕ )                   # [m/m³ * Wb/m] -> [Vs/m³]
    term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₕ_∂t )     # [m/m² * s/m * Wb/ms] -> [Vs/m³]
    term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₕ_∂t )           # [1/m * s²/m² * V/s]  -> [Vs/m³]
    term1  = ( (μ^-1) .* (term1a + term1b - term1c) )    # [Am/Vs * Vs/m³] -> [A/m²]
    
    # Calculate second term
    term2a = ( Jₑ ./ r_m^3 )                       # [A/m / m³] -> [A/m⁴]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₑ_∂t )   # [1/m² * s/m * A/ms] -> [A/m⁴]
    term2  = cross((term2a + term2b), Δr̄_m)        # [A/m⁴ * m] -> [A/m³]

    # Combine terms and return, in implicit units of [A/m²]
    return (term1 - term2)
end

"""
    _integrand_H_R2(r̄′::Coordinate; source::SurfaceSource{T}, media::PropagationMedia
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
function _integrand_H_R2(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time,
                         source::SurfaceSource{T}, media::PropagationMedia_Simple
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
    tr::Unitful.Time = tᵣ(r̄,t,r̄′,media)
    tr_s::T = ustrip(T, s, tr)

    # Source functions
    ρₕ = source.ρₕ(r̄′, tr_s)                                           # in [Wb m^-2]
    ∂ρₕ_∂t = ForwardDiff.derivative(t_s -> source.ρₕ(r̄′,t_s), tr_s)    # in [Wb m^-2 s^-1]
    Jₑ = source.Jₑ(r̄′, tr_s)                                          # in [A m^-1]
    ∂Jₑ_∂t = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), tr_s)   # in [A m^-1 s^-1]
    #Jₕ = source.Jₕ(r̄′, tr_s)                                          # in [V m^-1]
    ∂Jₕ_∂t = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), tr_s)    # in [V m^-1 s^-1]

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
    return (term1 - term2)
end

"""
    _integrand_H_R3(r̄′::Coordinate; source::VolumeSource{T}, media::PropagationMedia
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
function _integrand_H_R3(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time,
                         source::VolumeSource{T}, media::PropagationMedia_Simple
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
    tr::Unitful.Time = tᵣ(r̄,t,r̄′,media)
    tr_s::T = ustrip(T, s, tr)

    # Source functions
    ρₕ = source.ρₕ(r̄′, tr_s)                                           # in [Wb m^-3]
    ∂ρₕ_∂t = ForwardDiff.derivative(t_s -> source.ρₕ(r̄′,t_s), tr_s)    # in [Wb m^-3 s^-1]
    Jₑ = source.Jₑ(r̄′, tr_s)                                          # in [A m^-2]
    ∂Jₑ_∂t = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), tr_s)   # in [A m^-2 s^-1]
    #Jₕ = source.Jₕ(r̄′, tr_s)                                          # in [V m^-2]
    ∂Jₕ_∂t = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), tr_s)    # in [V m^-2 s^-1]

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
    return (term1 - term2)
end

#=

function 𝐈h(r̄′::Coordinate, source::SurfaceSource_Disk_General{T}; r̄::Coordinate, t::Unitful.Time, media::PropagationMedia_Simple) where {T<:AbstractFloat}
    r̄′_cart = CoordinateCartesian(r̄′)
    Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′_cart))
    r_m = norm(Δr̄_m)
    ρ_m = ustrip(T, m, UnitfulCoordinateSystems.ρ(r̄′))
    c = ustrip(T, m/s, media.c)
    μ = ustrip(T, V*s/(A*m), media.μ)

    # Calculate source-observer retarded time
    tr::Unitful.Time = tᵣ(r̄,t,r̄′,media)
    tr_s::T = ustrip(T, s, tr)        # retarded time in s

    # Source functions
    ρₕ(t::Real) = source.ρₕ(r̄′_cart, t)               # in Wb m^-2
    ∂ρₕ_∂t(t::Real) = ForwardDiff.derivative(ρₕ, t)   # in Wb m^-2 s^-1
    Jₕ(t::Real) = source.Jₕ(r̄′_cart, t)                # in V m^-1
    ∂Jₕ_∂t(t::Real) = ForwardDiff.derivative(Jₕ, t)    # in V m^-1 s^-1
    Jₑ(t::Real) = source.Jₑ(r̄′_cart, t)               # in A m^-1
    ∂Jₑ_∂t(t::Real) = ForwardDiff.derivative(Jₑ, t)   # in A m^-1 s^-1

    # Calculate first term
    term1a = ( (Δr̄_m ./ r_m^3) .* ρₕ(tr_s) )                # [m/m^3 * Wb/m^2]        -> [Vs/m^4]
    term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₕ_∂t(tr_s) )  # [m/m^2 * s/m * Wb/sm^2] -> [Vs/m^4]
    term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₕ_∂t(tr_s) )        # [1/m * s^2/m^2 * V/sm]  -> [Vs/m^4]
    term1  = ( (μ^-1) .* (term1a + term1b - term1c) )       # [Am/Vs * Vs/m^4] -> [A/m^3]
    
    # Calculate second term
    term2a = ( Jₑ(tr_s) ./ r_m^3 )                          # [A/m / m^3]          -> [A/m^4]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₑ_∂t(tr_s) )      # [1/m^2 * s/m * A/sm] -> [A/m^4]
    term2  = cross((term2a + term2b), Δr̄_m)                 # [A/m^4 * m]          -> [A/m^3]

    # Combine terms and apply integration factor
    return ( (term1 + term2) * ρ_m )  # [A/m^3 * m] -> [A/m^2]
end

function 𝐈h(r̄′::Coordinate, source::SurfaceSource_Disk_ElectricOnly{T}; r̄::Coordinate, t::Unitful.Time, media::PropagationMedia_Simple) where {T<:AbstractFloat}
    r̄′_cart = CoordinateCartesian(r̄′)
    Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′_cart))
    r_m = norm(Δr̄_m)
    ρ_m = ustrip(T, m, UnitfulCoordinateSystems.ρ(r̄′))
    c = ustrip(T, m/s, media.c)

    # Calculate source-observer retarded time
    tr::Unitful.Time = tᵣ(r̄,t,r̄′,media)
    tr_s::T = ustrip(T, s, tr)        # retarded time in s

    # Source functions
    Jₑ(t::Real) = source.Jₑ(r̄′_cart, t)               # in A m^-1
    ∂Jₑ_∂t(t::Real) = ForwardDiff.derivative(Jₑ, t)   # in A m^-1 s^-1
    
    # Calculate second term
    term2a = ( Jₑ(tr_s) ./ r_m^3 )                          # [A/m / m^3]          -> [A/m^4]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₑ_∂t(tr_s) )      # [1/m^2 * s/m * A/sm] -> [A/m^4]
    term2  = cross((term2a + term2b), Δr̄_m)                 # [A/m^4 * m]          -> [A/m^3]

    # Apply integration factor
    return ( term2 * ρ_m )  # [A/m^3 * m] -> [A/m^2]
end

function 𝐈h(r̄′::Coordinate, source::SurfaceSource_Disk_CurrentsOnly{T}; r̄::Coordinate, t::Unitful.Time, media::PropagationMedia_Simple) where {T<:AbstractFloat}
    r̄′_cart = CoordinateCartesian(r̄′)
    Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′_cart))
    r_m = norm(Δr̄_m)
    ρ_m = ustrip(T, m, UnitfulCoordinateSystems.ρ(r̄′))
    c = ustrip(T, m/s, media.c)
    μ = ustrip(T, V*s/(A*m), media.μ)

    # Calculate source-observer retarded time
    tr::Unitful.Time = tᵣ(r̄,t,r̄′,media)
    tr_s::T = ustrip(T, s, tr)        # retarded time in s

    # Source functions
    Jₕ(t::Real) = source.Jₕ(r̄′_cart, t)                # in V m^-1
    ∂Jₕ_∂t(t::Real) = ForwardDiff.derivative(Jₕ, t)    # in V m^-1 s^-1
    Jₑ(t::Real) = source.Jₑ(r̄′_cart, t)               # in A m^-1
    ∂Jₑ_∂t(t::Real) = ForwardDiff.derivative(Jₑ, t)   # in A m^-1 s^-1

    # Calculate first term
    term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₕ_∂t(tr_s) )        # [1/m * s^2/m^2 * V/sm]  -> [Vs/m^4]
    term1  = ( (μ^-1) .* (-term1c) )       # [Am/Vs * Vs/m^4] -> [A/m^3]
    
    # Calculate second term
    term2a = ( Jₑ(tr_s) ./ r_m^3 )                          # [A/m / m^3]          -> [A/m^4]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₑ_∂t(tr_s) )      # [1/m^2 * s/m * A/sm] -> [A/m^4]
    term2  = cross((term2a + term2b), Δr̄_m)                 # [A/m^4 * m]          -> [A/m^3]

    # Combine terms and apply integration factor
    return ( (term1 + term2) * ρ_m )  # [A/m^3 * m] -> [A/m^2]
end

=#
