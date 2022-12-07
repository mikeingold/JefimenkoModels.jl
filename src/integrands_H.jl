###########################################################################
#              SurfaceSource_Disk in PropagationMedia_Simple
###########################################################################

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