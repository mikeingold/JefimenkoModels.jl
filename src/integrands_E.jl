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