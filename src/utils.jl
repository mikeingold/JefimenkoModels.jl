###########################################################################
#                     OFF-THE-SHELF PARTS
###########################################################################

CLASSICAL_VACUUM = let
    ε₀ = uconvert((A*s)/(V*m), float(ε_0))
    μ₀ = uconvert((V*s)/(A*m), float(μ_0))
    c₀ = uconvert(m/s, float(c_0))
    PropagationMedia_Simple(ε₀, μ₀, c₀)
end

NULL_CHARGE(r̄::Meshes.Point, t_s::Real) = 0
NULL_CURRENT(r̄::Meshes.Point, t_s::Real) = StaticArrays.SVector(0, 0, 0)

###########################################################################
#                     RETARDED-TIME CALCULATIONS
###########################################################################

"""
    t′(r̄, t, r̄′, c)

Calculate the retarded-time at a source point `r̄′` for an observer at the space-time
point (`r̄`,`t`) through a medium with speed of light `c`.

# Arguments
- `r̄::Meshes.Point`: spatial location of the observation point
- `t::Unitful.Time`: time at the observation point
- `r̄′::Meshes.Point`: spatial location of the source point
- `c::Quantity`: Unitful speed of light in the medium between r̄′ and r̄
"""
function t′(r̄::AbstractCoordinate, t::Unitful.Time, r̄′::AbstractCoordinate, c::Quantity)
    return (t - (norm(r̄-r̄′)/c))
end

"""
    t′(r̄, t, r̄′, media)

Calculate the retarded-time at a source point `r̄′` for an observer at the space-time
point (`r̄`,`t`) through a `propagation medium`.

# Arguments
- `r̄::Meshes.Point`: spatial location of the observation point
- `t::Unitful.Time`: time at the observation point
- `r̄′::Meshes.Point`: spatial location of the source point
- `media::AbstractPropagationMedia`: properties of the medium between r̄′ and r̄
"""
function t′(r̄::Meshes.Point, t::Unitful.Time, r̄′::Meshes.Point, media::PropagationMedia_Simple)
    return t′(r̄, t, r̄′, media.c)
end

function t′(r̄::Meshes.Point, t::Unitful.Time, r̄′::Meshes.Point, media::PropagationMedia_DiagonallyAnisotropic)
    Δr̄ = SVector(r̄ - r̄′)
    Δt = norm(media.c^-1 * Δr̄) |> unit(t)
    return (t - Δt)
end
