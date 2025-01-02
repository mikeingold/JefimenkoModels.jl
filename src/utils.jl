###########################################################################
#                     OFF-THE-SHELF PARTS
###########################################################################

CLASSICAL_VACUUM = let
    ε₀ = uconvert((A*s)/(V*m), float(ε_0))
    μ₀ = uconvert((V*s)/(A*m), float(μ_0))
    c₀ = uconvert(m/s, float(c_0))
    SimpleMedia(ε₀, μ₀, c₀)
end

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
function t′(r̄::Meshes.Point, t::Unitful.Time, r̄′::Meshes.Point, c::Quantity)
    return (t - (LinearAlgebra.norm(r̄ - r̄′) / c))
end

"""
    t′(r̄, t, r̄′, media)

Calculate the retarded-time at a source point `r̄′` for an observer at the space-time
point (`r̄`,`t`) through a `propagation medium`.

# Arguments
- `r̄::Meshes.Point`: spatial location of the observation point
- `t::Unitful.Time`: time at the observation point
- `r̄′::Meshes.Point`: spatial location of the source point
- `media::SimpleMedia`: properties of the medium between r̄′ and r̄
"""
function t′(r̄::Meshes.Point, t::Unitful.Time, r̄′::Meshes.Point, media::SimpleMedia)
    return t′(r̄, t, r̄′, media.c)
end

#=
function t′(r̄::Meshes.Point, t::Unitful.Time, r̄′::Meshes.Point, media::PropagationMedia_DiagonallyAnisotropic)
    Δr̄ = r̄ - r̄′
    Δt = norm(media.c^-1 * Δr̄) |> unit(t)
    return (t - Δt)
end
=#
