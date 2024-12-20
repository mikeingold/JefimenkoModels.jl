################################################################################
#                             BUILD THE MODEL
################################################################################

using JefimenkoModels
using Meshes
using Unitful
using Unitful.DefaultSymbols: V, A, m, ns, s

model_disk = let
    # Surface disk with radius 0.5m
    ρ₀ = 0.5m
    origin = Point(0, 0, 0)
    ẑ = Vec(0, 0, 1)
    xy_plane = Plane(origin, ẑ)
    disk = Disk(origin, ρ₀)

    # Electric current only: spatially-uniform, x-directed, driven by a transient pulse
    (t₀, f₀, β₀) = (5.0ns, 500e6/s, 1.25)
    signal(t) = sin(2π * f₀ * t) * exp(-β₀ * (f₀ * t)^2)
    Je(r̄, t) = signal(t - t₀) .* x̂ .* A
    source = RadiationSource(disk, J_e = Je)

    metadata = Dict(
        :description => "Uniform current over a 0.5m disk, stimulated by transient pulse signal."
    )
    
    JefimenkoModel(CLASSICAL_VACUUM, source, metadata)
end

################################################################################
#                             RUN THE MODEL
################################################################################

# Observation location and time domain of interest
r = point(0.0m, 0.0m, 1.5m)
t = range(0.0ns, 20.0ns, length=800)

# Calculate the fields at r over the time domain
efield = map(t -> E(r, t, model_disk), t)
hfield = map(t -> H(r, t, model_disk), t)

################################################################################
#                             PLOT THE DATA
################################################################################

using Plots

# Accessor functions
e(i) = map(e -> getindex(e,i), efield)
h(i) = map(h -> getindex(h,i), hfield)

common_formatting = Dict(
    # Major grid
    :gridwidth => 1,
    :gridalpha => 0.2,
    :gridstyle => :dash,
    # Minor grid
    :minorgrid => true,
    :minorgridalpha => 0.15,
    :minorgridwidth => 1,
    :minorgridstyle => :dot,
)

jx(t) = model_disk.sources[1].Jₑ(r, t)[1]
jy(t) = model_disk.sources[1].Jₑ(r, t)[2]
jz(t) = model_disk.sources[1].Jₑ(r, t)[3]

# Plot the source current density
p1 = plot(t, [jx.(t), jy.(t), jz.(t)], label=["Jx" "Jy" "Jz"], linewidth=3,
        title="Source Current Density (Spatially-Uniform)",
        xlims=(0,20), xticks=0:4:20, xminorticks=4,
        ylims=(-1,1), yticks=-1:0.25:1, yminorticks=2,
        ; common_formatting...)
plot!(p1, xlabel="Time [ns]", ylabel="Magnitude [A/m]")
savefig(p1, joinpath(@__DIR__,"disk_fig_source.png"))

# Plot the electric field at the observation point
p2 = plot(t, [e(1), e(2), e(3)], label=["Ex" "Ey" "Ez"], linewidth=3,
        title="Electric Field (z = 1.5m)", framestyle=:zerolines,
        xlims=(0,20), xticks=0:4:20, xminorticks=4,
        ylims=(-200,150), yticks=-200:50:150, yminorticks=2,
        ; common_formatting...)
plot!(p2, xlabel="Time [ns]", ylabel="Magnitude [V/m]")
savefig(p2, joinpath(@__DIR__,"disk_fig_efield.png"))

# Plot the magnetic field at the observation point
p3 = plot(t, [h(1), h(2), h(3)], label=["Hx" "Hy" "Hz"], linewidth=3,
        title="Magnetic Field (z = 1.5m)", framestyle=:zerolines,
        xlims=(0,20), xticks=0:4:20, xminorticks=4,
        ylims=(-0.5,0.5), yticks=-0.5:0.1:0.5, yminorticks=2,
        ; common_formatting...)
plot!(p3, xlabel="Time [ns]", ylabel="Magnitude [A/m]")
savefig(p3, joinpath(@__DIR__,"disk_fig_hfield.png"))
