using TestItemRunner
using TestItems

@run_package_tests verbose=true

@testsnippet Setup begin
    using JefimenkoModels
    using Meshes
    using Unitful
    using Unitful.DefaultSymbols: m, s
end

@testitem "Construct a Model" setup=[Setup] begin
    # Geometry
    a = Point(-0.5m, 0.0m, 0.0m)
    b = Point( 0.5m, 0.0m, 0.0m)
    segment = Segment(a, b)

    # Signals
    f = 100e6 / s
    ρe(r̄, t) = cos(2π * f * t)
    Je(r̄, t) = cos(2π * f * t) .* x̂ .* u"A"
    source = RadiationSource(segment, rho_e = ρe, J_e = Je)

    @test source.geometry === segment
    @test source.rho_e === ρe
    @test source.rho_h === NULL_CHARGE
    @test source.J_e === Je
    @test source.J_h === NULL_CURRENT
end
