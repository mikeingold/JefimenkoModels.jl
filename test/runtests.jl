using TestItemRunner
using TestItems

@run_package_tests verbose=true

@testsnippet Setup begin
    using Unitful
    using Unitful.DefaultSymbols: m, ns
    using UnitfulCoordinateSystems
    using JefimenkoModels
end

@testitem "LineSource_Straight" setup=[Setup] begin
    a = CoordinateCartesian(-0.5m, 0.0m, 0.0m)
    b = CoordinateCartesian( 0.5m, 0.0m, 0.0m)
    rho(r̄::AbstractCoordinate, t_s::Real) = cos(2π * 100e6 * t_s)^2
    J(r̄::AbstractCoordinate, t_s::Real) = x̂ .* cos(2π * 100e6 * t_s)
    source = LineSource_Straight{Float64}(a, b, rho, rho, J, J)

	@test source.a.x ≈ -0.5m
    @test source.b.x ≈  0.5m
end
