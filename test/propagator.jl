module TestPropagator

using Test
using PIMD, LinearAlgebra

ω, Δt, m = 2.3, 0.1, 1.4
A = [0 inv(m); - m * ω^2 0]

@testset "Propagator" begin
    @test exp(A * Δt) ≈ [cos(ω * Δt) sin(ω * Δt) * inv(ω * m); -sin(ω * Δt) * (ω * m) cos(ω * Δt)]
    @test sqrtcayley(A * Δt) ≈ [1 Δt / m; -ω^2 * Δt * m 1] / √(1 + ω^2 * Δt^2)
    @test cayley(A * Δt) ≈ [1 - ω^2 * Δt^2 / 4 Δt / m; -ω^2 * Δt * m 1 - ω^2 * Δt^2 / 4] / (1 + ω^2 * Δt^2 / 4)
end

end