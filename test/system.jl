module TestSystem

using Test, PIMD

β, m, V, n = 1.4, 2.3, harmonicPotential, 3
ω2 = n / β^2
quantum = OneDimensionalQuantumSystem(β, m, V)

@testset "Primitive" begin
    primitive = Primitive(quantum, n)
    @test primitive.A ≈ [2 -1 -1; -1 2 -1; -1 -1 2] * m * n / β^2
end

@testset "Staging" begin
    staging = Staging(quantum, n)
    @test staging.A ≈ [0 0 0; 0 ω2 0; 0 0 ω2] * staging.M
end

@testset "Normal" begin
    normal = Normal(quantum, n)
    @test normal.A ≈ [0 0 0; 0 ω2 0; 0 0 ω2] * normal.M
end

end