module TestPI

using Test, PIMD

@testset "PI" begin
    β, m, V = 1., 1., harmonicPotential
    system = OneDimensionalQuantumSystem(β, m, V)
    @test system.potential(10.) ≈ 50.
    @test system.force(10.) ≈ -10.
    @test system.hessian(10.) ≈ 1.
end

end