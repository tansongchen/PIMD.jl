module MDTest

using Test, PIMD

@testset "MD" begin
    phase = HamiltonianPhase([1., 2., 3., 4.])
    phase.z .+= 1
    @test phase.q ≈ [2., 3.]
    phase.p .= 0
    @test phase.z ≈ [2., 3., 0., 0.]
end

end