module MDTest

using Test, PIMD

@testset "MD" begin
    phase = HamiltonianPhase([1., 2.], [3., 4.])
    @test phase.q ≈ [1., 2.]
    phase.p .= 0
    @test phase.p ≈ [0., 0.]
end

end