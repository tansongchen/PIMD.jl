module PIMD
"""
PIMD package implemented in Julia.
"""

using Distributions

export OneDimensionalQuantumSystem
export harmonicPotential, anharmonicPotential, quarticPotential
export ClassicalSystem, IsomorphicClassicalSystem
export Primitive, Staging, Normal
export HamiltonianPhase
export cayley, sqrtcayley
export md
export PotentialPropagator, ExternalPotentialPropagator, KineticPropagator, LangevinPropagator, GeneralizedCayleyPropagator
export Integrator
export StatisticalProperty
export kinetic, potential
export estimatorPotential, estimatorKinetic, estimatorKineticVirial, estimatorHessian
export MonteCarloSample
export statistics

include("propagator.jl")

MonteCarloSample(system::ClassicalSystem, nSample::Int) = begin
    space = 1000
    ninit = 10000
    momenta = MvNormal(system.M / system.β)
    uniform = Uniform(-1., 1.)
    phaseArray = []

    q = zeros(system.n)
    v = potential(system, q)
    for nStep = -ninit:(nSample * space)
        qTrial = copy(q)
        qTrial[rand(1:system.n)] += rand(uniform)
        vTrial = potential(system, qTrial)
        if rand() < min(1., exp(-system.β * (vTrial - v)))
            q .= qTrial
            v = vTrial
        end
        if nStep > 0 && nStep % space == 0
            push!(phaseArray, HamiltonianPhase(q, rand(momenta)))
        end
    end
    phaseArray
end

end
