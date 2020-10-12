using PIMD
using Test
using LinearAlgebra

include("pi.jl")
include("md.jl")
include("system.jl")
include("propagator.jl")

β, m, V, n = 1.2, 3.4, quarticPotential, 16
quantum = OneDimensionalQuantumSystem(β, m, V)
primitive = Primitive(quantum, n)
energy(phase::HamiltonianPhase, system::ClassicalSystem) = kinetic(system, phase.p) + potential(system, phase.q)
Δt = .01
phase = HamiltonianPhase(.01randn(n), .01randn(n))
initialEnergy = energy(phase, primitive)
properties = [energy]
isEnergyConserving(expected, actual) = abs(expected - actual) / expected < .01

integrator1 = begin
    B = PotentialPropagator(Δt / 2)
    A = KineticPropagator(Δt)
    Integrator(B, A, B)
end

integrator2 = begin
    B = ExternalPotentialPropagator(Δt / 2)
    C = GeneralizedCayleyPropagator(primitive, exp, Δt)
    Integrator(B, C, B)
end

integrator3 = begin
    B = ExternalPotentialPropagator(Δt / 2)
    C = GeneralizedCayleyPropagator(primitive, sqrtcayley, Δt)
    Integrator(B, C, B)
end

@testset "PIMD Energy Conservation" begin
    for integrator in (integrator1, integrator2, integrator3)
        estimatedEnergy = md(phase, primitive, integrator, (10, 100), properties)[energy]
        @test isEnergyConserving(initialEnergy, estimatedEnergy)
    end
end
