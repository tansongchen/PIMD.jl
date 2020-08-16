include("system.jl")

struct KineticPropagator <: Propagator
    dt::Real
end
(propagator!::KineticPropagator)(phase::Phase, system::ClassicalSystem) = begin
    phase.q .+= gradient(p -> kinetic(system, p), phase.p) * propagator!.dt
end

struct PotentialPropagator <: Propagator
    dt::Real
end
(propagator!::PotentialPropagator)(phase::Phase, system::ClassicalSystem) = begin
    phase.p .-= gradient(q -> potential(system, q), phase.q) * propagator!.dt
end

struct ExternalPotentialPropagator <: Propagator
    dt::Real
end
(propagator!::ExternalPotentialPropagator)(phase::Phase, system::ClassicalSystem) = begin
    phase.p .-= gradient(q -> potentialExternal(system, q), phase.q) * propagator!.dt
end

struct LangevinPropagator <: Propagator
    dissipation::Operator
    fluctuation::Operator
    LangevinPropagator(system::ClassicalSystem, dt::Real) = begin
        Γ = friction(system)
        dissipation = exp(-dt * Γ)
        fluctuation = sqrt(I - dissipation^2) * sqrt(system.M / system.β)
        new(dissipation, fluctuation)
    end
end
(propagator!::LangevinPropagator)(phase::Phase, system::ClassicalSystem) = begin
    phase.p .= propagator!.dissipation * phase.p + propagator!.fluctuation * randn(system.n)
end

friction(system::Primitive) = √(inv(system.M) * system.A)
friction(system::Union{Normal, Staging}) = √system.n / system.β

struct GeneralizedCayleyPropagator <: Propagator
    matrixRepresentation::Matrix
    GeneralizedCayleyPropagator(system::IsomorphicClassicalSystem, generalizedCayleyFunction::Function, dt::Real) = begin
        n = system.n
        B = [0I inv(system.M); -system.A 0I] * dt
        new(generalizedCayleyFunction(B))
    end
end
(propagator!::GeneralizedCayleyPropagator)(phase::Phase, ::ClassicalSystem) = begin
    phase.z .= propagator!.matrixRepresentation * phase.z
end

const cayley = A::Matrix -> (I + A / 2) / (I - A / 2)
const sqrtcayley = A::Matrix -> (real ∘ √)((I + A)/(I - A))
