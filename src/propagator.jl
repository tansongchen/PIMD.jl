include("system.jl")

@doc raw"""
    KineticPropagator(Δt)

Creates a propagator that updates positions according to the canonical equation

``\boldsymbol q\leftarrow\boldsymbol q+\frac{\partial H}{\partial\boldsymbol p}\Delta t``
"""
struct KineticPropagator <: Propagator
    Δt::Real
end
(propagator!::KineticPropagator)(phase::Phase, system::ClassicalSystem) = begin
    phase.q .+= gradient(p -> kinetic(system, p), phase.p) * propagator!.Δt
end

@doc raw"""
    KineticPropagator(Δt)

Creates a propagator that updates momenta according to the canonical equation

``\boldsymbol p\leftarrow\boldsymbol p-\frac{\partial H}{\partial\boldsymbol q}\Delta t``
"""
struct PotentialPropagator <: Propagator
    Δt::Real
end
(propagator!::PotentialPropagator)(phase::Phase, system::ClassicalSystem) = begin
    phase.p .-= gradient(q -> potential(system, q), phase.q) * propagator!.Δt
end

@doc raw"""
    ExternalKineticPropagator(Δt)

Creates a propagator that updates momenta according to the canonical equation, with PIMD external potential only

``\boldsymbol p\leftarrow\boldsymbol p-\frac{\partial V_{\rm ext}}{\partial\boldsymbol q}\Delta t``
"""
struct ExternalPotentialPropagator <: Propagator
    Δt::Real
end
(propagator!::ExternalPotentialPropagator)(phase::Phase, system::ClassicalSystem) = begin
    phase.p .-= gradient(q -> potentialExternal(system, q), phase.q) * propagator!.Δt
end

@doc raw"""
    LangevinPropagator(system, Δt)

Creates a propagator that thermostats momenta using Langevin equation

``\boldsymbol p\leftarrow e^{-\Gamma\Delta t}\boldsymbol p+\sqrt{1-e^{-2\Gamma\Delta t}}\sqrt{\frac{M}{\beta}}\boldsymbol{\xi}``

where ``\boldsymbol{\xi}`` is sampled from standard normal distribution and ``\Gamma`` is determined by the characteristic frequency.
"""
struct LangevinPropagator <: Propagator
    dissipation::Operator
    fluctuation::Operator
    LangevinPropagator(system::ClassicalSystem, Δt::Real) = begin
        Γ = friction(system)
        dissipation = exp(-Δt * Γ)
        fluctuation = sqrt(I - dissipation^2) * sqrt(system.M / system.β)
        new(dissipation, fluctuation)
    end
end
(propagator!::LangevinPropagator)(phase::Phase, system::ClassicalSystem) = begin
    phase.p .= propagator!.dissipation * phase.p + propagator!.fluctuation * randn(system.n)
end

friction(system::Primitive) = √(inv(system.M) * system.A)
friction(system::Union{Normal, Staging}) = Diagonal(√system.n / system.β * I, system.n)

@doc raw"""
    GeneralizedCayleyPropagator(system, generalizedCayleyFunction, Δt)

Creates a generalized Cayley propagator that  PIMD internal potential using a generalized cayley function ``\varphi``

``\boldsymbol z\leftarrow\varphi(\mathsf J\mathsf H)\boldsymbol z``

Examples of generalized Cayley functions are `exp()`, `cayley()` and `sqrtcayley()` (recommended).
"""
struct GeneralizedCayleyPropagator <: Propagator
    matrixRepresentation::Matrix
    GeneralizedCayleyPropagator(system::IsomorphicClassicalSystem, generalizedCayleyFunction::Function, Δt::Real) = begin
        n = system.n
        B = Matrix([0I inv(system.M); -system.A 0I]) * Δt
        new(generalizedCayleyFunction(B))
    end
end
(propagator!::GeneralizedCayleyPropagator)(phase::Phase, system::ClassicalSystem) = begin
    z = vcat(phase.q, phase.p)
    z .= propagator!.matrixRepresentation * z
    phase.q .= z[1:system.n]
    phase.p .= z[system.n+1:2(system.n)]
end

const cayley = A::Matrix -> (I + A / 2) / (I - A / 2)
const sqrtcayley = A::Matrix -> (real ∘ √)((I + A)/(I - A))
