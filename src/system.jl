include("pi.jl")
include("md.jl")

Operator = Union{Matrix, Diagonal, Symmetric, UniformScaling}

abstract type IsomorphicClassicalSystem <: ClassicalSystem end

kinetic(system::IsomorphicClassicalSystem, p::AbstractArray) = .5p' * inv(system.M) * p
potentialInternal(system::IsomorphicClassicalSystem, q::AbstractArray) = .5q' * system.A * q
potentialExternal(system::IsomorphicClassicalSystem, q::AbstractArray) = sum(system.quantum.potential, system.S * q) / system.n
potential(system::IsomorphicClassicalSystem, q::AbstractArray) = potentialInternal(system, q) + potentialExternal(system, q)

"""
    Primitive(quantum, n)

Creates an primitive isomorphic PIMD system from a quantum system, with beads number `n`.
"""
struct Primitive <: IsomorphicClassicalSystem
    n::Integer
    β::Real
    M::Operator
    A::Operator
    S::Operator
    quantum::OneDimensionalQuantumSystem
    Primitive(quantum::OneDimensionalQuantumSystem, n::Int) = begin
        M = Diagonal(quantum.m / n * I, n)
        A = stiffnessMatrix(quantum, n)
        S = I
        new(n, quantum.β, M, A, S, quantum)
    end
end

"""
    Staging(quantum, n)

Creates an staging isomorphic PIMD system from a quantum system, with beads number `n`.
"""
struct Staging <: IsomorphicClassicalSystem
    n::Integer
    β::Real
    M::Operator
    A::Operator
    S::Operator
    quantum::OneDimensionalQuantumSystem
    Staging(quantum::OneDimensionalQuantumSystem, n::Int) = begin
        M = Diagonal([quantum.m * i / (i - 1) for i in 1:n])
        M[1,1] = quantum.m
        K = stiffnessMatrix(quantum, n)
        S = zeros(n,n)
        S[:,1] .= 1.
        for i in 2:n
            S[i,i:n] = (i-1) ./ ((i-1):(n-1))
        end
        A = Diagonal(S' * K * S)
        new(n, quantum.β, M, A, S, quantum)
    end
end

"""
    Normal(quantum, n)

Creates an normal isomorphic PIMD system from a quantum system, with beads number `n`.
"""
struct Normal <: IsomorphicClassicalSystem
    n::Integer
    β::Real
    M::Diagonal
    A::Operator
    S::Operator
    quantum::OneDimensionalQuantumSystem
    Normal(quantum::OneDimensionalQuantumSystem, n::Int) = begin
        a = 0:n-1
        M = Diagonal([quantum.m * 4 * sin(i * π / n)^2 for i in a])
        M[1,1] = quantum.m
        K = stiffnessMatrix(quantum, n)
        S = zeros(n,n)
        S[:,1] .= sqrt(1 / n)
        for i in 1:(n-1)÷2
            c1, c2 = i+1, n+1-i
            @. S[:,c1] = sqrt(2 / n) * cos(2π * i * a / n)
            @. S[:,c2] = sqrt(2 / n) * sin(2π * i * a / n)
        end
        if n % 2 == 0
            @. S[:,n÷2+1] = sqrt(1 / n) * (-1) ^ a
        end
        A = Diagonal(S' * K * S)
        new(n, quantum.β, M, A, S, quantum)
    end    
end

@doc raw"""
    estimatorPotential(phase, system)

Estimates the potential energy for an isomorphic classical system under a phase.

``\hat V(\boldsymbol x)=\varphi(\boldsymbol x)``
"""
const estimatorPotential = (phase::HamiltonianPhase, system::IsomorphicClassicalSystem) -> sum(system.quantum.potential, system.S * phase.q) / system.n
@doc raw"""
    estimatorKinetic(phase, system)

Estimates the kinetic energy (using the primitive estimator) for an isomorphic classical system under a phase.

``\hat K_{\rm p}(\boldsymbol x)=\frac n{2\beta}-\frac12\boldsymbol x^T\mathsf K\boldsymbol x``
"""
const estimatorKinetic = (phase::HamiltonianPhase, system::IsomorphicClassicalSystem) -> system.n / 2system.β - phase.q' * system.A * phase.q / 2

@doc raw"""
    estimatorKinetic(phase, system)

Estimates the kinetic energy (using the virial estimator) for an isomorphic classical system under a phase.

``\hat K_{\rm v}(\boldsymbol x)=\frac 1{2\beta}+\frac1{2n}(\bar x-\boldsymbol x)\partial{\varphi(x)}{\boldsymbol x}``
"""
const estimatorKineticVirial = (phase::HamiltonianPhase, system::IsomorphicClassicalSystem) -> begin
    x = system.S * phase.q
    centroid = sum(x) / length(x)
    1 / 2system.β + (centroid .- x)' * system.quantum.force.(x) / system.n / 2.
end

@doc raw"""
    estimatorHessian(phase, system)

Estimates the hessian for an isomorphic classical system under a phase.

``\hat H(\boldsymbol x)=\frac1n\sum_iV''(x_i)``
"""
const estimatorHessian = (phase::HamiltonianPhase, system::IsomorphicClassicalSystem) -> sum(system.quantum.hessian, system.S * phase.q) / system.n
