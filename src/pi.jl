using ForwardDiff: derivative, gradient

"""
    OneDimensionalQuantumSystem(β, m, potential)

Creates an 1D quantum system with temperature β and Hamiltonian ``H=p^2/2m+V(x)``.
"""
struct OneDimensionalQuantumSystem
    β::Real
    m::Real
    potential::Function
    force::Function
    hessian::Function
    OneDimensionalQuantumSystem(β::Real, m::Real, potential::Function) = begin
        force = x -> -derivative(potential, x)
        hessian = x -> -derivative(force, x)
        new(β, m, potential, force, hessian)
    end
end

const harmonicPotential = x::Real -> .5x^2
const anharmonicPotential = x::Real -> .5x^2 + .1x^3 + .01x^4
const quarticPotential = x::Real -> .25x^4

const stiffnessMatrix = (system::OneDimensionalQuantumSystem, n::Integer) -> begin
    tmp = zeros(n, n)
    for i in 1:n
        tmp[i,i] += 1
        tmp[mod1(i+1, n), mod1(i+1, n)] += 1
        tmp[i, mod1(i+1, n)] -= 1
        tmp[mod1(i+1, n), i] -= 1
    end
    Symmetric(tmp * system.m * n / system.β^2)
end
