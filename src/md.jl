using LinearAlgebra
using Statistics

abstract type ClassicalSystem end
abstract type Propagator end
abstract type Phase end
abstract type Property end

function kinetic end
function potential end
function collect!(::Property, ::Phase, ::ClassicalSystem) end
function compute(::Property, ::ClassicalSystem) end

struct HamiltonianPhase <: Phase
    z::Vector
    q::SubArray
    p::SubArray
    HamiltonianPhase(z::Vector) = begin
        n = length(z) ÷ 2
        q = @view z[1:n]
        p = @view z[n+1:2n]
        new(z, q, p)
    end
end

struct Integrator
    propagators::Tuple{Vararg{Propagator}}
    Integrator(propagators...) = begin
        new(propagators)
    end
end
(integrator!::Integrator)(phase::Phase, system::ClassicalSystem) = begin
    for propagator! in integrator!.propagators
        propagator!(phase, system)
    end
end

md(initPhase::Phase, system::ClassicalSystem, integrator!::Integrator, timetable::Tuple{<:Integer, <:Integer}, properties::Vector{<:Property}) = begin
    phase = deepcopy(initPhase)
    nEquil, nSample = timetable
    for i in 1:nEquil
        integrator!(phase, system)
    end
    for i in 1:nSample
        integrator!(phase, system)
        for property in properties
            collect!(property, phase, system)
        end
    end
    [compute(property, system) for property in properties]
end

struct StatisticalProperty <: Property
    estimator::Function
    estimates::Vector{Real}
    StatisticalProperty(estimator) = new(estimator, [])
end

collect!(property::StatisticalProperty, phase::Phase, system::ClassicalSystem) = begin
    push!(property.estimates, property.estimator(phase, system))
end

compute(property::StatisticalProperty, ::ClassicalSystem) = begin
    mean(property.estimates)
end