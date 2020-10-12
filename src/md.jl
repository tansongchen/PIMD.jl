using LinearAlgebra
using Statistics

abstract type ClassicalSystem end
abstract type Propagator end
abstract type Phase end

function kinetic end
function potential end

struct HamiltonianPhase <: Phase
    q::Vector
    p::Vector
end

"""
    Integrator(propagators...)

Constructs an integrator from a series of propagator, so that within a time step these propagators act sequentially on the phase.
"""
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

"""
    md(initPhase, system, integrator!, timetable, properties)

Launches an MD simulation of `system` with the initial conditions `initPhase` and collect `properties`. During the simulation, `integrator!` is used to integrate the equations of motion. With `timetable = (equilibrationSteps, sampleSteps)`, the system is first equilibrated for `equilibrationSteps` and then sampled for `sampleSteps`.
"""
md(initPhase::Phase, system::ClassicalSystem, integrator!::Integrator, timetable::Tuple{<:Integer, <:Integer}, properties::Vector{<:Function}) = begin
    phase = deepcopy(initPhase)
    nEquil, nSample = timetable
    for i in 1:nEquil
        integrator!(phase, system)
    end
    data = Dict(property => [] for property in properties)
    for i in 1:nSample
        integrator!(phase, system)
        for property in properties
            push!(data[property], property(phase, system))
        end
    end
    Dict(property => mean(data[property]) for property in properties)
end

statistics(properties::Vector{<:Function}, data::Vector{<:Dict}) = begin
    Dict(
        property => (
            mean=mean([item[property] for item in data]),
            std=std([item[property] for item in data]) / √(length(data))
        )
        for property in properties
    )
end
