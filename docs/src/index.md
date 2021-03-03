# Introduction

`PIMD.jl` is a comprehensive Julia implementation of path integral molecular dynamics for 1D model systems.

For the purpose of developing new 

## Defining a quantum system

To define a one-dimensional quantum system we need temperature, mass and potential function. Pre-defined potential functions include `harmonicPotential`, `weaklyAnharmonicPotential` and `quarticPotential`. For example,

```julia
β, m, V = 1.2, 3.4, quarticPotential
quantum = OneDimensionalQuantumSystem(β, m, V)
```

## Defining a classical isomorphic system

Classical isomorphic system is a system that correspond to an one-dimensional quantum system and gives the same statistical property in the infinite-bead limit. This package supports three types of classical isomorphic systems, namely `Primitive`, `Normal` and `Staging`.

```julia
primitive = Primitive(quantum, 16)
```

## Defining a integrator

An integrator is a way to discretize the continuous-time evolution that is suitable for numerical integration. An integrator is further decomposed into a series of propagators that updates a specific part of the phase (e.g. the momentum, position or else). The simplest one is the Verlet, or second-order symplectic integrator:

```julia
integrator = begin
    B = PotentialPropagator(Δt / 2)
    A = KineticPropagator(Δt)
    Integrator(B, A, B)
end
```

## Defining properties and run simulation

We may choose several properties of interest, such as `estimatorPotential`, `estimatorKinetic`, `estimatorKineticVirial` to collect during our simulation. To run a simulation we simply need to initialize a `phase` and specify a list of `properties`, and call the main function `md`.

```julia
properties = [estimatorPotential, estimatorKinetic]
data = md(phase, primitive, integrator, (10, 100), properties)
```

`data` is a dict object from which the properties can be extracted as

```julia
potential = data[estimatorPotential]
```

---

Please also view the next page for more detailed description of APIs.
