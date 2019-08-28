module ScalarSolvers

using JuMP
using NLPModels
using LinearAlgebra
using Printf
using State, Stopping
using Polynomials

include("includes.jl")

include("solvers.jl")

end # module
