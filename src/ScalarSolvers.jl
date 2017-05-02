module ScalarSolvers

using JuMP
using NLPModels
using Optimize
#using PolynomialRoots
using Polynomials

include("includes.jl")

include("solvers.jl")

end # module
