module ScalarSolvers

using JuMP
using NLPModels
using Optimize, LSDescentMethods
using PolynomialRoots #Polynomials

include("includes.jl")

include("solvers.jl")

end # module
