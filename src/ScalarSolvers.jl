module ScalarSolvers

using JuMP
using NLPModels
using Optimize, LSDescentMethods
using Polynomials

include("includes.jl")

include("solvers.jl")

end # module
