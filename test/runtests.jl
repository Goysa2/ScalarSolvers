using ScalarSolvers
using Test

using JuMP, NLPModels, NLPModelsJuMP
using State, Stopping
using OptimizationProblems
using Polynomials

for solver in scalar_solvers
    printstyled("we solve AMPGO02 using $solver \n", color = :green)
    nlp = MathProgNLPModel(AMPGO02());
    nlp_at_x = NLPAtX([3.7])
    nlp_stop = NLPStopping(nlp, Stopping.unconstrained, nlp_at_x)
    optimal, stop = eval(solver)(nlp, nlp_stop, verbose = true)
    @show optimal
    @show typeof(stop)
    status(stop) == :Optimal ? color_status = :light_green : color_status = :light_yellow
    printstyled("status = $(String(status(stop))) \n", color = color_status)
    @test status(stop) == :Optimal
end
