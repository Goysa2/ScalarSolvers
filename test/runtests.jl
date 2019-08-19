using ScalarSolvers
using Test

using JuMP, NLPModels, NLPModelsJuMP
using State, Stopping
using OptimizationProblems

nlp1 = MathProgNLPModel(AMPGO02());
TR_Nwt(nlp1, verbose = true)

nlp2 = MathProgNLPModel(AMPGO02());
nlp_at_x = NLPAtX([2.7])
stop_meta = StoppingMeta(max_iter = 5)
nlp_stop = NLPStopping(nlp2, Stopping.unconstrained, nlp_at_x)
TR_Nwt_Stop(nlp2, nlp_stop, verbose = true)

# for solver in scalar_solvers
#   println(solver)
#   (topt, iter) = solver(nlp, verbose = true)
#   nftot = nlp.counters.neval_obj +
#           nlp.counters.neval_grad +
#           nlp.counters.neval_hess
#   println("Total functions and derivatives (and hessian if necessary): ",
#            nftot, "  iterations: ", iter)
#   println("(topt, iter)=", (topt, iter))
#   reset!(nlp)
# end
