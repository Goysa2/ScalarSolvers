using ScalarSolvers
using Test

using JuMP, NLPModels, NLPModelsJuMP
using State, Stopping
using OptimizationProblems

# nlp1 = MathProgNLPModel(AMPGO02());
# TR_Nwt(nlp1, verbose = true)

stop_solvers = [:ARC_Nwt_Stop, :TR_SecA_Stop, :TR_Sec_Stop, :TR_Nwt_Stop]

for solver in stop_solvers
    printstyled("we solve AMPGO02 using $(String(solver)) \n", color = :green)
    nlp2 = MathProgNLPModel(AMPGO02());
    nlp_at_x = NLPAtX([3.7])
    # stop_meta = StoppingMeta(max_iter = 5)
    nlp_stop = NLPStopping(nlp2, Stopping.unconstrained, nlp_at_x)#, meta = stop_meta)
    optimal, stop = eval(solver)(nlp2, nlp_stop, verbose = true)
    status(stop) == :Optimal ? color_status = :light_green : color_status = :light_yellow
    printstyled("status = $(String(status(stop))) \n", color = color_status)
end


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
