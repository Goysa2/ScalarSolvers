using ScalarSolvers
using Base.Test

using JuMP, NLPModels, Optimize, Polynomials
using OptimizationProblems

nlp = MathProgNLPModel(AMPGO02())
h = LineModel(nlp, [0.0], [1.0]);

for solver in scalar_solvers
  println(solver)
  (topt, iter) = solver(h, 2.7, 7.5, verbose = true)
  nftot = nlp.counters.neval_obj +
          nlp.counters.neval_grad +
          nlp.counters.neval_hess
  println("Total functions and derivatives (and hessian if necessary): ",
           nftot, "  iterations: ", iter)
  println("(topt, iter)=", (topt, iter))
  reset!(nlp)
end
