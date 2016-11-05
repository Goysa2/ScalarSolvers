using ScalarSolvers
using Base.Test

using JuMP
using NLPModels
using Optimize
using Roots

include("AMPGO02.jl")

nbSolver=0
nlp=MathProgNLPModel(AMPGO02())
h=C1LineFunction(nlp,[0.0],[1.0])
Ivl=[2.7 7.5]

for solver in all_solvers
  nbSolver += 1
  println(nbSolver,"  ",solver)
  (topt,iter)=solver(h,Ivl[1],verbose=false)
  nftot = nlp.counters.neval_obj + nlp.counters.neval_grad
  println("Total functions and derivatives: ",nftot, "  iterations: ",iter)
  println("(topt, iter)=",(topt, iter))
  reset!(nlp)
end
