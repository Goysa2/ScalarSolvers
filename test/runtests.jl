using ScalarSolvers
using Test

using JuMP, NLPModels, NLPModelsJuMP
using State, Stopping
using OptimizationProblems

stop_solvers = [:bissect_secA_stop, :bissect_sec_stop, :bissect_nwt_stop, :bissect_stop, :ARC_SecA_Stop, :ARC_Sec_Stop, :ARC_Nwt_Stop, :TR_SecA_Stop, :TR_Sec_Stop, :TR_Nwt_Stop]

for solver in stop_solvers
    printstyled("we solve AMPGO02 using $(String(solver)) \n", color = :green)
    nlp2 = MathProgNLPModel(AMPGO02());
    nlp_at_x = NLPAtX([3.7])
    nlp_stop = NLPStopping(nlp2, Stopping.unconstrained, nlp_at_x)
    optimal, stop = eval(solver)(nlp2, nlp_stop, verbose = true)
    status(stop) == :Optimal ? color_status = :light_green : color_status = :light_yellow
    printstyled("status = $(String(status(stop))) \n", color = color_status)
    @test status(stop) == :Optimal
end
