export TR_Nwt_Stop

function TR_Nwt_Stop(h :: AbstractNLPModel,
					 nlpstop :: AbstractStopping;
					 kwargs...)

					 
	(t, f, opt_res, iter, optimality, tired, status, hf, hg,hh) = TR_generic_Stop(h, nlpstop; direction = :Nwt, kwargs...)
	return (t, f, opt_res, iter, optimality, tired, status, hf, hg,hh)
end
