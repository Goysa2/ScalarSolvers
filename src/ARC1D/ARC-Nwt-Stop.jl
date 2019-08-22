export ARC_Nwt_Stop

function ARC_Nwt_Stop(h :: AbstractNLPModel,
					  nlpstop :: AbstractStopping;
				 	  kwargs...)

	(t, f, opt_res, iter, optimality, tired, status, hf, hg,hh) = ARC_generic_stop(h, nlpstop; direction = :Nwt, kwargs...)
	return (t, f, opt_res, iter, optimality, tired, status, hf, hg,hh)
end
