export TR_Nwt_Stop

function TR_Nwt_Stop(h :: AbstractNLPModel,
					 nlpstop :: AbstractStopping;
					 kwargs...)

	optimal, stop = TR_generic_Stop(h, nlpstop; direction = :Nwt, kwargs...)
	return optimal, stop
end
