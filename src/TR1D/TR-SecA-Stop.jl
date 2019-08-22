export TR_SecA_Stop

function TR_SecA_Stop(h :: AbstractNLPModel,
					 nlpstop :: AbstractStopping;
					 kwargs...)


	optimal, stop = TR_generic_Stop(h, nlpstop; direction = :SecA, kwargs...)
	return optimal, stop
end
