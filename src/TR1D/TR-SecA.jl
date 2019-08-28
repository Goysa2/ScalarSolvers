export TR_SecA

function TR_SecA(h :: AbstractNLPModel,
				 nlpstop :: AbstractStopping;
				 kwargs...)


	optimal, stop = TR_generic(h, nlpstop; direction = :SecA, kwargs...)
	return optimal, stop
end
