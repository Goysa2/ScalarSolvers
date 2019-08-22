export TR_Sec_Stop

function TR_Sec_Stop(h :: AbstractNLPModel,
					 nlpstop :: AbstractStopping;
					 kwargs...)


	optimal, stop = TR_generic_Stop(h, nlpstop; direction = :Sec, kwargs...)
	return optimal, stop
end
