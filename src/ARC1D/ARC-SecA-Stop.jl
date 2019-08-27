export ARC_SecA_Stop

function ARC_SecA_Stop(h :: AbstractNLPModel,
					  nlpstop :: AbstractStopping;
				 	  kwargs...)

	OK, stop = ARC_generic_stop(h, nlpstop; direction = :SecA, kwargs...)
	return OK, stop
end
