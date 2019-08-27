export ARC_Sec_Stop

function ARC_Sec_Stop(h :: AbstractNLPModel,
					  nlpstop :: AbstractStopping;
				 	  kwargs...)

	OK, stop = ARC_generic_stop(h, nlpstop; direction = :Sec, kwargs...)
	return OK, stop
end
