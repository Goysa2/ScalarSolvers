export ARC_Nwt_Stop

function ARC_Nwt_Stop(h :: AbstractNLPModel,
					  nlpstop :: AbstractStopping;
				 	  kwargs...)

	OK, stop = ARC_generic_stop(h, nlpstop; direction = :Nwt, kwargs...)
	return OK, stop
end
