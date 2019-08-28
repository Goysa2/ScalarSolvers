export ARC_Sec

function ARC_Sec(h :: AbstractNLPModel,
					  nlpstop :: AbstractStopping;
				 	  kwargs...)

	OK, stop = ARC_generic(h, nlpstop; direction = :Sec, kwargs...)
	return OK, stop
end
