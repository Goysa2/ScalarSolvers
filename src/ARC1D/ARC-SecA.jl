export ARC_SecA

function ARC_SecA(h :: AbstractNLPModel,
					  nlpstop :: AbstractStopping;
				 	  kwargs...)

	OK, stop = ARC_generic(h, nlpstop; direction = :SecA, kwargs...)
	return OK, stop
end
