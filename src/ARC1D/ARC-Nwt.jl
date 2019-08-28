export ARC_Nwt

function ARC_Nwt(h :: AbstractNLPModel,
					  nlpstop :: AbstractStopping;
				 	  kwargs...)

	OK, stop = ARC_generic(h, nlpstop; direction = :Nwt, kwargs...)
	return OK, stop
end
