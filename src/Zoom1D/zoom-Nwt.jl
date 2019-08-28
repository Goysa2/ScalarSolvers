export zoom_Nwt

function zoom_Nwt(h :: AbstractNLPModel,
                  nlpstop :: AbstractStopping;
                  verbose :: Bool = true,
				  kwargs...)

	optimal, nlpstop = Trouve_IntervalleA(h, nlpstop; direction = :Nwt, verbose = verbose, kwargs...)
	return optimal, nlpstop
end
