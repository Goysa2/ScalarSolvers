export zoom_Sec

function zoom_Sec(h :: AbstractNLPModel,
                  nlpstop :: AbstractStopping;
                  verbose :: Bool = true,
				  kwargs...)

	optimal, nlpstop = Trouve_IntervalleA(h, nlpstop; direction = :Sec, verbose = verbose, kwargs...)
	return optimal, nlpstop
end
