export zoom_SecA

function zoom_SecA(h :: AbstractNLPModel,
                   nlpstop :: AbstractStopping;
                   verbose :: Bool = true,
				   kwargs...)

	optimal, nlpstop = Trouve_IntervalleA(h, nlpstop; direction = :SecA, verbose = verbose, kwargs...)
	return optimal, nlpstop
end
