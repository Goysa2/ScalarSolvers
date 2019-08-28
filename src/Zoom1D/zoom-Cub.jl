export zoom_Cub

function zoom_Cub(h :: AbstractNLPModel,
                  nlpstop :: AbstractStopping;
                  verbose :: Bool = true,
				  kwargs...)

	optimal, nlpstop = Trouve_IntervalleA(h, nlpstop; direction = :Cub, verbose = verbose, kwargs...)
	return optimal, nlpstop
end
