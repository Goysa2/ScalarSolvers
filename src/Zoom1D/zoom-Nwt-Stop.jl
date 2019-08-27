export zoom_Nwt_stop

function zoom_Nwt_stop(h :: AbstractNLPModel,
                       nlpstop :: AbstractStopping;
                       verbose :: Bool = true,
				  	   kwargs...)

	optimal, nlpstop = trouve_intervalleA_stop(h, nlpstop; direction = :Nwt, kwargs...)
	return optimale, nlpstop
end
