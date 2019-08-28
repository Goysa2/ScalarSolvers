export TR_Sec

function TR_Sec(h :: AbstractNLPModel,
				nlpstop :: AbstractStopping;
				kwargs...)


	optimal, stop = TR_generic(h, nlpstop; direction = :Sec, kwargs...)
	return optimal, stop
end
