export TR_Nwt

function TR_Nwt(h :: AbstractNLPModel,
				nlpstop :: AbstractStopping;
				kwargs...)

	optimal, stop = TR_generic(h, nlpstop; direction = :Nwt, kwargs...)
	return optimal, stop
end
