export zoom_Sec

function zoom_Sec(h :: AbstractNLPModel;
				  kwargs...)

	(t, fₖ, normg, iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess) = trouve_intervalleA(h; direction = "Sec", kwargs...)
	return (t, fₖ, normg, iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
end
