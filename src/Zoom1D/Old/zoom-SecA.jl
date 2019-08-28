export zoom_SecA

function zoom_SecA(h :: AbstractNLPModel;
				  kwargs...)

	(t, fₖ, normg, iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess) = trouve_intervalleA(h; direction = "SecA", kwargs...)
	return (t, fₖ, normg, iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
end
