export zoom_Cub

function zoom_Cub(h :: AbstractNLPModel;
				  kwargs...)
	(t, fₖ, normg, iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess) = trouve_intervalleA(h; direction = "Cub", kwargs...)
	return (t, fₖ, normg, iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
end
