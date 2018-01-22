export TR_SecA

function TR_SecA(h :: AbstractNLPModel;
				 kwargs...)

	(t, f, opt_res, iter, optimality, tired, status, hf, hg,hh) = TR_generic(h; direction = "SecA", kwargs...)
	return (t, f, opt_res, iter, optimality, tired, status, hf, hg,hh)
end
