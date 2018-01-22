export TR_Sec

function TR_Sec(h :: AbstractNLPModel;
				kwargs...)

	(t, f, opt_res, iter, optimality, tired, status, hf, hg,hh) = TR_generic(h; direction = "Sec", kwargs...)
	return (t, f, opt_res, iter, optimality, tired, status, hf, hg,hh)
end
