export ARC_Sec

function ARC_Sec(h :: AbstractNLPModel;
				 kwargs...)

	(t, f, opt_res, iter, optimality, tired, status, hf, hg,hh) = ARC_generic(h; direction = "Sec", kwargs...)
	return (t, f, opt_res, iter, optimality, tired, status, hf, hg,hh)
end
