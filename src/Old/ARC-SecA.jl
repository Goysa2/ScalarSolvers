export ARC_SecA

function ARC_SecA(h :: AbstractNLPModel;
				  kwargs...)

	(t, f, opt_res, iter, optimality, tired, status, hf, hg,hh) = ARC_generic(h; direction = "SecA", kwargs...)
	return (t, f, opt_res, iter, optimality, tired, status, hf, hg,hh)
end
