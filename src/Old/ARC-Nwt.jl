export ARC_Nwt

function ARC_Nwt(h :: AbstractNLPModel;
				 kwargs...)

	(t, f, opt_res, iter, optimality, tired, status, hf, hg,hh) = ARC_generic(h; direction = "Nwt", kwargs...)
	return (t, f, opt_res, iter, optimality, tired, status, hf, hg,hh)
end
