export TR_Nwt

function TR_Nwt(h :: AbstractNLPModel;
				kwargs...)

	(t, f, opt_res, iter, optimality, tired, status, hf, hg,hh) = TR_generic(h; direction = "Nwt", kwargs...)
	return (t, f, opt_res, iter, optimality, tired, status, hf, hg,hh)
end
