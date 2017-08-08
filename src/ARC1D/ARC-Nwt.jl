export ARC_Nwt

function ARC_Nwt(h :: LineModel,
		     	 t₀ :: Float64,
				 tₘ :: Float64;
				 kwargs...)

	(t, iter) = ARC_generic(h, t₀, tₘ; direction = "Nwt", kwargs...)
	return (t, iter)
end
