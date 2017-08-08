export ARC_Sec

function ARC_Sec(h :: LineModel,
		     	 t₀ :: Float64,
				 tₘ :: Float64;
				 kwargs...)

	(t, iter) = ARC_generic(h, t₀, tₘ; direction = "Sec", kwargs...)
	return (t, iter)
end
