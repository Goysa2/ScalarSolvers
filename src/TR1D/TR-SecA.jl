export TR_SecA

function TR_SecA(h :: LineModel,
				 t₀ :: Float64,
				 tₘ :: Float64;
				 kwargs...)

	(t, iter) = TR_generic(h, t₀, tₘ; direction = "SecA", kwargs...)
	return (t, iter)
end
