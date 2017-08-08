export ARC_SecA

function ARC_SecA(h :: LineModel,
		     	  t₀ :: Float64,
				  tₘ :: Float64;
				  kwargs...)

	(t, iter) = ARC_generic(h, t₀, tₘ; direction = "SecA", kwargs...)
	return (t, iter)
end
