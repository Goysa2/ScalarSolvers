export TR_Sec

function TR_Sec(h :: LineModel,
				t₀ :: Float64,
				tₘ :: Float64;
				kwargs...)

	(t, iter) = TR_generic(h, t₀, tₘ; direction = "Sec", kwargs...)
	return (t, iter)
end
