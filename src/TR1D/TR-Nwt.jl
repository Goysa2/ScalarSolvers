export TR_Nwt

function TR_Nwt(h :: LineModel,
				t₀ :: Float64,
				tₘ :: Float64;
				kwargs...)

	(t, iter) = TR_generic(h, t₀, tₘ; direction = "Nwt", kwargs...)
	return (t, iter)
end
