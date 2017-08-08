export zoom_Nwt

function zoom_Nwt(h :: LineModel,
				  t₀ :: Float64,
				  tₘ :: Float64;
				  kwargs...)

	(topt, iter) = trouve_intervalleA(h, t₀, tₘ; direction = "Nwt", kwargs...)
	return (topt, iter)
end
