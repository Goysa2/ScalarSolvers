export zoom_Sec

function zoom_Sec(h :: LineModel,
				  t₀ :: Float64,
				  tₘ :: Float64;
				  kwargs...)

	(topt, iter) = trouve_intervalleA(h, t₀, tₘ; direction = "Sec", kwargs...)
	return (topt, iter)
end
