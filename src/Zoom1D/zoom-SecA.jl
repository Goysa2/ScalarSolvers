export zoom_SecA

function zoom_SecA(h :: LineModel,
				  t₀ :: Float64,
				  tₘ :: Float64;
				  kwargs...)

	(topt, iter) = trouve_intervalleA(h, t₀, tₘ; direction = "SecA", kwargs...)
	return (topt, iter)
end
