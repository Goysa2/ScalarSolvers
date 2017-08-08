export zoom_Cub

function zoom_Cub(h :: LineModel,
				  t₀ :: Float64,
				  tₘ :: Float64;
				  kwargs...)

	(topt, iter) = trouve_intervalleA(h, t₀, tₘ; direction = "Cub", kwargs...)
	return (topt, iter)
end
