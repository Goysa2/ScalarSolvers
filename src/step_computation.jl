export step_computation

function step_computation(direction :: Symbol,
                          h :: AbstractNLPModel,
                          t, d, gₖ, fₖ, ftestTR, gtestTR)

  T = eltype(t)

  if direction == :Nwt
    t = t + d; fₖ = ftestTR; gₖ = gtestTR; H = hess(h, t)
    return (t, fₖ, gₖ, H)
  end
  if direction == :Sec
	tpred = t; gₖₘ₁ = gₖ; t = t + d; fₖ = ftestTR
	gₖ = gtestTR; s = t - tpred; y = gₖ - gₖₘ₁; secₖ = y / s
	return (t, fₖ, gₖ, secₖ)
  end
  if direction == :SecA
	tpred = t; gkm1 = gₖ; fkm1 = fₖ
	t = t .+ d; gₖ = gtestTR; fₖ = ftestTR
 	s = t .- tpred; y = gₖ .- gkm1
	Γ = 3 .* (gₖ .+ gkm1) .* s .- 6 .* (fₖ .- fkm1)
	if (y .* s .+ Γ) < eps(T) .* (s .^ 2) 
	  yt = y
	else
	  yt = y .+ Γ ./ s
	end
	secₖ = yt ./ s
  	return (t, fₖ, gₖ, secₖ)
  end
end
