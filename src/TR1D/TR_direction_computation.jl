export TR_step_computation

function TR_step_computation(h :: Float64,
                             g :: Float64,
                             dN :: Float64,
                             Δ :: Float64)

  if h > 0.0
    if g > 0.0
      d = max(-Δ, dN)
    else
      d = min(dN, Δ)
    end
  else
    if g > 0.0
      d = -Δ
    else
      d = Δ
    end
  end

  return d

end
