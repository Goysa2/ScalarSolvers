export TR_step_computation

function TR_step_computation(h  :: Any,
                             g  :: Any,
                             dN :: Any,
                             Δ  :: Any)
  if true in (h .> 0.0)
    if true in (g .> 0.0)
      d = max.(-Δ, dN)
    else
      d = min.(dN, Δ)
    end
  else
    if true in (g .> 0.0)
      d = .- Δ
    else
      d = Δ
    end
  end

  d = d[1]

  return [d]
end
