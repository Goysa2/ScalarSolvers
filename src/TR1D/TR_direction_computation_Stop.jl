export TR_step_computation_stop

# function TR_step_computation(h :: Any,
#                              # h :: Union{T, Array{T,1}},
#                              g :: Union{T, Array{T,1}},
#                              dN :: Union{T, Array{T,1}},
#                              Δ :: Union{T, Array{T,1}}) where T

function TR_step_computation_stop(h  :: Any,
                             g  :: Any,
                             dN :: Any,
                             Δ  :: Any)
  # printstyled("on est dans TR_step_computation \n", color = :bold)
  if true in (h .> 0.0)
    # @show dN
    # @show Δ
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
  # d = [d]
  # @show d
  # @show [d]
  # printstyled("on sort de TR_step_computation \n", color = :bold)

  return [d]
end
