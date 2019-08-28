export ARC_step_computation

function ARC_step_computation(h, g, Δ; kwargs...)

    h2     = ((h) .^ 2)
    gdelta = (4 .* (g ./ Δ))
    discr1 = h2 .- gdelta
    discr2 = h2 .+ gdelta

    if true in (h .> .0)
      if true in (g .> .0)
        d = (-2 .* g) ./ (h .+ sqrt.(discr2))
      else
        d = (-2 .* g) ./ (h .+ sqrt.(discr1))
      end
    else
      if true in (g .< .0)
        d = (-h .+ sqrt.(discr1)) ./ (2 ./ Δ)
      else
        d = (-h .+ sqrt.(discr2)) ./ (-2 ./ Δ)
      end
    end

    return vec(d)
  end
