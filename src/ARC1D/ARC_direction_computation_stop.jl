export ARC_step_computation_stop

# function ARC_step_computation_stop(h::Float64,
#                               g::Float64,
#                               Δ::Float64;
#                               kwargs...)
function ARC_step_computation_stop(h, g, Δ; kwargs...)

    discr1 = (h) .^ 2 .- 4 .* (g ./ Δ)
    discr2 = (h) .^ 2 .+ 4 .* (g ./ Δ)

    if h .> .0
      if g > .0
        d = (-2 .* g) ./ (h .+ sqrt.(discr2))
      else
        d = (-2 .* g) ./ (h .+ sqrt.(discr1))
      end
    else
      if g < .0
        d = (-h .+ sqrt.(discr1)) ./ (2 . / Δ)
      else
        d = (-h .+ sqrt.(discr2)) ./ (-2 ./ Δ)
      end
    end

    return d
  end
