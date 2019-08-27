export ARC_step_computation_stop

# function ARC_step_computation_stop(h::Float64,
#                               g::Float64,
#                               Δ::Float64;
#                               kwargs...)
function ARC_step_computation_stop(h, g, Δ; kwargs...)

    h2     = ((h) .^ 2)
    gdelta = (4 .* (g ./ Δ))
    discr1 = h2 .- gdelta
    # @show discr1
    discr2 = h2 .+ gdelta
    # @show discr2

    # @show h
    # @show g
    # @show h .> .0
    # println("ici!")
    # @show g .> .0
    # println("ici 2 !!")

    if true in (h .> .0)
      if true in (g .> .0)
        # printstyled("if 1 \n", color = :bold)
        d = (-2 .* g) ./ (h .+ sqrt.(discr2))
      else
        # printstyled("else 1 \n", color = :bold)
        d = (-2 .* g) ./ (h .+ sqrt.(discr1))
      end
    else
      if true in (g .< .0)
        # printstyled("if 2 \n", color = :bold)
        d = (-h .+ sqrt.(discr1)) ./ (2 ./ Δ)
      else
        # printstyled("else 2 \n", color = :bold)
        d = (-h .+ sqrt.(discr2)) ./ (-2 ./ Δ)
      end
    end

    return vec(d)
  end
