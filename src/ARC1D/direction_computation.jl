export ARC_step_computation

function ARC_step_computation(h::Float64,
                               g::Float64,
                               Δ::Float64)

    discr1 = (h)^2-4(g/Δ)
    discr2 = (h)^2+4(g/Δ)

    if g>0.0
      #println("gₖ>0")
      d=(-h+sqrt(discr2))/(-2/Δ)
    else
      #println("gₖ<0")
      d=(-h+sqrt(discr1))/(2/Δ)
    end

    return d
  end
