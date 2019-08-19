export pred_ared_computation_stop
# function pred_ared_computation(gₖ :: Float64,
#                                fₖ :: Float64,
#                                dersec :: Float64,
#                                d :: Float64,
#                                ftestTR :: Float64,
#                                gtestTR :: Float64;
#                                seuil :: Float64=-1e-10)

function pred_ared_computation_stop(gₖ, fₖ, dersec, d, ftestTR, gtestTR; seuil :: Float64=-1e-10)
  pred = gₖ .* d + 0.5 .* dersec .* d .^ 2
  # @show pred
  # @show pred .> seuil

  if true in (pred .> seuil)
    ared = (gₖ .+ gtestTR) .* d ./ 2.0
  else
    #print_with_color(:green,"!!")
    ared = ftestTR .- fₖ
  end

  ratio = ared ./ pred
  # @show ratio

  return (pred[1], ared[1], ratio[1])
end
