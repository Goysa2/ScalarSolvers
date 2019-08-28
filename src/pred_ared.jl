export pred_ared_computation

function pred_ared_computation(gₖ, fₖ, dersec, d, ftestTR, gtestTR;
                               seuil = -1e-10)
  pred = gₖ .* d + 0.5 .* dersec .* d .^ 2

  if true in (pred .> seuil)
    ared = (gₖ .+ gtestTR) .* d ./ 2.0
  else
    ared = ftestTR .- fₖ
  end

  ratio = ared ./ pred

  return (pred[1], ared[1], ratio[1])
end
