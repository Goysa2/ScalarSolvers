export SecA_computation
function SecA_computation(t::Float64,
                         gₖ::Float64,
                         fₖ::Float64,
                         d::Float64,
                         ftestTR::Float64,
                         gtestTR::Float64)
  tpred=t
  gkm1=gₖ
  fkm1=fₖ

  t=t+d
  gₖ=gtestTR
  fₖ=ftestTR

  #secant improved two points interpolation
  s=t-tpred
  y=gₖ-gkm1

  Γ=3*(gₖ+gkm1)*s-6*(fₖ-fkm1)
  if (y*s+Γ)<eps(Float64)*(s^2) #correction trop petite
    yt=y
  else
    yt=y+Γ/s
  end

  secₖ=yt/s

  return (t,fₖ,gₖ,s,y,secₖ)
end
