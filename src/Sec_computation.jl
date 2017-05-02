export Sec_computation
function Sec_computation(t::Float64,
                         gₖ::Float64,
                         d::Float64,
                         ftestTR::Float64,
                         gtestTR::Float64)
  tpred = t
  gₖₘ₁ = gₖ

  t = t + d
  fₖ = ftestTR
  gₖ = gtestTR

  s = t-tpred
  y = gₖ - gₖₘ₁
  secₖ = y/s

  return (t,fₖ,gₖ,s,y,secₖ)
end
