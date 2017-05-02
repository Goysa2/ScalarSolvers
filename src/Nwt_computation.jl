export Nwt_computation
function Nwt_computation(t::Float64,
                         d::Float64,
                         gtestTR::Float64,
                         ftestTR::Float64,
                         h::AbstractLineFunction)

  t=t+d
  gₖ=gtestTR
  fₖ=ftestTR
  hₖ=hess(h,t)

  return (t,gₖ,fₖ,hₖ)

end
