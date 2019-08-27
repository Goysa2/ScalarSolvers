export ARC_Cub_stop
function ARC_Cub_stop(h       :: AbstractNLPModel,
                      nlpstop :: AbstractStopping;
                      t₀ = h.meta.x0, tₘ = 100.0,
                      tol = 1e-7, verbose = false)
    # ARC parameters
    eps1 = 0.2; eps2 = 0.8
    red = 0.5; aug = 20.0
    Δ = [1.0]
    t = t₀

    iter = 0;
    fₖ = obj(h, t)
    gₖ = grad(h, t)

    secₖ = [1.0]
    dN = -gₖ #pour la première iter
    A = [0.5]; B = [0.0]
    q(d) = fₖ .+ gₖ .* d .+ A .* d .^ 2 .+ B .* d .^ 3

    OK  = update_and_start!(nlpstop, x = t, fx = fₖ, gx = gₖ, g0 = copy(gₖ))

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter, t[1], gₖ[1], Δ[1])

    while !OK

      #Step computation
      Quad(t) = fₖ .+ gₖ .* t .+ A .* t .^ 2 .+ B .* t .^ 3 .+ (1 ./ (4.0 .* Δ)) .* t .^ 4
      dQuad(t) = gₖ .+ 2.0 .* A .* t .+ 3.0 .* B .* t .^ 2 .+ (1 ./ Δ) .* t .^ 3


      p = Poly([gₖ[1], 2.0 .* A[1], 3.0 .* B[1], 1.0 ./ Δ[1]])
      dR = roots(p)
      vmin = Inf; dN = 0
      for i=1:length(dR)
        rr = dR[i]
        if isreal(rr)
          rrr = real(rr)
          vact = fₖ .+ gₖ .* rrr .+ A .* rrr .^2 .+ B .* rrr .^ 3 .+ (1 ./ (4.0 .* Δ)) .* rrr .^ 4
          if true in (rrr .* gₖ .< 0.0)
            if true in (vact .< vmin)
              dN = rrr
              vmin = vact
            end
          end
        end
      end

      d = dN

      if true in (d .== 0.0)
        return OK, nlpstop
      end

      ftestTR=obj(h, t .+ d)
      gtestTR=grad(h, t .+ d)

      pred = gₖ .* d .+ A .* d .^ 2 .+ B .* d .^ 3
      if true in (pred .> -1e-10)
        ared = (gₖ .+ gtestTR) .* d ./ 2.0
      else
        ared = ftestTR .- fₖ
      end

      ratio = ared ./ pred
      if true in (ratio .< eps1)
        Δ = red .* Δ
      else
        #two point memory
        tpred = t
        gkm1 = gₖ
        fkm1 = fₖ

        t = t .+ d
        gₖ = gtestTR
        fₖ = ftestTR

        #Cubic two points interpolation
        s = t .- tpred
        y = gₖ .- gkm1
        α = -s
        z = gₖ .+ gkm1 .+ 3.0 .* ( fₖ .- fkm1) ./ α
        discr= z .^ 2 .- gₖ .* gkm1
        denom = gₖ .+ gkm1 .+ 2.0 .* z
        B = (1 ./ 3) .* (gₖ .+ gkm1 .+ 2.0 .* z) ./ ( α .* α)
        A = -(gₖ .+ z) ./ α
        if true in (ratio .> eps2)
          Δ = aug .* Δ
        end
      end
      OK  = update_and_start!(nlpstop, x = t, fx = fₖ, gx = gₖ, g0 = copy(gₖ));
      iter += 1
      verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter, t[1], gₖ[1], Δ[1], pred[1], ared[1])
    end
    optimal = OK
    return optimal, nlpstop

end
