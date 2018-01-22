export ARC_Cub
function ARC_Cub(h :: AbstractNLPModel;
                 t₀ :: Float64 = -10.0,
                 tₘ :: Float64 = 100.0,
                 tol :: Float64 = 1e-7,
                 maxiter :: Int = 50,
                 verbose :: Bool = false)

    # Trust region parameters
    eps1 = 0.2
    eps2 = 0.8
    red = 0.5
    aug = 20
    Δ = 1.0
    t = t₀

    iter = 0;
    fₖ = obj(h,[t])[1]
    gₖ = grad(h, [t])[1]

    secₖ = 1.0
    dN=-gₖ #pour la première iter
    A=0.5
    B=0.0
    q(d) = fₖ + gₖ * d + A * d^2 + B * d^3

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter, t, gₖ, Δ)

    while ((abs(gₖ[1])>tol) & (iter < maxiter)) || (iter == 0)

      #Step computation
      Quad(t) = fₖ + gₖ * t + A * t^2 + B * t^3 + (1 / (4 * Δ)) * t^4
      dQuad(t) = gₖ + 2 * A * t + 3 * B * t^2 + (1 / Δ) * t^3

      #dR=PolynomialRoots.roots([gₖ,2*A,3*B,(1/Δ)])
      #TECHNIQUEMENT LE CALCUL DES RACINES PAR PolynomialRoots À LA MÊME PRÉCISION QUE SCILAB, MAIS IL DONNE UN SEG FAULT (comme pour TR_Cub)
      p = Poly([gₖ, 2 * A, 3 * B, (1 / Δ)])
      dR = roots(p)
      vmin=Inf
      dN=0
      for i=1:length(dR)
        rr = dR[i]
        if isreal(rr)
          rr = real(rr)
          vact = Quad(rr)
          if rr * gₖ[1] < 0.0
            if vact < vmin
              dN = rr
              vmin = vact
            end
          end
        end
      end

      d = dN

      if d == 0.0
        return (t, fₖ, norm(gₖ, Inf), iter + 1, false, false, :FailedDirectionComputation, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
        #L'OUTIL DE CALCUL DES RACINES DU PACKAGE Polynomials ARRONDIE À 0 CERTAINES RACINES D'OÙ LA NOUVELLE CONDITION
      end

      ftestTR=obj(h, [t + d])[1]
      gtestTR=grad(h, [t + d])[1]

      pred = gₖ[1] * d + A *d ^2 + B * d^3
      if pred > -1e-10
        ared=(gₖ[1] + gtestTR[1]) * d / 2
      else
        ared = ftestTR - fₖ
      end

      ratio = ared / pred
      if ratio < eps1
        Δ = red * Δ
      else
        #two point memory
        tpred = t
        gkm1 = gₖ
        fkm1 = fₖ

        t = t + d
        gₖ = gtestTR
        fₖ = ftestTR

        #Cubic two points interpolation
        s = t - tpred
        y = gₖ - gkm1
        α = -s
        z = gₖ + gkm1 + 3 *( fₖ - fkm1) / α
        discr= z^2 - gₖ * gkm1
        denom = gₖ + gkm1 + 2 * z
        B = (1 / 3) * (gₖ + gkm1 + 2 * z) / ( α * α)
        A = -(gₖ + z) / α
        if ratio > eps2
          Δ = aug * Δ
        end
      end
      iter += 1
      verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter, t, gₖ, Δ, pred, ared)
    end
    if maxiter <= iter
        tired = true
    else
        tired = false
    end
    if (abs(gₖ) > tol)
        optimal = false
    else
        optimal = true
    end
    status = :tmp
    return (t, fₖ, norm(gₖ, Inf), iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)

end
