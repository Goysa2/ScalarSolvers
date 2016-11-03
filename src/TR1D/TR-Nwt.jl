function TR_Nwt(hh :: C2LineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true)

    # Trust region parameters
    eps1 = 0.2
    eps2 = 0.8
    red = 0.5
    aug = 2
    Δ = 1.0
    t = t₀

    iter = 0;
    hₖ = obj(hh,t)
    gₖ = grad(hh, t)
    kₖ = hess(hh,t)

    q(d) = hₖ + gₖ*d + 0.5*kₖ*d^2

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) || (iter == 0)

        dN = -gₖ/kₖ; # direction de Newton

        if q(Δ)<q(-Δ)
            d=Δ
        else
            d=-Δ
        end
        if (abs(dN)<Δ) & (q(d)>q(dN))
            d=dN
        end

        # Numerical reduction computation
        htestTR = obj(hh,t+d)
        gtestTR = grad(hh,t+d)

        pred = gₖ*d + 0.5*kₖ*d^2

        # truc numérique pour éviter les erreurs d'arrondi lorsque la réduction est faible
        if pred > - 1e-10
            ared = (gₖ + gtestTR)*d/2
        else
            ared = htestTR-hₖ
        end

        ratio = ared / pred
        if (ratio < eps1)
            Δ = red*Δ
        else
            t=t+d

            gₖ=gtestTR
            hₖ=htestTR
            kₖ=hess(hh,t)

            if ratio > eps2
                Δ = aug * Δ
            end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
