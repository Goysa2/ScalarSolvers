export TR_Sec
function TR_Sec(h :: AbstractLineFunction,
               t₀ :: Float64,
               tₘ :: Float64;
               tol :: Float64=1e-7,
               maxiter :: Int=50,
               verbose :: Bool=true)
    #print("on entre dans ")
    # Trust region parameters
    eps1 = 0.2
    eps2 = 0.8
    red = 0.5
    aug = 2
    Δ = 1.0
    t = t₀

    iter = 0;
    hₖ = obj(h,t)
    gₖ = grad(h, t)
    secₖ = 1.0

    q(d) = hₖ + gₖ*d + 0.5*secₖ*d^2

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)

        dS = -gₖ/secₖ; # direction de secante

        if q(Δ)<q(-Δ)
            d=Δ
        else
            d=-Δ
        end
        if (abs(dS)<Δ) & (q(d)>q(dS))
            d=dS
        end

        # Numerical reduction computation
        htestTR = obj(h,t+d)
        gtestTR = grad(h,t+d)

        pred = gₖ*d + 0.5*secₖ*d^2

        # truc numérique pour éviter les erreurs d'arrondi lorsque la réduction est faible
        if pred > - 1e-10
            ared = (gₖ + gtestTR)*d/2
        else
            ared = htestTR-hₖ
        end

        ratio = ared / pred
        if (ratio < eps1) & (abs(d)==Δ)
            Δ = red*Δ
        else
            # two points memory
            tpred = t
            gₖₘ₁ = gₖ

            t = t + d
            hₖ = htestTR
            gₖ = gtestTR

            s = t-tpred
            y = gₖ - gₖₘ₁
            secₖ = y/s

            if ratio > eps2
                Δ = aug * Δ
            end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    nftot=h.nlp.counters.neval_obj+h.nlp.counters.neval_grad+h.nlp.counters.neval_hprod
    ht=obj(h,t)
    return (t, true,ht,nftot)
end
