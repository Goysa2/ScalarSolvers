export new_ARC_Sec
function new_ARC_Sec(h :: AbstractLineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true)

    #print("on entre dans ")
    # Trust region parameters
    eps1 = 0.1
    eps2 = 0.75
    red = 0.5
    aug = 2
    α = 0.5
    t = t₀

    iter = 0;
    hₖ = obj(h,t)
    gₖ = grad(h, t)


    secₖ = 1.0

    #q(d) = hₖ + gₖ*d + 0.5*secₖ*d^2

    verbose && @printf(" iter  t         gₖ          α        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,α)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)

        #step computation
      #   discr=secₖ^2-4*(gₖ/α)
      #   if discr<0
      #     discr=secₖ^2+4*(gₖ/α)
      #   end
      #
      # #println("on a discr")
      #
      # if seck<0
      #   dNp=(-secₖ+sqrt(discr))/(2/Δ) #direction de Newton
      # else
      #   dNp=-2*gₖ/(secₖ+sqrt(discr))
      # end
      #
      #   dNn=(-secₖ-sqrt(discr))/(2/α) #direction de Newton
      #   #println("on a dNn")
      #
      #   if q(dNp)<q(dNn)
      #       d=dNp
      #   else
      #       d=dNn
      #   end

      d=ARC_step_computation(hₖ,gₖ,α)

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
        if (ratio < eps1)
            α = red*α
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
                α = aug * α
            end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,α,pred,ared)
    end

    return (t, iter)
end
