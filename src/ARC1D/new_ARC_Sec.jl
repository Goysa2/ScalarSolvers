export new_ARC_Sec
function new_ARC_Sec(h :: AbstractLineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true,
                eps1 :: Float64=0.1,
                eps2 :: Float64=0.75,
                red :: Float64=0.5,
                aug :: Float64=2.0)

    #print("on entre dans ")
    # Trust region parameters
    #eps1 = 0.25
    #eps2 = 0.75
    #red = 0.5
    #aug = 2
    α = 0.5
    t = t₀

    iter = 0;
    fₖ = obj(h,t)
    gₖ = grad(h, t)


    secₖ = 1.0

    #q(d) = fₖ + gₖ*d + 0.5*secₖ*d^2 + (1/3*(α))*abs(d)^3

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

      d=ARC_step_computation(secₖ,gₖ,α)

        # Numerical reduction computation
        ftestTR = obj(h,t+d)
        gtestTR = grad(h,t+d)

        pred = gₖ*d + 0.5*secₖ*d^2

        # truc numérique pour éviter les erreurs d'arrondi lorsque la réduction est faible
        if pred > - 1e-10
            ared = (gₖ + gtestTR)*d/2
        else
            ared = ftestTR-fₖ
        end

        ratio = ared / pred
        if (ratio < eps1)
          #print_with_color(:red,"ratio<eps1")
            α = red*α
        else
            # two points memory
            #print_with_color(:yellow,"ratio>eps1")
            tpred = t
            gₖₘ₁ = gₖ

            t = t + d
            fₖ = ftestTR
            gₖ = gtestTR

            s = t-tpred
            y = gₖ - gₖₘ₁
            secₖ = y/s

            if ratio > eps2
                α = aug * α
            end
        end

        iter += 1
        #println("secₖ=",secₖ)
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,α,pred,ared)
    end

    return (t, iter)
end
