export new_ARC_Nwt
function new_ARC_Nwt(hh :: AbstractLineFunction,
                 t₀ :: Float64,
                 tₘ :: Float64;
                 tol :: Float64=1e-7,
                 maxiter :: Int=50,
                 verbose :: Bool=true,
                 eps1 :: Float64=0.1,
                 eps2 :: Float64=0.75,
                 red :: Float64=0.1,
                 aug :: Float64=5.0)

    # Trust region parameters
    # eps1 = 0.1
    # eps2 = 0.75
    # red = 0.1
    # aug = 5.0
    Δ = 1.0
    t = t₀

    iter = 0;
    succ=0
    unsucc=0
    fₖ = obj(hh,t)
    gₖ = grad(hh,t)
    hₖ = hess(hh,t)

    #q(d) = fₖ + gₖ*d + 0.5*hₖ*d^2 + (1/3*Δ)*abs(d)^3

    #println("a l'iteration:",iter,"on a f=",fₖ," g=",gₖ," h=",hₖ)

    verbose && @printf(" iter  t         gₖ          Δ \n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)
        # #step computation
        # discr1= (hₖ)^2-4(gₖ/Δ)
        # println("discr1=",discr1)
        # discr2= (hₖ)^2+4(gₖ/Δ)
        # println("discr2=",discr2)
        # if hₖ>0
        #   #println("hₖ>0")
        #   if gₖ>0
        #     #println("gₖ>0")
        #     #print_with_color(:yellow,"CAS 1 \n")
        #     d=(-hₖ+sqrt(discr2))/(-2/Δ)
        #   else
        #     #println("gₖ<0")
        #     #print_with_color(:yellow,"CAS 2 \n")
        #     d=(-hₖ+sqrt(discr1))/(2/Δ)
        #   end
        # else
        #   #println("hₖ<0")
        #   d₁=(-hₖ+sqrt(discr2))/(-2/Δ)
        #   d₂=(-hₖ+sqrt(discr1))/(2/Δ)
        #   if q(d₁)<q(d₂)
        #     #println("q(d₁)")
        #     #print_with_color(:yellow,"CAS 3 \n")
        #     d=d₁
        #   else
        #     #print_with_color(:yellow,"CAS 4 \n")
        #     #println("q(d₂)")
        #     d=d₂
        #   end
        # end

        # if gₖ>0.0
        #   #println("gₖ>0")
        #   d=(-hₖ+sqrt(discr2))/(-2/Δ)
        # else
        #   #println("gₖ<0")
        #   d=(-hₖ+sqrt(discr1))/(2/Δ)
        # end

        #println("d=",d)

        d=ARC_step_computation(hₖ,gₖ,Δ)

        #Numerical reduction computation
        ftestTR=obj(hh,t+d)
        gtestTR=grad(hh,t+d)
        pred= gₖ*d + 0.5*hₖ*d^2
        if pred>-1e-10
          ared = (gₖ + gtestTR)*d/2
        else
          ared = ftestTR-fₖ
        end

        ratio=ared/pred
        if ratio<eps1
          Δ=red*Δ

          unsucc=unsucc+1
        else
          t=t+d

          gₖ=gtestTR
          fₖ=ftestTR
          hₖ=hess(hh,t)

          if ratio>eps2
            Δ=aug*Δ
          end
          succ=succ+1
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e\n", iter,t,gₖ,Δ)
    end

    return (t, iter)
end
