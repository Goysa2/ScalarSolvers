export new_ARC_generic
function new_ARC_generic(h :: AbstractLineFunction2,
                 t₀ :: Float64,
                 tₘ :: Float64;
                 tol :: Float64=1e-7,
                 maxiter :: Int=50,
                 verbose :: Bool=true,
                 eps1 :: Float64=0.25,
                 eps2 :: Float64=0.75,
                 red :: Float64=0.2,
                 aug :: Float64=5.0,
                 Δ :: Float64=0.5,
                 direction :: String="Nwt")
    t = t₀

    iter = 0;
    fₖ = obj(h,t)
    gₖ = grad(h, t)

    if direction=="Nwt"
      hₖ=hess(h,t)
    elseif direction=="Sec" || direction=="SecA"
      secₖ=1.0
    end

    #q(d) = fₖ + gₖ*d + 0.5*secₖ*d^2 + (1/3*(α))*abs(d)^3

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)

      if direction=="Nwt"
        d=ARC_step_computation(hₖ,gₖ,Δ)
      elseif direction=="Sec" || direction=="SecA"
        d=ARC_step_computation(secₖ,gₖ,Δ)
        #println("d=",d)
      end

        # Numerical reduction computation
        ftestTR = obj(h,t+d)
        gtestTR = grad(h,t+d)

        if direction=="Nwt"
          (pred,ared,ratio)=pred_ared_computation(gₖ,fₖ,hₖ,d,ftestTR,gtestTR)
        elseif direction=="Sec" || direction=="SecA"
          (pred,ared,ratio)=pred_ared_computation(gₖ,fₖ,secₖ,d,ftestTR,gtestTR)
        end

        #println("(pred,ared,ratio)=(",pred,",",ared,",",ratio,")")

        if (ratio < eps1)
          #print_with_color(:red,"ratio<eps1 \n")
            Δ = red*Δ
            #println("on réduit l'intervalle")
        else
            #println("on est dans le else")
            if direction=="Nwt"
              (t,gₖ,fₖ,hₖ)=Nwt_computation(t,d,gtestTR,ftestTR,h)
            elseif direction=="Sec"
              (t,fₖ,gₖ,s,y,secₖ)  =Sec_computation(t,gₖ,d,ftestTR,gtestTR)
            elseif direction=="SecA"
              (t,fₖ,gₖ,s,y,secₖ)  =SecA_computation(t,gₖ,fₖ,d,ftestTR,gtestTR)
              #println("apres SecA_computation (t,fₖ,gₖ,s,y,secₖ)=(",t,",",fₖ,",",gₖ,",",s,",",y,",",secₖ,")")
            end

            if ratio > eps2
                Δ = aug * Δ
            end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
