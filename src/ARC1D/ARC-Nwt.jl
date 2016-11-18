export ARC_Nwt
function ARC_Nwt(hh :: AbstractLineFunction,
                 t₀ :: Float64,
                 tₘ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true)

    nf=0
    ng=0
    nh=0

    # Trust region parameters
    eps1 = 0.1
    eps2 = 0.75
    red = 0.1
    aug = 5.0
    Δ = 1.0
    t = t₀

    iter = 0;
    succ=0
    unsucc=0
    hₖ = obj(hh,t)
    gₖ = grad(hh,t)
    kₖ = hess(hh,t)


    q(d) = hₖ + gₖ*d + 0.5*kₖ*d^2

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)
        #step computation
        discr=kₖ^2-4*(gₖ/Δ)
        if discr<0
          discr=kₖ^2+4*(gₖ/Δ)
        end
        if kₖ<0
          dNp=(-kₖ+sqrt(discr))/(2/Δ) #direction de Newton
        else
          dNp=-2*gₖ/(kₖ+sqrt(discr))
        end
        dNn=(-kₖ-sqrt(discr))/(2/Δ) #direction de Newton

        if q(dNp)<q(dNn)
          d=dNp
        else
          d=dNn
        end

        #Numerical reduction computation
        htestTR=obj(hh,t+d)
        gtestTR=grad(hh,t+d)
        pred= gₖ*d + 0.5*kₖ*d^2
        if pred>-1e-10
          ared = (gₖ + gtestTR)*d/2
        else
          ared = htestTR-hₖ
        end

        ratio=ared/pred
        if ratio<eps1
          Δ=red*Δ
          unsucc=unsucc+1
        else
          t=t+d

          gₖ=gtestTR
          hₖ=htestTR
          kₖ=hess(hh,t)

          if ratio>eps2
            Δ=aug*Δ
          end
          succ=succ+1
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
