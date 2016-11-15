export ARC_Cub
function ARC_Cub(h :: AbstractLineFunction,
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
    aug = 20
    Δ = 1.0
    t = t₀

    iter = 0;
    hₖ = obj(h,t)
    gₖ = grad(h, t)


    secₖ = 1.0
    dN=-gₖ #pour la première iter
    A=0.5
    B=0.0
    q(d) = hₖ + gₖ*d + A*d^2 + B*d^3

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) || (iter == 0)
      #println("premier &")

      #Step computation
      Quad(t)=hₖ+gₖ*t+A*t^2+B*t^3+(1/(4*Δ))*t^4
      dQuad(t)=gₖ+2*A*t+3*B*t^2+(1/Δ)*t^3

      #println(Quad(t))
      #println(dQuad(t))

      dR=roots(dQuad)
      #println("longeur de dR:")
      #println(length(dR))
      vmin=Inf
      #dN=0
      for i=1:length(dR)
        rr=dR[i]
        if isreal(rr)
          rr=real(rr)
          vact=Quad(rr)
          if rr*gₖ<0
            if vact<vmin
              #println("deuxième &")
              dN=rr
              vmin=vact
            end
          end
        end
      end

      d=dN

      #numerical reduction computation
      htestTR=obj(h,t+d)
      gtestTR=grad(h,t+d)
      pred=gₖ*d + A*d^2 + B*d^3
      if pred>-1e-10
        ared=(gₖ+gtestTR)*d/2
      else
        ared=htestTR-hₖ
      end

      ratio=ared/pred
      if ratio<eps1
        Δ=red*Δ
      else
        #two point memory
        tpred=t
        gkm1=gₖ
        hkm1=hₖ

        t=t+d
        gₖ=gtestTR
        hₖ=htestTR

        #Cubic two points interpolation
        s=t-tpred
        y=gₖ-gkm1

        α=-s
        z=gₖ+gkm1+3*(hₖ-hkm1)/α
        discr=z^2-gₖ*gkm1
        denom=gₖ+gkm1+2*z
        B=(1/3)*(gₖ+gkm1+2*z)/(α*α)
        A=-(gₖ+z)/α

        if ratio>eps2
          Δ=aug*Δ
        end
      end

      iter += 1
      verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
