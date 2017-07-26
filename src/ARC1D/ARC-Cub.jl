export ARC_Cub
function ARC_Cub(h :: LineModel,
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
    fₖ = obj(h,t)
    gₖ = grad(h, t)

    secₖ = 1.0
    dN=-gₖ #pour la première iter
    A=0.5
    B=0.0
    q(d) = fₖ + gₖ*d + A*d^2 + B*d^3

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) || (iter == 0)

      #Step computation
      Quad(t)=fₖ+gₖ*t+A*t^2+B*t^3+(1/(4*Δ))*t^4
      dQuad(t)=gₖ+2*A*t+3*B*t^2+(1/Δ)*t^3

      #dR=PolynomialRoots.roots([gₖ,2*A,3*B,(1/Δ)])
      #TECHNIQUEMENT LE CALCUL DES RACINES PAR PolynomialRoots À LA MÊME PRÉCISION QUE SCILAB, MAIS IL DONNE UN SEG FAULT (comme pour TR_Cub)
      p=Poly([gₖ,2*A,3*B,(1/Δ)])
      dR=roots(p)
      # println("dR=",dR)

      # discr=(18*(1/Δ)*(3*B)*(2*A)*gₖ)-4*((3*B)^3)*gₖ+((3*B)^2)*((2*A)^2)-4*(1/Δ)*((2*A)^3)-27*((1/Δ)^2)*((gₖ)^2)
      # discr₀=(3*B)^2-3*(1/Δ)*(2*A)
      # discr₁=2*((3*B)^3)-9*(1/Δ)*(3*B)*(2*A)+27*((1/Δ)^2)*gₖ
      # C₁=cbrt((discr₁-sqrt(complex(discr₁^2-4*(discr₀^3)))/2))
      # C₂=cbrt((discr₁+sqrt(discr₁^2-4*(discr₀^3)))/2)

      #ζ=-0.5+(0.5*sqrt(3))im

      #x₀=-(1/(3*(1/Δ)))*(3*B+((ζ)^0)*C₁+(discr₀)/(((ζ)^0)*C₁))
      #x₀₀=-(1/(3*(1/Δ)))*(3*B+((ζ)^0)*C₂+(discr₀)/(((ζ)^0)*C₂))
      #x₁=-(1/(3*(1/Δ)))*(3*B+((ζ)^1)*C₁+(discr₀)/(((ζ)^1)*C₁))
      #x₂=-(1/(3*(1/Δ)))*(3*B+((ζ)^2)*C₁+(discr₀)/(((ζ)^2)*C₁))

      #dR=[x₀ x₁ x₂]

      vmin=Inf
      dN=0
      for i=1:length(dR)
        rr=dR[i]
        if isreal(rr)
          rr=real(rr)
          vact=Quad(rr)
          if rr*gₖ<0
            if vact<vmin
              dN=rr
              vmin=vact
            end
          end
        end
      end

      d=dN
      #println("d=",d)

      if d==0.0
        return (t,iter+1)
        #L'OUTIL DE CALCUL DES RACINES DU PACKAGE Polynomials ARRONDIE À 0 CERTAINES RACINES D'OÙ LA NOUVELLE CONDITION
      end

      ftestTR=obj(h,t+d)
      gtestTR=grad(h,t+d)

      #println("gtestTR=",gtestTR)

      pred=gₖ*d + A*d^2 + B*d^3

      if pred>-1e-10
        ared=(gₖ+gtestTR)*d/2
      else
        ared=ftestTR-fₖ
      end

      ratio=ared/pred
      if ratio<eps1
        Δ=red*Δ
      else
        #two point memory
        tpred=t
        gkm1=gₖ
        fkm1=fₖ

        t=t+d
        #println("t=",t)
        gₖ=gtestTR
        #println("gₖ=",gₖ)
        fₖ=ftestTR

        #Cubic two points interpolation
        s=t-tpred
        y=gₖ-gkm1

        α=-s
        z=gₖ+gkm1+3*(fₖ-fkm1)/α
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
