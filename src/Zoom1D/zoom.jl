function zoom(h :: C2LineFunction,
                t₀ :: Float64,
                t₁ :: Float64;
                c₁ :: Float64=0.01,
                tol :: Float64=1e-7,
                ϵ :: Float64=1e-10,
                maxiter :: Int=50,
                verbose :: Bool=false)

        if obj(h,t₀)<obj(h,t₁)
          tl=t₀
          th=t₁
        else
          tl=t₁
          th=t₀
        end

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        hl=obj(h,tl)
        dhl=grad(h,tl)
        hh=obj(h,th)
        dhh=grad(h,th)

        i=0

        ti= (tl+th)/2
        hi=obj(h,ti)
        dhi=grad(h,ti)

        #println("on a les paramètre de zoom")



        verbose && @printf(" iter        tₗ        tₕ         tᵢ        hₗ        hₕ         hᵢ         dhᵢ\n")
        verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)


        while ((i<maxiter) & (abs(dhi)>tol)) || i==0
          #println("on est dans le while de zoom")
          if (hi>h₀+c₁*ti*dh₀) || (hl<=hi)
            # if (hi>h₀+c₁*ti*dh₀)
            #   println("(hi>h₀+c₁*ti*dh₀)")
            # else
            #   println("(hl<=hi)")
            # end
            th=ti
          else
            #println("else")
            if (abs(dhi)<ϵ)
              #print(abs(dhi),"<",-0.99*dh₀)
              #println("1er if dans le else zoom")
              topt=ti
              iter=i
              return (topt,iter)
            elseif dhi*(th-tl)>=0
              #println("elseif dans le else de zoom")
              th=tl
            end
            tl=ti
          end

          hl=obj(h,tl)
          dhl=grad(h,tl)
          hh=obj(h,th)
          dhh=grad(h,th)
          ti= (tl+th)/2
          hi=obj(h,ti)
          dhi=grad(h,ti)



          i=i+1
          #println("on continue d'itérer")
          verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)

        end
        topt=ti
        iter=i

        return (topt,iter)

end
