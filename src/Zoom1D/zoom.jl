export zoom
function zoom(h :: AbstractLineFunction,
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

        verbose && @printf(" iter        tₗ        tₕ         tᵢ        hₗ        hₕ         hᵢ         dhᵢ\n")
        verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)


        while ((i<maxiter) & (abs(dhi)>tol)) || i==0
          if (hi>h₀+c₁*ti*dh₀) || (hl<=hi)
            th=ti
          else
            if (abs(dhi)<ϵ)
              topt=ti
              iter=i
              return (topt,iter)
            elseif dhi*(th-tl)>=0
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
          verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)

        end
        topt=ti
        iter=i

        return (topt,iter)

end
