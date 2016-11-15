export zoom_Cub
function zoom_Cub(h :: AbstractLineFunction,
                t₀ :: Float64,
                t₁ :: Float64;
                c₁ :: Float64=0.01,
                tol :: Float64=1e-7,
                ϵ :: Float64=1e-10,
                maxiter :: Int=100,
                verbose :: Bool=false)

        if obj(h,t₀)<obj(h,t₁)
          tl=t₀
          th=t₁
        else
          tl=t₁
          th=t₀
        end

        if t₀<t₁
          tmin=t₀
          tmax=t₁
        else
          tmin=t₁
          tmax=t₀
        end

        γ=0.8

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        hl=obj(h,tl)
        dhl=grad(h,tl)
        hh=obj(h,th)
        dhh=grad(h,th)

        t=t₁
        ti=t
        tim1=t₀
        tqnp=t₀
        him1=0
        dhim1=0
        hi=0
        dhi=0
        i=0

        hi=obj(h,t)
        dhi=grad(h,t)
        him1=obj(h,tqnp)
        dhim1=grad(h,tqnp)

        verbose && @printf(" iter      tₗ       tₕ        tᵢ        hₗ        hₕ        hᵢ        dhᵢ\n")
        verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)


        while ((i<50) & (abs(dhi)>tol)) || (i==0)

          if (i>0) & ((hi>h₀+c₁*ti*dh₀) || (hl<=hi))
            th=ti
          else
            if (abs(dhi)<ϵ) & (i>0)
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

          s=ti-tqnp
          y=dhi-dhim1
          α=-s
          z=dhi+dhim1+3*(hi-him1)/α
          discr=z^2-dhi*dhim1
          denom=dhi+dhim1+2*z

          if (discr>0) & (abs(denom)>eps(Float64))
            w=sqrt(discr)
            dC=-s*(dhi+z+sign(α)*w)/(denom)
          else
            dC=-dhi*s/y
          end

          tiB=(tim1+ti)/2
          tiC=ti+dC

          him1=hi
          dhim1=dhi

          if ((tmin<tiC) & (tiC<tmax)) & (obj(h,tiC)<hl)
            ti=tiC
            hi=obj(h,tiC)
            dhi=grad(h,tiC)
            verbose && print_with_color(:green,"C")
          else
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            verbose && print_with_color(:green,"B")
          end

          if t>ti
            if dhi<0
              tim1=t
              tqnp=t
              t=ti
            else
              tqnp=t
              t=ti
            end
          else
            if dhi>0
              tim1=t
              tqnp=t
              t=ti
            else
              tqnp=t
              t=ti
            end
          end

          i=i+1
          verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)
        end
        topt=ti
        iter=i


        return (topt,iter)

end
