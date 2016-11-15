export trouve_intervalleACub
function trouve_intervalleACub(h :: AbstractLineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                ϵ :: Float64=1e-10,
                verbose :: Bool=false)

        ϵ₂=1e-5

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        t=tₘ
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

        verbose && @printf("iter tim1      dhim1      him1       ti      dhi      hi\n")
        verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)

        while i==0 || ((abs(dhi)>ϵ) & (i<50))

          if ((hi>h₀+0.01*t*dh₀) ||(hi>him1)) & (i>1)
            (topt,iter)=zoom_Cub(h,tim1,ti)
            return (topt,iter)
          end

          if (i>0) & (abs(dhi)<=ϵ)
            iter=i
            topt=ti
            return (topt,iter)
          end

          if (dhi)>=0
            (topt,iter)=zoom_Cub(h,ti,tim1)
            return (topt,iter)
          end

          s=t-tqnp
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


          tiB=(ti+t)/2
          tiC=t+dC

          him1=hi
          dhim1=dhi

          if (tiC<t₀) || (tiC>tₘ)
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            verbose && print_with_color(:green,"B")
          else
            ti=tiC
            hi=obj(h,tiC)
            dhi=grad(h,tiC)
            verbose && print_with_color(:green,"C")
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

          verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)


        end
return (ti,i)

end
