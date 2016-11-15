export trouve_intervalleASec
function trouve_intervalleASec(h :: AbstractLineFunction,
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

        verbose && @printf("iter tim1      dhim1      him1       ti      dhi      hi\n")
        verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)

        while i==0 || ((abs(dhi)>ϵ) & (i<50))

          if ((hi>h₀+0.01*t*dh₀) ||(hi>him1)) & (i>1)
            (topt,iter)=zoom_Sec(h,tim1,ti)
            return (topt,iter)
          end

          if (i>0) & (abs(dhi)<=ϵ)
            iter=i
            topt=ti
            return (topt,iter)
          end

          if (dhi)>=0
            (topt,iter)=zoom_Sec(h,ti,tim1)
            return (topt,iter)
          end

          s=t-tqnp
          y=dhi-dhim1

          dS=-dhi*s/y

          tiB=(ti+t)/2
          tiS=t+dS

          him1=hi
          dhim1=dhi

          if (tiS<t₀) || (tiS>tₘ)
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            verbose && print_with_color(:green,"B")
          else
            ti=tiS
            hi=obj(h,tiS)
            dhi=grad(h,tiS)
            verbose && print_with_color(:green,"S")
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
