export trouve_intervalleANwt
function trouve_intervalleANwt(h :: AbstractLineFunction,
                t₀ :: Float64,
                #inc0 :: Float64;
                tₘ :: Float64;
                ϵ :: Float64=1e-10,
                verbose :: Bool=false)

        ϵ₂=1e-5

        tim1=t₀
        him1=obj(h,tim1)
        dhim1=grad(h,tim1)
        ddhim1=hess(h,tim1)

        dN=-dhim1/ddhim1

        tiN=tim1+dN
        tiB=(tim1+tₘ)/2
        i=0

        distt0=abs(tₘ-t₀)
        disttiN=abs(tₘ-tiN)

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)
        ddh₀=hess(h,t₀)



        if (tiN<tim1) || (tiN>tₘ)
          ti=tiB
          hi=obj(h,tiB)
          dhi=grad(h,tiB)
          ddhi=hess(h,tiB)
          verbose && print_with_color(:green,"B")
        else
          ti=tiN
          hi=obj(h,tiN)
          dhi=grad(h,tiN)
          ddhi=hess(h,tiN)
          verbose && print_with_color(:green,"N")
        end

        verbose && @printf("iter tim1      dhim1      ddhim1      him1       ti      dhi      hi\n")
        verbose && @printf("%4d %7.2e %7.2e  %7.2e %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,ddhim1,him1,ti,dhi,hi)

        if (abs(dhi)<ϵ) & (ddhi>=0)
          topt=ti
          iter=i
          return (topt,iter)
        end

        while (ti<tₘ) || (abs(dhi)>ϵ)

          if ((hi>h₀+0.01*ti*dh₀) ||(hi>him1)) & (i>1)
            (topt,iter)=zoom_Nwt(h,tim1,ti)
            return (topt,iter)
          end

          if (abs(dhi)<=ϵ)
            iter=i
            topt=ti
            return (topt,iter)
          end

          if (dhi)>=0
            (topt,iter)=zoom_Nwt(h,ti,tim1)
            return (topt,iter)
          end

          tim1=ti
          him1=hi
          dhim1=dhi
          ddhim1=ddhi

          dN=-dhim1/ddhim1

          tiN=tim1+dN
          tiB=(tim1+tₘ)/2

          if (tiN<tim1) || (tiN>tₘ)
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            ddhi=hess(h,tiB)
            verbose && print_with_color(:green,"B")
          else
            ti=tiN
            hi=obj(h,tiN)
            dhi=grad(h,tiN)
            ddhi=hess(h,tiN)
            verbose && print_with_color(:green,"N")
          end

          i=i+1

          verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)

        end
return (ti,i)

end
