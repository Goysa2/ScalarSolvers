function trouve_intervalleANwt(h :: C2LineFunction,
                t₀ :: Float64,
                #inc0 :: Float64;
                tₘ :: Float64;
                ϵ :: Float64=1e-10,
                verbose :: Bool=true)

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

        #println("tiN=",tiN," tiB=",tiB)

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)
        ddh₀=hess(h,t₀)



        if (tiN<tim1) || (tiN>tₘ)
          #println("on choisit tiB")
          ti=tiB
          hi=obj(h,tiB)
          dhi=grad(h,tiB)
          ddhi=hess(h,tiB)
          print_with_color(:green,"B")
        else
          #println("on choisit tiN")
          ti=tiN
          hi=obj(h,tiN)
          dhi=grad(h,tiN)
          ddhi=hess(h,tiN)
          #print_with_color(:green,"N")
        end

        #println("on a les paramètres de trouve intervalle")

        verbose && @printf("iter tim1      dhim1      ddhim1      him1       ti      dhi      hi\n")
        verbose && @printf("%4d %7.2e %7.2e  %7.2e %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,ddhim1,him1,ti,dhi,hi)

        if (abs(dhi)<ϵ) & (ddhi>=0)
          topt=ti
          iter=i
          return (topt,iter)
        end

        while (ti<tₘ) || (abs(dhi)>ϵ)
          # if i==0
          #  println("on est dans le while de trouve intervalle")
          # end

          if ((hi>h₀+0.01*ti*dh₀) ||(hi>him1)) & (i>1)
            #println("1er if de trouve_intervalleANwt")
            (topt,iter)=zoom_Nwt(h,tim1,ti)
            # if obj(h,topt)<obj(h,t₀)
            #  println("h(0)<h(t*)")
            # end
            return (topt,iter)
          end

          if (abs(dhi)<=ϵ)
            #println("2e if de trouve_intervalleANwt ")
            iter=i
            topt=ti
            # if obj(h,topt)<obj(h,t₀)
            #   println("h(0)<h(t*)")
            # end
            return (topt,iter)
          end

          #if (dhi*dhim1<=0) #pas la vrai condition...
          if (dhi)>=0
            #println("3e if de trouve_intervalleANwt")
            (topt,iter)=zoom_Nwt(h,ti,tim1)
            # if obj(h,topt)<obj(h,t₀)
            #   println("h(t*)<h(t₀)")
            # end
            return (topt,iter)
          end

          tim1=ti
          him1=hi
          dhim1=dhi
          ddhim1=ddhi


          dN=-dhim1/ddhim1

          #println("dN=", dN)

          tiN=tim1+dN
          tiB=(tim1+tₘ)/2
          #println("tiN=",tiN," tiB=",tiB)

          if (tiN<tim1) || (tiN>tₘ)
            #println("on choisi tiB")
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            ddhi=hess(h,tiB)
            print_with_color(:green,"B")
            #println("ti=",ti,"   hi=",hi,"   dhi=",dhi, "   ddhi=",ddhi)
          else
            #println("on choisi tiN")
            ti=tiN
            hi=obj(h,tiN)
            dhi=grad(h,tiN)
            ddhi=hess(h,tiN)
            print_with_color(:green,"N")
            #println("ti=",ti,"   hi=",hi,"   dhi=",dhi, "   ddhi=",ddhi)
          end

          i=i+1

          verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)


        end
return (ti,i)

end
