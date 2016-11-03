function zoom_Nwt(h :: C2LineFunction,
                t₀ :: Float64,
                t₁ :: Float64;
                c₁ :: Float64=0.01,
                tol :: Float64=1e-7,
                ϵ :: Float64=1e-10,
                maxiter :: Int=100,
                verbose :: Bool=false)

        #println("on est dans zoom")

        if obj(h,t₀)<obj(h,t₁)
          tl=t₀
          th=t₁
        else
          tl=t₁
          th=t₀
        end

        println("t₀=",t₀," t₁=",t₁)

        if t₀<t₁
          tmin=t₀
          tmax=t₁
        else
          tmin=t₁
          tmax=t₀
        end
        println("tmin=",tmin," tmax=",tmax)

        γ=0.8

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        hl=obj(h,tl)
        dhl=grad(h,tl)
        ddhl=hess(h,tl)
        hh=obj(h,th)
        dhh=grad(h,th)
        ddhh=hess(h,th)

        dN=-dhl/ddhl

        i=0

        tiN= tl + dN
        tiB=(tl+th)/2

        if (tiN>tmin) & (tiN<tmax)
          #println("on est rendu ici")
          ti=tiN
          hi=obj(h,tiN)
          dhi=grad(h,tiN)
          ddhi=hess(h,tiN)
          verbose && print_with_color(:green,"N")
          #println("ti=",ti,"   hi=",hi,"   dhi=",dhi, "   ddhi=",ddhi)
        else
          #println("on est rendu ici 2")
          ti=tiB
          hi=obj(h,tiB)
          dhi=grad(h,tiB)
          ddhi=hess(h,tiB)
          verbose && print_with_color(:green,"B")
          #println("ti=",ti,"   hi=",hi,"   dhi=",dhi, "   ddhi=",ddhi)
        end

        #println("on a les paramètre de zoom")

        verbose && @printf(" iter      tₗ       tₕ        tᵢ        hₗ        hₕ        hᵢ        dhᵢ\n")
        verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)


        while ((i<50) & (abs(dhi)>tol)) || i==0
          # if i==0
          #  println("on est dans le while de zoom")
          # end

          if (hi>h₀+c₁*ti*dh₀) || (hl<=hi)
            # println("if de zoomNwt")
            # if (hi>h₀+c₁*ti*dh₀)
            #  println("(hi>h₀+c₁*ti*dh₀)")
            # else
            #  println("(hl<=hi)")
            # end
            th=ti
            tlast=th
          else
            #println("else zoom_Nwt")
            if (abs(dhi)<ϵ)
               #println(abs(dhi),"<",-0.99*dh₀)
               #println("1er if dans le else zoom")
              topt=ti
              iter=i
              return (topt,iter)
            elseif dhi*(th-tl)>=0
              #println("elseif dans le else de zoom")
              th=tl
            end
            tl=ti
            tlast=ti
          end


          hl=obj(h,tl)
          dhl=grad(h,tl)
          ddhl=hess(h,tl)
          hh=obj(h,th)
          dhh=grad(h,th)
          ddhh=hess(h,th)

          hlast=obj(h,tlast)
          dhlast=grad(h,tlast)
          ddhlast=hess(h,tlast)

          dN=-dhlast/ddhlast

          tiN=ti+dN
          tiB=(tl+th)/2
          #println("tiN=",tiN," tiB=",tiB," tmin=",tmin," tmax=",tmax)


          if ((tmin<tiN) & (tiN<tmax)) & (obj(h,tiN)<hl)
          #if ((ti-tl)*dN>0) & (dN/(ti-tl)<0.8)
            #println("on choisi tiN")
            ti=tiN
            hi=obj(h,tiN)
            dhi=grad(h,tiN)
            ddhi=hess(h,tiN)
            print_with_color(:green,"N")
            #println("ti=",ti,"   hi=",hi,"   dhi=",dhi, "   ddhi=",ddhi)
          else
            #println("on choisi tiB")
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            ddhi=hess(h,tiB)
            print_with_color(:green,"B")
            #println("ti=",ti,"   hi=",hi,"   dhi=",dhi, "   ddhi=",ddhi)
          end



          i=i+1

          verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)
          #println("th-ti=",th-ti)
          # if th==tl || abs(th-tl)<1e-15
          #   print_with_color(:red,"BREAK!")
          #   break
          # end

        end
        topt=ti
        iter=i


        return (topt,iter)

end
