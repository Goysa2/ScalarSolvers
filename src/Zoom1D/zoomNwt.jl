export zoom_Nwt
function zoom_Nwt(h :: AbstractLineFunction,
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
        ddhl=hess(h,tl)
        hh=obj(h,th)
        dhh=grad(h,th)
        ddhh=hess(h,th)

        dN=-dhl/ddhl

        i=0

        tiN= tl + dN
        tiB=(tl+th)/2

        if (tiN>tmin) & (tiN<tmax)
          ti=tiN
          hi=obj(h,tiN)
          dhi=grad(h,tiN)
          ddhi=hess(h,tiN)
          verbose && print_with_color(:green,"N")
        else
          ti=tiB
          hi=obj(h,tiB)
          dhi=grad(h,tiB)
          ddhi=hess(h,tiB)
          verbose && print_with_color(:green,"B")
        end

        verbose && @printf(" iter      tₗ       tₕ        tᵢ        hₗ        hₕ        hᵢ        dhᵢ\n")
        verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)

        while ((i<50) & (abs(dhi)>tol)) || i==0
          if (hi>h₀+c₁*ti*dh₀) || (hl<=hi)
            th=ti
            tlast=th
          else
            if (abs(dhi)<ϵ)
              topt=ti
              iter=i
              return (topt,iter)
            elseif dhi*(th-tl)>=0
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


          if ((tmin<tiN) & (tiN<tmax)) & (obj(h,tiN)<hl)
            ti=tiN
            hi=obj(h,tiN)
            dhi=grad(h,tiN)
            ddhi=hess(h,tiN)
            verbose && print_with_color(:green,"N")
          else
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            ddhi=hess(h,tiB)
            verbose && print_with_color(:green,"B")
          end

          i=i+1

          verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)

        end
        topt=ti
        iter=i


        return (topt,iter)

end
