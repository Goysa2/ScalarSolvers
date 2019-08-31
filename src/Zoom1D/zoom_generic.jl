export zoom_generic

function zoom_generic(h :: AbstractNLPModel,
					  nlpstop :: AbstractStopping,
                  	  t₀ = 2.7, t₁ = 7.5;
					  direction = :Nwt, tol = 1e-7, verbose = false,
					  kwargs...)

    if true in (obj(h, t₀) .< obj(h, t₁))
    	tl = t₀; th = t₁
    else
    	tl = t₁; th = t₀
    end

    hl = obj(h, tl); dhl = grad(h, tl)
    hh = obj(h, th); dhh = grad(h, th)
    h₀ = obj(h, t₀); dh₀ = grad(h, t₀)

    γ = 0.8; t = t₁; tp = t₀; tqnp = t₀; hk = 0; hkm1 = 0; gkm1 = 0; hp = 0
    i = 0

    hk   = obj(h, t);    gk   = grad(h, t)
    hkm1 = obj(h, tqnp); gkm1 = grad(h, tqnp)

	OK = stop!(nlpstop)

    verbose &&
        @printf(" iter        tₗ        tₕ         t        hₗ        hₕ")
    verbose &&
        @printf("         hk         gk\n")
    verbose &&
        @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e\n",
                i, tl[1], th[1], t[1], hl[1], hh[1], hk[1], gk[1])


	while !OK

		if true in (hk .> h₀) || true in (hl .<= hk)
			th = t; tlast = th; hlast = hh
			dhlast = dhh; hh = hk; dhh = gk
		else
			if true in (abs.(gk) .< tol)
				topt = t; iter = i
				return (topt, iter)
			elseif  true in (gk .* (th - tl) .>= 0)
				th = tl; hh = hl; dhh = dhl
			end
		tl = t; tlast = tl; hlast = hl
		dhlast = dhl; hl = hk; dhl = gk
		end

		if direction == :Nwt
			kₖ = hess(h, t)
			dN = vec(- gk ./ kₖ) #direction de Newton
		elseif direction == :Sec
			s = t .- tqnp; y = gk .- gkm1
			dN =-gk .* s ./ y
		elseif direction == :SecA
			s = t .- tqnp; y = gk .- gkm1
			Γ = 3.0 .* (gk .+ gkm1) .* s .- 6.0 .* (hk .- hkm1)
			if y .* s .+ Γ < eps(Float64) .* (s .^ 2)
				yt = y
			else
				yt = y .+ Γ ./ s
			end
			dN = -gk .* s ./ yt
		elseif direction == :Cub
			s = t .- tqnp; y = gk .- gkm1; α = -s
			z = gk .+ gkm1 .+ 3.0 .* (hk .- hkm1) ./ α
			discr = z .^ 2 .- gk .* gkm1; denom = gk .+ gkm1 .+ 2.0 .* z
			if true in (discr .> 0.0) && true in (abs.(denom) .> eps(Float64))
				w = sqrt.(discr)
				dN = -s .* (gk .+ z .+ sign.(α) .* w) ./ (denom)
			else
				dN=-gk .* s ./ y
			end
		end

		if true in ((tp .- t) .* dN .> 0.0) && true in  (dN ./ (tp .- t) .< γ)
      		tplus = t .+ dN
      		hplus = obj(h, tplus); gplus = grad(h, tplus)
      		verbose && println("N")
    	else
      		tplus = (t .+ tp) ./ 2.0
      		hplus = obj(h, tplus); gplus = grad(h, tplus)
      		verbose && println("B")
    	end

		if true in (t .> tp)
      		if true in (gplus .< 0.0)
        		tp = t; tqnp = t; t = tplus
      		else
        		tqnp = t; t = tplus
      		end
    	else
      		if true in (gplus .> 0.0)
        		tp = t; tqnp = t; t = tplus
      		else
        		tqnp = t; t = tplus
      		end
    	end

    	#mise à jour des valeurs
    	hkm1 = hk; gkm1 = gk; hk = hplus; gk = gplus
		OK = update_and_stop!(nlpstop, x = t, fx = hk, gx = gk)
    	i = i + 1
    	verbose &&
        	@printf(" %7.2e %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e %7.2e\n",
                	i, tl[1], th[1], t[1], hl[1], hh[1], hk[1], gk[1])

    end

    return nlpstop.meta.optimal, nlpstop
end #function
