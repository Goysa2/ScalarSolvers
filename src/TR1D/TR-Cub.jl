export TR_Cub
function TR_Cub(h :: AbstractNLPModel,
                nlpstop :: AbstractStopping;
                t₀ = h.meta.x0, verbose = true, eps1 = 0.2, eps2 = 0.8,
                red = 0.5, aug = 2.0, Δ = 1.0)
    Δ  = [Δ]
    t  = t₀; iter = 0;
    fₖ = obj(h, t); gₖ = grad(h, t);

    s = t .- t₀; y = gₖ # just for first iteration 

    dN = -gₖ # for the first iteration
    A = [0.5]; B = [0.0];

    q(d) = fₖ .+ gₖ .* d .+ A .* d .^ 2 .+ B .* d .^ 3

    OK  = update_and_start!(nlpstop, x = t, fx = fₖ, gx = gₖ, g0 = copy(gₖ))

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter, t[1], gₖ[1], Δ[1])

    while !OK

        if true in (q(Δ) .< q(-Δ))
            d = Δ
        else
            d = -Δ
        end

        cub(t) = fₖ .+ gₖ .* t .+ A .* t .^ 2 .+ B .* t .^ 3

        dR = roots(Poly([gₖ[1], 2.0 .* A[1], 3.0 .* B[1]]))

        if isreal(dR)
          dR = real(dR)
          dN2 = dR[1]
          if length(dR) > 1
            if true in ((cub(dR[2]) .< cub(dR[1])))
              dN2 = dR[2]
            end
          end
        else
          dN2 = -gₖ .* s ./ y
        end

        if true in ((abs.(dN2) .< Δ)) && true in ((q(d) .> q(dN2)))
            d = dN2
        end

        if true in (d .== 0.0)
            return OK, nlpstop  # still useful?
        end

        ftestTR = obj(h, t .+ d)
        gtestTR = grad(h, t .+ d)

        pred = gₖ .* d .+ A .* d .^ 2 .+ B .* d .^ 3

        if true in (pred .> - 1e-10)
            ared = (gₖ .+ gtestTR) .* d ./ 2
        else
            ared = ftestTR .- fₖ
        end

        ratio = ared ./ pred
        if true in ((ratio .< eps1)) && true in ((abs.(d) .== Δ))
            Δ = red * Δ
        else
          #two points memory
          tpred = t; gkm1 = gₖ; fkm1 = fₖ
          t = t .+ d; gₖ = gtestTR; fₖ = ftestTR

          #cubic two points interpolation
          s = t .- tpred; y = gₖ .- gkm1

          α = -s
          z = gₖ .+ gkm1 .+ 3. .* (fₖ .- fkm1) ./ α
          discr = z .^ 2 .- gₖ .* gkm1
          denom = gₖ .+ gkm1 .+ 2 .* z
          B = 1 ./ 3 .* (gₖ .+ gkm1 .+ 2. .* z) ./ (α .* α)
          A = -(gₖ .+ z) ./ α

          if true in (ratio .> eps2)
            Δ = aug * Δ
          end
        end

        OK  = update_and_start!(nlpstop, x = t, fx = fₖ, gx = gₖ, g0 = copy(gₖ));

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter, t[1], gₖ[1], Δ[1], pred[1], ared[1])
    end

    optimal = OK
    return optimal, nlpstop
end
