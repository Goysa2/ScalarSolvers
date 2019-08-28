export scalar_solvers

scalar_solvers = []

#ARC methods
push!(scalar_solvers, ARC_Nwt)
push!(scalar_solvers, ARC_Cub)
push!(scalar_solvers, ARC_Sec)
push!(scalar_solvers, ARC_SecA)

#TR methods
push!(scalar_solvers, TR_Nwt)
push!(scalar_solvers, TR_Cub)
push!(scalar_solvers, TR_Sec)
push!(scalar_solvers, TR_SecA)

#zoom methods
push!(scalar_solvers, zoom_Nwt)
push!(scalar_solvers, zoom_Cub)
push!(scalar_solvers, zoom_Sec)
push!(scalar_solvers, zoom_SecA)

#find interval and bissection methods
push!(scalar_solvers, bissect)
push!(scalar_solvers, bissect_Cub)
push!(scalar_solvers, bissect_nwt)
push!(scalar_solvers, bissect_sec)
push!(scalar_solvers, bissect_secA)
