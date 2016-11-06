export all_solvers

all_solvers=[]

#ARC methods
push!(all_solvers,ARC_Cub)
push!(all_solvers,ARC_Nwt)
push!(all_solvers,ARC_Sec)
push!(all_solvers,ARC_SecA)

#TR methods
push!(all_solvers,TR_Cub)
push!(all_solvers,TR_Nwt)
push!(all_solvers,TR_Sec)
push!(all_solvers,TR_SecA)

#zoom methods
push!(all_solvers,trouve_intervalleA)
push!(all_solvers,trouve_intervalleACub)
push!(all_solvers,trouve_intervalleANwt)
push!(all_solvers,trouve_intervalleASec)
push!(all_solvers,trouve_intervalleASecA)

#find interval and bissection methods
push!(all_solvers,bissect)
push!(all_solvers,bissect_Cub)
push!(all_solvers,bissect_nwt)
push!(all_solvers,bissect_sec)
push!(all_solvers,bissect_secA)
