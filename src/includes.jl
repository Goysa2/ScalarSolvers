#include ARC methods
include("ARC1D/ARC-Cub.jl")
include("ARC1D/ARC-Nwt.jl")
include("ARC1D/ARC-Sec.jl")
include("ARC1D/ARC-SecA.jl")

include("ARC1D/new_ARC_Nwt.jl")
include("ARC1D/new_ARC_Sec.jl")
include("ARC1D/new_ARC_SecA.jl")
include("ARC1D/ARC_direction_computation.jl")

#include Trust Region methods
include("TR1D/TR-Cub.jl")
include("TR1D/TR-Nwt.jl")
include("TR1D/TR-Sec.jl")
include("TR1D/TR-SecA.jl")

include("TR1D/new_TR-Nwt.jl")
include("TR1D/new_TR-Sec.jl")
include("TR1D/new_TR-SecA.jl")
include("TR1D/TR_direction_computation.jl")

#include zoom
include("Zoom1D/Trouve_IntervalleA.jl")
include("Zoom1D/Trouve_IntervalleACub.jl")
include("Zoom1D/Trouve_IntervalleANwt.jl")
include("Zoom1D/Trouve_IntervalleASec.jl")
include("Zoom1D/Trouve_IntervalleASecA.jl")
include("Zoom1D/zoom.jl")
include("Zoom1D/zoomCub.jl")
include("Zoom1D/zoomNwt.jl")
include("Zoom1D/zoomSec.jl")
include("Zoom1D/zoomSecA.jl")

#include bissection
include("Bissection/Trouve_Intervalle.jl")
include("Bissection/Bissect.jl")
include("Bissection/BissectCub.jl")
include("Bissection/BissectNwt.jl")
include("Bissection/BissectSec.jl")
include("Bissection/BissectSecA.jl")
