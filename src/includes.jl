# include generic methods
include("pred_ared.jl")
include("step_computation.jl")

include("pred_ared_stop.jl")
include("step_computation_stop.jl")

# # include ARC methods
# include("ARC1D/ARC-Cub.jl")
# include("ARC1D/ARC-Nwt.jl")
include("ARC1D/ARC-Nwt-Stop.jl")
include("ARC1D/ARC-Sec-Stop.jl")
include("ARC1D/ARC-SecA-Stop.jl")
# include("ARC1D/ARC-Sec.jl")
# include("ARC1D/ARC-SecA.jl")
#
# include("ARC1D/ARC_generic.jl")
# include("ARC1D/ARC_direction_computation.jl")
include("ARC1D/ARC_generic_Stop.jl")
include("ARC1D/ARC_direction_computation_stop.jl")

# include Trust Region methods
# include("TR1D/TR-Cub.jl")
include("TR1D/TR-Nwt.jl")
include("TR1D/TR-Nwt-Stop.jl")
include("TR1D/TR-Sec-Stop.jl")
include("TR1D/TR-SecA-Stop.jl")
# include("TR1D/TR-Sec.jl")
# include("TR1D/TR-SecA.jl")

include("TR1D/TR_direction_computation.jl")
include("TR1D/TR_generic.jl")
include("TR1D/TR_generic_Stop.jl")

include("TR1D/TR_direction_computation_Stop.jl")
#
# # include zoom
# include("Zoom1D/Trouve_IntervalleA.jl")
# # include("Zoom1D/zoom.jl")
# include("Zoom1D/zoom_generic.jl")
# include("Zoom1D/zoom-Cub.jl")
# include("Zoom1D/zoom-Nwt.jl")
# include("Zoom1D/zoom-Sec.jl")
# include("Zoom1D/zoom-SecA.jl")
#
# # include bissection
# include("Bissection/Trouve_Intervalle.jl")
# include("Bissection/Bissect.jl")
# include("Bissection/BissectCub.jl")
# include("Bissection/BissectNwt.jl")
# include("Bissection/BissectSec.jl")
# include("Bissection/BissectSecA.jl")
