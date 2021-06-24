using Gridap
import Gridap: ∇
using Arpack
using Roots
using SparseArrays
using Interpolations
using UnicodePlots

# Add packages
include("dispersion.jl")

include("nonLocal.jl")

include("FEMSolve.jl")

include("refcoeff.jl")

include("solver.jl")

include("coarse2fine.jl")
