module FL2D

# --- Standard libs ---
using LinearAlgebra
using Printf
using Random

# --- Numerics / solvers ---
using FFTW, MAT
using IterativeSolvers
using Subroutines
using ApproxFun

# --- Special functions used internally ---
using SpecialFunctions: gamma
using HypergeometricFunctions: pFq

# --- Plotting (keep if you want plot flags in-package) ---
using GLMakie
using LaTeXStrings
using ColorSchemes
using LinearMaps
# --- Finding time taken by FL solver ---
using BenchmarkTools
using WriteVTK, VTKBase

# --------------------------
# Include INTERNAL code files
# --------------------------
# Data / domains / compression
include("domains/domain.jl")

include("FLdata.jl")
include("compressvars.jl")

# Core assembly pieces
include("domprop.jl")
include("bvec.jl")
include("precomps.jl")

#Boundary Integral files
include("DLPfuncs.jl")
include("SLPfuncs.jl")

# Matrix assembly kernels
include("matrixvecprod/Axintpth.jl")
include("matrixvecprod/Axbdpth.jl")
include("matrixvecprod/Axbdop.jl")
include("matrixvecprod/Ax.jl")

# --------------------------
# Public API (black-box)
# --------------------------
include("api.jl")

end # module
