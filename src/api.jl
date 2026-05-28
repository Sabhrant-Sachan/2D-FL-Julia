# =========================
# api.jl  (public black-box)
# =========================

export Problem, Options, Result, CoreResult, SolveView,
   solveFL, solveFL_core, solveFL_post, plot_result,
   save_uapp_mat, save_uapp_bin, load_vec_bin

using LinearAlgebra
using Printf
using MAT

#-------------------------
#User-facing configuration
#-------------------------

Base.@kwdef struct Options
   plot::Bool = false
   benchmark::Bool = false
   solver::Symbol = :gmres
   cond_num::Bool = false

   # solver controls, gmres
   reltol::Float64 = 3e-15
   abstol::Float64 = 0.0
   restart::Int = 450
   matrixfree::Bool = false
   s_small::Bool = false

   # quadrature controls
   nr::Int = 32          # regular volume integration: nr × nr
   nbd::Int = 128        # boundary-touching volume integration: nbd × nr
   nr_bdy::Int = 64      # DLP mode-0 regular boundary quadrature
   ns_near::Int = 128    # DLP mode-1 near-boundary split quadrature
   pbd::Int = 2          # DLP boundary smoothing order
end

"""
Problem(; N, nₚᵣ, s, p, δv_near, δv_close, δ_near, δ_intp, f!, dom, uex=nothing)
- `δv_near`: volume near/projection cutoff.
- `δv_close`: volume close/inverse cutoff.
- `δ_near`: boundary DLP near cutoff.
- `δ_intp`: boundary DLP interpolation cutoff.
"""
Base.@kwdef struct Problem
   N::Int
   nₚᵣ::Int
   s::Float64
   p::Int
   δv_near::Float64
   δv_close::Float64
   δ_near::Float64
   δ_intp::Float64
   f!::Function
   dom::abstractdomain
   #This uex is not multiplied by d^s vector
   uex::Union{Nothing,Function,String,Vector{Float64}} = nothing
end

Base.@kwdef struct ErrorMetrics
   max::Float64 = NaN
   relmax::Float64 = NaN
   l2::Float64 = NaN
   rmse::Float64 = NaN
end

Base.@kwdef struct SolveInfo
   solver::Symbol = :gmres
   iters::Int = 0
   converged::Bool = false
   reltol::Float64 = 0
end

Base.@kwdef struct Result
   dp
   d
   problem::Problem
   options::Options
   Uapp::Vector{Float64}                   # Before d^s multiplication
   cond_num::Float64 = -1
   # evaluated fields (for plotting/errors/saving)
   uappv::Vector{Float64}
   uexv::Union{Vector{Float64},Nothing} = nothing
   err::Union{Vector{Float64},Nothing} = nothing
   errors::ErrorMetrics = ErrorMetrics()

   info::SolveInfo = SolveInfo()
end

#CoreResult holds only the essential computation outputs
#that is, stuff to be benchmarked. Plots/error NOT benchmarked!
Base.@kwdef struct CoreResult
   dp
   d
   IntS
   A::Union{Matrix{Float64},Nothing}
   b::Vector{Float64}
   Uapp::Vector{Float64}
   info::SolveInfo
   bench = nothing
end

"""
Save a Float64 vector to a raw .bin file:
Format: n Float64 values when n is length of the vector
"""
function save_uapp_bin(path::String, v::Vector{Float64})
   open(path, "w") do io
      write(io, v)
   end
   return nothing
end

"Load the raw .bin format written by save_uapp_bin."
function load_vec_bin(path::String)::Vector{Float64}
   open(path, "r") do io
      n = filesize(io) ÷ sizeof(Float64)
      v = Vector{Float64}(undef, n)
      read!(io, v)
      return v
   end
end
# -------------------------
# Internal computation helpers
# -------------------------

"Compute uappv from vector Uapp."
function _compute_uappv(dp, d, Uapp::Vector{Float64}, N::Int, s::Float64)
   Np = N^2
   M = d.Npat
   uappv = Vector{Float64}(undef, M * Np)

   @inbounds for i in 1:(M*Np)
      ℓ = cld(i, Np)
      j = i - (ℓ - 1) * Np
      _, r = divrem(j - 1, N)
      j1 = r + 1
      uappv[i] = Uapp[i] * dfunc(d, ℓ, 2 * sinpi((2 * j1 - 1) / (4N))^2, s)
   end

   return uappv
end

"Compute exact vector from a pointwise function uex(x,y) (if given)."
function _compute_uexv_from_function(dp, uex::Function, n::Int)
   uexv = Vector{Float64}(undef, n)
   @inbounds for i in 1:n
      x = dp.tgtpts[1, i]
      y = dp.tgtpts[2, i]
      uexv[i] = uex(x, y)
   end
   return uexv
end

"Compute error vector and metrics."
function _compute_errors(uappv::Vector{Float64}, uexv::Vector{Float64})
   err = Float64.(uappv) .- Float64.(uexv)
   maxerr = maximum(abs.(err))
   relmax = maxerr / maximum(abs.(uexv))
   l2err = norm(err, 2)
   rmse = sqrt(mean(err .^ 2))
   return err, ErrorMetrics(max=maxerr, relmax=relmax, l2=l2err, rmse=rmse)
end

"Load exact vector uexv from .bin path."
function _load_uexv(path::String)
   if endswith(lowercase(path), ".bin")
      return load_vec_bin(path)
   else
      error("Unknown uex file extension: $path. Use .bin.")
   end
end

function _assemble_matrix(dp, d, IntS, s, IV; s_small::Bool=false)
   (; N, Np, M, Mbd, Ni) = IV.IV1

   Lp = Ni + Mbd * N

   A = zeros(Float64, Lp, Lp)

   if s_small
      for k in 1:M
         @views v = A[:, (1+Np*(k-1)):(Np*k)]
         if k in d.kd
            Axbdpth_small!(v, k, IntS, d, dp, s, IV)
         else
            Axintpth_small!(v, k, IntS, d, dp, s, IV)
         end
      end

   else
      for k in 1:M
         @views v = A[:, (1+Np*(k-1)):(Np*k)]
         if k in d.kd
            Axbdpth!(v, k, IntS, d, dp, s, IV)
         else
            Axintpth!(v, k, IntS, d, dp, s, IV)
         end
      end
   end

   for k in 1:Mbd
      @views v = A[:, Ni+N*(k-1)+1:Ni+N*k]
      Axbdop!(v, k, d, dp, s, IV)
   end

   return A
end

# -------------------------
# Phase 1: CORE solve (benchmark this part)
# -------------------------
function solveFL_core(prob::Problem; opts::Options=Options())
   d = prob.dom
   n = prob.nₚᵣ
   dp = domprop(prob.N,
      prob.δv_near, prob.δv_close,
      prob.δ_near, prob.δ_intp, d)

   #b vec computed first and then precomps
   b = bvec(d, dp, prob.s, prob.f!)

   #In this case, a direct solver is always due to high
   #condition number of the discretized matrix
   if opts.s_small

      IntS = precompsLs(d, dp, prob.s, prob.p; n=n)

      IV = compress_vars(d, dp, prob.s, prob.p;
         matrix_form=true, nr=opts.nr, nbd=opts.nbd,
         nr_bdy=opts.nr_bdy, ns_near=opts.ns_near, pbd=opts.pbd)

      A = _assemble_matrix(dp, d, IntS, prob.s, IV; s_small=true)

   else
      IntS = prob.s >= 0.5 ?
             precompsH(d, dp, prob.s, prob.p; n=n) :
             precompsL(d, dp, prob.s, prob.p; n=n)

      if opts.matrixfree
         #Matrix free approach is not for domains with holes in it
         IV = compress_vars(d, dp, prob.s, prob.p;
            matrix_form=false, nr=opts.nr, nbd=opts.nbd,
            nr_bdy=opts.nr_bdy, ns_near=opts.ns_near, pbd=opts.pbd)

         Uapp = copy(b)
         (; N, Np, M, Mbd, Ni) = IV.IV1
         Nt = Ni + Mbd * N
         restart_eff = min(opts.restart, Nt)

         #===
         The argument is an in-place function
             (v, x) -> Ax!(v, x, IntS, d, dp, prob.s, IV)
         which means:
         given an input vector x, compute A*x and write the result into v.
         The two Nt arguments specify the size of the linear operator:
             size(Aop) == (Nt, Nt)
         This is necessary because gmres! needs to know the dimensions of the
         linear system, even though A is not stored explicitly.
         The keyword ismutating=true tells LinearMap that the supplied function
         is an in-place/mutating matvec of the form f!(v, x), rather than a
         function of the form f(x) that returns a new vector.
         ===#
         Aop = LinearMap{Float64}((v, x) -> Ax!(v, x, IntS, d, dp, prob.s, IV),
            Nt, Nt; ismutating=true)

         Uapp, ch = gmres!(Uapp, Aop, b;
            reltol=opts.reltol, abstol=opts.abstol, restart=restart_eff, log=true)

         iters = hasproperty(ch, :iters) ? ch.iters : 0
         conv = hasproperty(ch, :isconverged) ? ch.isconverged : false
         info = SolveInfo(solver=:gmres, iters=iters, converged=conv, reltol=opts.reltol)

         return CoreResult(dp=dp, d=d, IntS=IntS, A=nothing, b=b, Uapp=Uapp, info=info)

      end

      IV = compress_vars(d, dp, prob.s, prob.p;
         matrix_form=true, nr=opts.nr, nbd=opts.nbd,
         nr_bdy=opts.nr_bdy, ns_near=opts.ns_near, pbd=opts.pbd)

      A = _assemble_matrix(dp, d, IntS, prob.s, IV; s_small=false)
   end

   if size(IntS, 2) > 100_000
      IntS = nothing
      GC.gc()
   end

   # solve
   info = SolveInfo(solver=opts.solver)

   if opts.solver == :direct
      if d.nh == 0
         Uapp = A \ b
         info = SolveInfo(solver=:direct, iters=0, converged=true, reltol=0)

      else
         # There are holes in the domain.
         # A is rank-deficient by d.nh. We use SLP basis functions
         # to modify the RHS so that it is in the image of A.

         (; N, Np, M, Mbd, Ni) = IV.IV1

         # bZ : The single layer potential on all the
         #      target points (including the boundary).
         #      It's a rectangular matrix of size
         #      (Npat*N*N + Mbd*N) * nh where nh is the
         #      number of holes and Npat*N*N + Mbd*N is total
         #      number of target points.
         bZ = SLPeval(d, dp, IV)

         if !opts.s_small && prob.s >= 0.5
            bZ[(Ni+1):(Ni+Mbd*N), 1:d.nh] .= 0.0
         end

         Lp = Ni + Mbd * N
         indx = (Lp-d.nh+1):Lp

         # Use pivoted QR for rank-deficient A.
         F = qr(A, ColumnNorm())

         Qa = Matrix(F.Q)
         Ra = F.R
         Pa = Matrix(F.P)

         if size(A, 2) > 30_000
            A = nothing
            GC.gc()
         end


         B_SLP = Qa' * b
         Beta_SLP = Qa' * bZ

         if size(Qa, 2) > 30_000
            Qa = nothing
            GC.gc()
         end

         # Choose SLP coefficients so bottom compatibility equations vanish:
         # Q₂ᵀ(b + bZ*R) = 0.
         RNh = -(Beta_SLP[indx, :] \ B_SLP[indx])

         beq = B_SLP + Beta_SLP * RNh

         # Enforce exact numerical compatibility in the artificial rows.
         beq[indx] .= 0.0

         for j in 1:d.nh
            Ra[indx[j], indx[j]] = 1.0
         end

         Uapp = Pa * (Ra \ beq)

         info = SolveInfo(solver=:direct, iters=0, converged=true, reltol=0.0)
      end

   elseif opts.solver == :gmres

      if d.nh == 0
         Uapp = copy(b)
         (; N, Np, M, Mbd, Ni) = IV.IV1
         Lp = Ni + Mbd * N
         restart_eff = min(opts.restart, Lp)

         Uapp, ch = gmres!(Uapp, A, b;
            reltol=opts.reltol, abstol=opts.abstol, restart=restart_eff, log=true)

         iters = hasproperty(ch, :iters) ? ch.iters : 0
         conv = hasproperty(ch, :isconverged) ? ch.isconverged : false
         info = SolveInfo(solver=:gmres, iters=iters, converged=conv, reltol=opts.reltol)
      else
         error("solver=$(opts.solver) not possible. Use direct solver.")
      end
   else
      error("Unknown solver=$(opts.solver). Use :gmres or :direct.")
   end

   #Too much RAM usage for modern laptops. ~15GB
   if A !== nothing && size(A, 2) > 40_000
      A = nothing
      GC.gc()
   end

   return CoreResult(dp=dp, d=d, IntS=IntS, A=A, b=b, Uapp=Uapp, info=info)

end

# -------------------------
# Phase 2: Postprocessing (NOT benchmarking this part)
# -------------------------

function solveFL_post(prob::Problem, core::CoreResult; opts::Options=Options())
   dp, d, A, Uapp = core.dp, core.d, core.A, core.Uapp

   # Evaluate solution on target points
   uappv = _compute_uappv(dp, d, Uapp, prob.N, prob.s)

   # Exact / errors
   uexv = nothing
   err = nothing
   em = ErrorMetrics()

   if opts.cond_num && A !== nothing
      CN = cond(A)
   else
      CN = -1
   end

   if prob.uex === nothing
      # no errors
   elseif prob.uex isa Function
      uexv = _compute_uexv_from_function(dp, prob.uex, length(uappv))
      err, em = _compute_errors(uappv, uexv)
   elseif prob.uex isa AbstractVector
      uexv = _compute_uappv(dp, d, prob.uex, prob.N, prob.s)
      err, em = _compute_errors(uappv, uexv)
   elseif prob.uex isa String
      Uexv = _load_uexv(prob.uex)
      uexv = _compute_uappv(dp, d, Uexv, prob.N, prob.s)
      err, em = _compute_errors(uappv, uexv)
   else
      error("Unsupported uex type: $(typeof(prob.uex)). Use nothing, Function, or String path.")
   end

   res = Result(d=core.d, dp=core.dp, problem=prob, options=opts, Uapp=Uapp, uappv=uappv,
      uexv=uexv, err=err, errors=em, info=core.info, cond_num=CN)

   return res
end

#This function takes in A finer solution, domain and N
#and then return a solution on the coarser mesh 
function build_ref_coarse(Uf::Vector{Float64}, df::abstractdomain,
   Nf::Int, dc::abstractdomain, Nc::Int)

   Mf, Mc = df.Npat, dc.Npat

   Npc, Npf = Nc * Nc, Nf * Nf

   Lₚᵪ = Mc * Npc

   chebcoeff = zeros(Float64, Mf * Npf)

   for k in 1:Mf
      idx = (k-1)*Npf+1:k*Npf
      fv = reshape(@views(Uf[idx]), Nf, Nf)

      # columns first (dim = 2), DCT-II
      c = (1 / Nf) .* FFTW.r2r(fv, FFTW.REDFT10, 2)
      c[:, 1] ./= 2

      # rows next (dim = 1), DCT-II
      c = (1 / Nf) .* FFTW.r2r(c, FFTW.REDFT10, 1)
      c[1, :] ./= 2

      chebcoeff[idx] .= vec(c)
   end

   Uc = Vector{Float64}(undef, Lₚᵪ)
   Tu = Vector{Float64}(undef, Nf)
   Tv = Vector{Float64}(undef, Nf)

   for col in 1:Lₚᵪ
      # Get patch of point in coarser pat.
      kc = cld(col, Npc)

      jj = col - (kc - 1) * Npc

      q, r = divrem(jj - 1, Nc)

      # Point in parametric space of kc patch
      x1 = cospi((2 * r + 1) / (2Nc))
      x2 = cospi((2 * q + 1) / (2Nc))

      # Convert point in coarser mesh to region pts
      x1r, x2r, r = ptconv(dc, x1, x2, kc, "to_reg")

      # Convert region point to patch pt
      u, v, k = ptconv(df, x1r, x2r, r, "to_pth")

      cf = reshape(chebcoeff[(k-1)*Npf+1:k*Npf], Nf, Nf)
      ChebyTN!(Tu, Nf, u)
      ChebyTN!(Tv, Nf, v)

      Uc[col] = dot(Tu, cf, Tv)

   end

   return Uc
end

# -------------------------
# Public wrapper: black-box call
# Benchmarks ONLY core solve
# -------------------------

function solveFL(prob::Problem; opts::Options=Options())
   bench = nothing
   core = nothing

   if opts.benchmark

      opts_core = Options(
         plot=false,
         benchmark=false,
         solver=opts.solver,
         cond_num=opts.cond_num,
         reltol=opts.reltol,
         abstol=opts.abstol,
         restart=opts.restart,
         matrixfree=opts.matrixfree,
         s_small=opts.s_small,
         nr=opts.nr,
         nbd=opts.nbd,
         nr_bdy=opts.nr_bdy,
         ns_near=opts.ns_near,
         pbd=opts.pbd)

      core = solveFL_core(prob; opts=opts_core)  # compute once for actual result

      # benchmarking only the core solve
      bench = BenchmarkTools.@benchmark solveFL_core($prob; opts=$opts_core) evals = 1 samples = 3 seconds = 10_000

   else
      core = solveFL_core(prob; opts=opts)

   end

   return CoreResult(dp=core.dp, d=core.d, IntS=core.IntS, A=core.A, b=core.b, Uapp=core.Uapp, info=core.info, bench=bench)
end

# -------------------------
# Plot convenience wrapper
# -------------------------

import Base: show

struct SolveView
   prob::Problem
   opts::Options
   core::CoreResult
end

function show(io::IO, ::MIME"text/plain", v::SolveView)
   println(io, "============ FractionalLaplacian2D ============")

   # build a postprocessed Result using current opts
   res = solveFL_post(v.prob, v.core; opts=v.opts)

   println(io, "\n----------- Domain -----------")
   try
      show(io, MIME"text/plain"(), v.core.d)
      println(io)
   catch
      println(io, v.core.d)
   end

   prob = res.problem
   println(io, "\n----------- Parameters -----------")
   @printf(io, "N=%d, nₚᵣ=%d, s=%.6g, p=%d,\n",
      prob.N, prob.nₚᵣ, prob.s, prob.p)
   @printf(io, "δv_near=%.6g, δv_close=%.6g\n",
      prob.δv_near, prob.δv_close)
   if v.opts.s_small
      @printf(io, "--- small s case used --- \n")
   else
      @printf(io, "Precomps nodes: r = nₚᵣ, θ = 2nₚᵣ,\n")
   end
   @printf(io, "Lᵢₙ=%d, δ_near=%.6g, δ_intp=%.6g\n",
      v.core.dp.Lᵢₙ, prob.δ_near, prob.δ_intp)

   println(io, "\n----------- Quadrature -----------")
   @printf(io, "Vol. regular patches  : nr = %d  (%d×%d grid)\n",
      v.opts.nr, v.opts.nr, v.opts.nr)

   @printf(io, "Vol. boundary patches : nbd = %d (%d×%d grid)\n",
      v.opts.nbd, v.opts.nbd, v.opts.nr)

   @printf(io, "Boundary DLP mode-0   : nr_bdy = %d\n", v.opts.nr_bdy)
   @printf(io, "Boundary DLP mode-1   : ns_near = %d\n", v.opts.ns_near)
   @printf(io, "Boundary DLP smoothing: pbd = %d\n", v.opts.pbd)

   println(io, "\n----------- System -----------")
   n = length(v.core.b)
   println(io, "Target points (Vars.) : $n")
   if v.core.IntS === nothing
      println(io, "IntS columns       : freed after solve")
   end

   println(io, "Stored near/close pts : ", size(v.core.dp.invpts, 2))
   println(io, "IntS layout columns   : ", v.core.dp.pthgo[end] - 1)

   if res.cond_num != -1
      if res.cond_num > 1000
         @printf(io, "Cond. num 2-norm: %.2e\n", res.cond_num)
      else
         @printf(io, "Cond. num 2-norm: %f\n", res.cond_num)
      end
   end

   println(io, "\n----------- Solve -----------")
   info = res.info
   println(io, "Solver : ", info.solver)
   if info.solver == :gmres
      println(io, "Iters  : ", info.iters)
      println(io, "Conv   : ", info.converged)
      println(io, "Reltol : ", info.reltol)
   end

   println(io, "\n----------- Error -----------")
   em = res.errors
   if isnan(em.max)
      println(io, "(No exact solution provided)")
   else
      @printf(io, "Max Error   : %.2e\n", em.max)
      @printf(io, "Max Rel err : %.2e\n", em.relmax)
      @printf(io, "L^2 error   : %.2e\n", em.l2)
      @printf(io, "Root MSE err: %.2e\n", em.rmse)
   end

   if v.core.bench !== nothing
      println(io, "\n----------- Benchmark -----------")
      try
         show(io, MIME"text/plain"(), v.core.bench)
         println(io)
      catch
         println(io, v.core.bench)
      end
   end

   if res.options.plot
      #Everytime show is called, the paraview files are changed!
      u_to_paraview(v.core.dp, v.core.d, v.core.Uapp, v.prob.s)
   end

   println(io, "===============================================")
end

show(io::IO, v::SolveView) = show(io, MIME"text/plain"(), v)
