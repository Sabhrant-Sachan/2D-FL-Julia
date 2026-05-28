using Revise, Dates, FL2D
using FL2D.FLdata

function getfield_default(x, name::Symbol, default)
   return name in propertynames(x) ? getproperty(x, name) : default
end

function run_one_ellipse!(io, case; s::Float64, p::Int, nJac::Int=2)

   f!, uex, _ = makeellipsefuex(nJac, s)

   # New domprop controls
   δv_near = getfield_default(case, :δv_near, 0.15)
   δv_close = getfield_default(case, :δv_close, 8e-3)
   δ_near = getfield_default(case, :δ_near, 0.15)
   δ_intp = getfield_default(case, :δ_intp, 5e-3)

   plot_u = getfield_default(case, :plot, false)
   solver = getfield_default(case, :solver, :direct)
   matrixfree = getfield_default(case, :matrixfree, false)
   s_small = getfield_default(case, :s_small, false)
   benchmark = getfield_default(case, :benchmark, false)
   cond_num = getfield_default(case, :cond_num, false)

   # Quadrature controls
   nr = getfield_default(case, :nr, 32)
   nbd = getfield_default(case, :nbd, 128)
   nr_bdy = getfield_default(case, :nr_bdy, 64)
   ns_near = getfield_default(case, :ns_near, 128)
   pbd = getfield_default(case, :pbd, 2)

   prob = Problem(;
      N=case.N,
      nₚᵣ=case.nₚᵣ,
      s=s,
      p=p,
      δv_near=δv_near,
      δv_close=δv_close,
      δ_near=δ_near,
      δ_intp=δ_intp,
      (f!)=f!,
      uex=uex,
      dom=case.dom
   )

   opts = Options(;
      plot=plot_u,
      solver=solver,
      cond_num=cond_num,
      benchmark=benchmark,
      matrixfree=matrixfree,
      s_small=s_small,
      nr=nr,
      nbd=nbd,
      nr_bdy=nr_bdy,
      ns_near=ns_near,
      pbd=pbd
   )

   core = solveFL(prob; opts=opts)

   println(io, SolveView(prob, opts, core))
   flush(io)

   return nothing
end

# Shared ellipse domains
d1 = FL2D.ellipse(
   b=[1, 1, 1, 1, 1],
   R1=1.0, R2=2.0, L1=2.0, L2=0.8);

d2 = FL2D.ellipse(
   b=[2, 2, 2, 2, 2],
   a=[1, 1, 1, 1, 1],
   R1=1.0, R2=2.0, L1=2.0, L2=0.8);

d3 = FL2D.ellipse(
   b=[4, 3, 4, 3, 3],
   a=[2, 2, 2, 2, 2],
   R1=1.0, R2=2.0, L1=2.0, L2=0.8);

d4 = FL2D.ellipse(
   b=[5, 4, 5, 4, 5],
   a=[3, 3, 3, 3, 4],
   R1=1.0, R2=2.0, L1=2.0, L2=0.8);

d5 = FL2D.ellipse(
   b=[6, 6, 6, 6, 6],
   a=[3, 3, 3, 3, 4],
   R1=1.0, R2=2.0, L1=2.0, L2=0.8);

# Ellipse runs
open("ellipse_outputs0750_direct.txt", "w") do io

   s, p = 0.75, 4

   println(io, "Run started: ", Dates.now())
   println(io, "s = ", s)
   println(io, "p = ", p)
   flush(io)

   common_opts = (
      solver=:direct,
      matrixfree=false,
      benchmark=false,
      cond_num=false,
      plot=false,

      # New domprop defaults
      δv_near=0.15,
      δv_close=8e-3,
      δ_near=0.15,
      δ_intp=5e-3,

      # Quadrature defaults
      nr=32,
      nbd=128,
      nr_bdy=64,
      ns_near=128,
      pbd=2)

   cases = [
      # d1
      (dom=d1, N=3, nₚᵣ=32),
      (dom=d1, N=4, nₚᵣ=32),
      (dom=d1, N=5, nₚᵣ=32),
      (dom=d1, N=6, nₚᵣ=32),
      (dom=d1, N=7, nₚᵣ=64),
      (dom=d1, N=8, nₚᵣ=64),
      (dom=d1, N=9, nₚᵣ=64),
      (dom=d1, N=10, nₚᵣ=128),
      (dom=d1, N=12, nₚᵣ=128),

      # d2
      (dom=d2, N=10, nₚᵣ=128),
      (dom=d2, N=12, nₚᵣ=128),

      # d3
      (dom=d3, N=10, nₚᵣ=128),
      (dom=d3, N=12, nₚᵣ=128),

      # d4
      (dom=d4, N=10, nₚᵣ=128),
      (dom=d4, N=12, nₚᵣ=256),

      # d5
      #(dom=d5, N=10, nₚᵣ=128),
      #(dom=d5, N=12, nₚᵣ=256),
   ]
   cases = [(; c..., common_opts...) for c in cases]

   for case in cases
      run_one_ellipse!(io, case; s=s, p=p, nJac=2)
   end

   println(io, "\n==End of Ellipse==")
   println(io, "Run finished: ", Dates.now())
   flush(io)
end

# \begin{table}[H]
# \centering
# \renewcommand{\arraystretch}{1.25}
# \begin{tabular}{|c|c|c|c|c|}
# \hline
# Target pts & Max rel error & Root MSE error & noc & time (median) \\ \hline
# $57$   & $6.21 \times 10^{-2}$  & $6.25 \times 10^{-3}$  & --   & $0.1946$ \\
# $96$   & $3.40 \times 10^{-2}$  & $3.94 \times 10^{-3}$  & $2.3$  & $0.4286$ \\
# $204$  & $3.43 \times 10^{-3}$  & $2.52 \times 10^{-4}$  & $6.1$  & $0.9152$ \\
# $273$  & $8.30 \times 10^{-4}$  & $4.32 \times 10^{-5}$  & $9.7$  & $2.216$ \\
# $1080$ & $1.43 \times 10^{-8}$  & $4.96 \times 10^{-10}$ & $16.0$ & $29.41$ \\
# $3540$ & $2.73 \times 10^{-10}$ & $5.11 \times 10^{-12}$ & $6.7$  & $142.0$ \\
# $5064$ & $1.33 \times 10^{-10}$ & $2.19 \times 10^{-12}$ & $4.0$  & $222.1$ \\
# $7580$ & $7.99 \times 10^{-12}$ & $1.73 \times 10^{-13}$ & $13.9$ & $393.0$ \\
# $9840$ & $5.42 \times 10^{-12}$ & $6.12 \times 10^{-14}$ & $3.0$  & $606.9$ \\
# \hline
# \end{tabular}
# \end{table}