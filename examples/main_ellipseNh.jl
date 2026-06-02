using Revise, Dates, FL2D
using FL2D.FLdata

function getfield_default(x, name::Symbol, default)
   return name in propertynames(x) ? getproperty(x, name) : default
end

function run_one_ellipseNh!(io, case; s::Float64, p::Int, nJac::Int=5)

   f!, uex = makeellipseNhfuex(nJac, s)

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
      s=s, p=p,
      δv_near=δv_near,
      δv_close=δv_close,
      δ_near=δ_near,
      δ_intp=δ_intp,
      (f!)=f!, uex=uex,
      dom=case.dom)

   opts = Options(;
      plot=plot_u,
      solver=solver,
      cond_num=cond_num,
      benchmark=benchmark,
      matrixfree=matrixfree,
      s_small=s_small,
      nr=nr, nbd=nbd,
      nr_bdy=nr_bdy,
      ns_near=ns_near,
      pbd=pbd)

   core = solveFL(prob; opts=opts)

   println(io, SolveView(prob, opts, core))
   flush(io)

   return nothing
end

open("test_temp.txt", "w") do io

   s, p = 0.75, 4

   println(io, "s = ", s)
   println(io, "p = ", p)
   flush(io)

   case = (dom=FL2D.ellipseNh(b=[1, 1, 1, 1, 1, 1, 1, 1]),
      N=12, nₚᵣ=128,
      δv_near=0.15, δv_close=8e-3,
      δ_near=0.15, δ_intp=5e-3,
      nr=32,
      s_small=false,
      benchmark=false,
      cond_num=false,
      plot=true)

   run_one_ellipseNh!(io, case; s=s, p=p, nJac=5)

   flush(io)
end

