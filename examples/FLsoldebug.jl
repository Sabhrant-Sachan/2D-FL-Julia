using Revise, FL2D

import FL2D.FLdata as FLdata
import FL2D.GLMakie as GLMakie
# ------------------------------------------------------------
# Plot helper
# ------------------------------------------------------------
function plot_mtx_err(Errs::AbstractMatrix{Float64};
   title::String="Matrix Error")

   row_indices = 1:size(Errs, 1)
   col_indices = 1:size(Errs, 2)

   fig = GLMakie.Figure(size=(1300, 600), fontsize=16)

   ax2d = GLMakie.Axis(fig[1, 1],
      title="2D Error Heatmap: $title",
      xlabel="Matrix Row Index",
      ylabel="Matrix Column Index")

   ax3d = GLMakie.Axis3(fig[1, 3],
      title="3D Error Landscape: $title",
      xlabel="Matrix Row Index",
      ylabel="Matrix Column Index",
      zlabel="log10(Error)",
      aspect=(1, 1, 0.8),
      perspectiveness=0.5)

   cmap = :inferno

   hm = GLMakie.heatmap!(ax2d, row_indices, col_indices, Errs, colormap=cmap)
   surface!(ax3d, row_indices, col_indices, Errs, colormap=cmap)

   GLMakie.Colorbar(fig[1, 2], hm,
      label="log10(Relative Error + floor)",
      width=20,
      tickalign=1)

   GLMakie.display(fig)

   return nothing
end

# ------------------------------------------------------------
# Manufactured-solution selector
# ------------------------------------------------------------
function get_fuex(kind::Symbol, nJac::Int, s::Float64)
   if kind == :disc
      return FLdata.makediscfuex(nJac, s)
   elseif kind == :kite
      return FLdata.makekitefuex(nJac, s)
   elseif kind == :squircle
      return FLdata.makesquirclefuex(nJac, s, 4)
   elseif kind == :ellipse
      return FLdata.makeellipsefuex(nJac, s)
   elseif kind == :peanut
      return FLdata.makepeanutfuex(nJac, s)
   elseif kind == :bean
      return FLdata.makebeanfuex(nJac, s)
   elseif kind == :star
      return FLdata.makestarfuex(nJac, s)
   elseif kind == :annulus
      return FLdata.makeannulusfuex(nJac, s)
   elseif kind == :ellipseNh
      return FLdata.makeellipseNhfuex(nJac, s)
   else
      error("Unknown manufactured-solution kind: $kind")
   end
end

# ------------------------------------------------------------
# Common debug runner
# ------------------------------------------------------------
function debug_domain(d;
   kind::Symbol,
   N::Int=12,
   Lᵢₙ::Int=5,
   δv_near::Float64=0.15,
   δv_close::Float64=8e-3,
   δ_near::Float64=0.15,
   δ_intp::Float64=5e-3,
   draw_domain::Bool=true,
   bdquadtest::Bool=false,
   invtest::Bool=true,
   dlptest::Bool=true,
   bvectest::Bool=false,
   precompstest::Bool=false,
   slptest::Bool=false,
   plot_bvec::Bool=true,
   plot_precomp_err::Bool=true,
   s_bvec::Float64=0.75,
   s_small::Float64=0.01,
   s_pre::Float64=0.1,
   p::Int=4,
   nref_bvec::Int=256,
   nref_pre::Int=256,
   ntest_bvec=(16, 32, 64, 128),
   ntest_pre=(16, 32, 64, 128),
   nJac::Int=2)

   println("================================================")
   println("Debugging domain kind = ", kind)
   println("Domain type           = ", typeof(d))
   println("N                     = ", N)
   println("================================================")

   # ------------------------------------------------------------
   # 1. Geometry / mapping tests
   # ------------------------------------------------------------
   if draw_domain
      FL2D.draw(d, 1)
      FL2D.drawbd(d, 1)
   end

   display(d)

   FL2D.chk_map(d)

   dp = FL2D.domprop(N, δv_near, δv_close, δ_near, δ_intp, d; Lᵢₙ=Lᵢₙ)

   display(dp)

   #FL2D.memory_report(dp)

   if bdquadtest
      for k in 1:d.Npat
         FL2D.plotns(dp, d, k)
      end
   end

   if invtest
      FL2D.chkinvpts(dp, d; flag=true)
   end

   if dlptest
      FL2D.testDLP(d, dp; nr=64)
   end

   # ------------------------------------------------------------
   # 2. bvec convergence
   # ------------------------------------------------------------
   if bvectest
      f!, _ = get_fuex(kind, nJac, s_bvec)

      println()
      println("bvec convergence test")
      println("  s = ", s_bvec)
      println("  nJac = ", nJac)
      println("  reference n = ", nref_bvec)

      bex = FL2D.bvec(d, dp, s_bvec, f!; n=nref_bvec)

      if plot_bvec
         FL2D.plotfunc(dp, d, bex)
      end

      for n in ntest_bvec
         b = FL2D.bvec(d, dp, s_bvec, f!; n=n)
         err = maximum(abs.(b .- bex))
         println("  n = ", n, "  max abs err = ", err)
      end

      # println()
      # println("bvec_small convergence test")
      # println("  s = ", s_small)
      # println("  reference n = ", nref_bvec)

      # fsmall!, _ = get_fuex(kind, nJac, s_small)

      # bex_small = FL2D.bvec_small(d, dp, s_small, fsmall!; n=nref_bvec)

      # if plot_bvec
      #    FL2D.plotfunc(dp, d, bex_small)
      # end

      # for n in ntest_bvec
      #    b = FL2D.bvec_small(d, dp, s_small, fsmall!; n=n)
      #    err = maximum(abs.(b .- bex_small))
      #    println("  n = ", n, "  max abs err = ", err)
      # end
   end

   # ------------------------------------------------------------
   # 3. precomps convergence
   # ------------------------------------------------------------
   if precompstest
      println()
      println("precompsLs convergence test")
      println("  s = ", s_pre)
      println("  p = ", p)
      println("  reference n = ", nref_pre)

      IntSex = FL2D.precompsLs(d, dp, s_pre, p; n=nref_pre)

      for n in ntest_pre
         IntS = FL2D.precompsLs(d, dp, s_pre, p; n=n)
         err = maximum(abs.(IntS .- IntSex))
         println("  n = ", n, "  max abs err = ", err)
      end

      IntS = FL2D.precompsLs(d, dp, s_pre, p; n=maximum(ntest_pre))
      total_err = maximum(abs.(IntS .- IntSex))
      println("  final tested n = ", maximum(ntest_pre), "  max abs err = ", total_err)

      println()
      println("Near/close blocks by patch:")
      for k in 1:d.Npat
         S = dp.pthgo[k] + dp.N^2
         E = dp.pthgo[k+1] - 1

         if S <= E
            err = maximum(abs.(IntS[:, S:E] .- IntSex[:, S:E]))
            println("  patch ", k, "  near/close abs err = ", err)
         end
      end

      println()
      println("Singular/self blocks by patch:")
      for k in 1:d.Npat
         S = dp.pthgo[k]
         E = dp.pthgo[k] + dp.N^2 - 1

         abs_err = maximum(abs.(IntS[:, S:E] .- IntSex[:, S:E]))
         rel_err = abs_err / maximum(abs.(IntSex[:, S:E]))

         println("  patch ", k, "  abs = ", abs_err, "  rel = ", rel_err)
      end

      if plot_precomp_err
         Errs = log10.(abs.(IntS .- IntSex) ./ maximum(abs.(IntSex)) .+ 1e-16)
         plot_mtx_err(Errs; title="precompsLs, s=$s_pre")
      end
   end

   # ------------------------------------------------------------
   # 4. SLP test for hole domains
   # ------------------------------------------------------------
   if slptest
      if d.nh > 0
         println()
         println("SLP test")

         IV = FL2D.compress_vars(d, dp, s_pre, p; matrix_form=true)
         FL2D.SLPtest(d, dp, IV)
      else
         println()
         println("SLP test skipped: domain has no holes.")
      end
   end

   println("Done.")
   return nothing
end

# ------------------------------------------------------------
# Disc
# ------------------------------------------------------------
function debug_disc(domain::Symbol=:coarse; kwargs...)
   d = if domain == :coarse
      FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)
   elseif domain == :fine
      FL2D.disc(b=[6, 6, 6, 6, 5], a=[4, 4, 4, 4, 4], L1=0.8, L2=0.8)
   else
      error("Unknown disc domain = $domain. Use :coarse or :fine.")
   end

   return debug_domain(d; kind=:disc, kwargs...)
end

# ------------------------------------------------------------
# Squircle
# ------------------------------------------------------------
function debug_squircle(domain::Symbol=:coarse; P::Float64=4.0, kwargs...)
   d = if domain == :coarse
      FL2D.squircle(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8, P=P)
   elseif domain == :fine
      FL2D.squircle(b=[6, 6, 6, 6, 5], a=[4, 4, 4, 4, 4], L1=0.8, L2=0.8, P=P)
   else
      error("Unknown squircle domain = $domain. Use :coarse or :fine.")
   end

   return debug_domain(d; kind=:squircle, p=Int(round(P)), kwargs...)
end

# ------------------------------------------------------------
# Ellipse
# ------------------------------------------------------------
function debug_ellipse(domain::Symbol=:coarse; kwargs...)
   d = if domain == :coarse
      FL2D.ellipse(
         b=[2, 2, 2, 2, 2],
         a=[1, 1, 1, 1, 1],
         R1=1.0, R2=2.0, L1=2.0, L2=0.8
      )
   elseif domain == :fine
      FL2D.ellipse(
         b=[4, 3, 4, 3, 4],
         a=[2, 2, 2, 2, 3],
         R1=1.0, R2=2.0, L1=2.0, L2=0.8)
   else
      error("Unknown ellipse domain = $domain. Use :coarse or :fine.")
   end

   return debug_domain(d; kind=:ellipse, kwargs...)
end

# ------------------------------------------------------------
# Annulus
# ------------------------------------------------------------
function debug_annulus(domain::Symbol=:coarse; kwargs...)
   d = if domain == :coarse
      FL2D.annulus(b=[1, 1, 1, 1, 1, 1, 1, 1])
   elseif domain == :fine
      FL2D.annulus(b=[4, 4, 4, 4, 4, 4, 4, 4])
   else
      error("Unknown annulus domain = $domain. Use :coarse or :fine.")
   end

   return debug_domain(d; kind=:annulus, slptest=true, kwargs...)
end

# ------------------------------------------------------------
# EllipseNh
# ------------------------------------------------------------
function debug_ellipseNh(domain::Symbol=:coarse; kwargs...)
   d = if domain == :coarse
      FL2D.ellipseNh(b=[1, 1, 1, 1, 1, 1, 1, 1])
   else
      error("Unknown ellipseNh domain = $domain. Add a :fine constructor if needed.")
   end

   return debug_domain(d; kind=:ellipseNh, slptest=true, kwargs...)
end

# ------------------------------------------------------------
# Peanut
# ------------------------------------------------------------
function debug_peanut(domain::Symbol=:coarse; kwargs...)
   d = if domain == :coarse
      FL2D.peanut(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
   elseif domain == :fine
      dd = FL2D.peanut(b=[6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6])
      FL2D.refine!(dd, 2, 2, [67, 96, 187, 216])
      dd
   else
      error("Unknown peanut domain = $domain. Use :coarse or :fine.")
   end

   return debug_domain(d; kind=:peanut, kwargs...)
end

# ------------------------------------------------------------
# Kite
# ------------------------------------------------------------
function debug_kite(domain::Symbol=:coarse; kwargs...)
   d = if domain == :coarse
      FL2D.kite(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
   elseif domain == :fine
      dd = FL2D.kite(
         b=[4, 5, 9, 9, 5, 4, 5, 5, 5, 5, 4, 4],
         a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 3, 3]
      )
      FL2D.refine!(dd, 1, 2, [9, 12, 13, 16, 181, 184, 185, 188])
      dd
   else
      error("Unknown kite domain = $domain. Use :coarse or :fine.")
   end

   return debug_domain(d; kind=:kite, kwargs...)
end

# ------------------------------------------------------------
# Bean
# ------------------------------------------------------------
function debug_bean(domain::Symbol=:coarse; kwargs...)
   d = if domain == :coarse
      FL2D.bean(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
   elseif domain == :fine
      dd = FL2D.bean(
         b=[4, 5, 5, 5, 5, 4, 5, 9, 9, 5, 4, 4],
         a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 4, 4]
      )
      FL2D.refine!(dd, 1, 2, [9, 12, 13, 16, 125, 128, 129, 132])
      dd
   else
      error("Unknown bean domain = $domain. Use :coarse or :fine.")
   end

   return debug_domain(d; kind=:bean, kwargs...)
end

# ------------------------------------------------------------
# Star
# ------------------------------------------------------------
function debug_star(domain::Symbol=:coarse; kwargs...)
   d = if domain == :coarse
      FL2D.star(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
   elseif domain == :petal4
      FL2D.star(
         b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         P=4, tht1=0.24, tht2=0.225, L1=0.5, L2=0.5
      )
   elseif domain == :petal6
      FL2D.star(
         b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         P=6, tht1=0.14, tht2=0.15, L1=0.3, L2=0.52
      )
   else
      error("Unknown star domain = $domain. Use :coarse, :sharp4, or :sharp6.")
   end

   return debug_domain(d; kind=:star, kwargs...)
end

# d = FL2D.disc(b=[6, 6, 6, 6, 5], a=[4, 4, 4, 4, 4], L1=0.8, L2=0.8)

# dp = FL2D.domprop(12, 8e-3, 0.15, 0.15, 5e-3, d; Lᵢₙ=5)

# f!, _ = get_fuex(:disc, 9, 0.5)

# FL2D.plotfunc(dp, d, f!)

debug_disc(:coarse; N=12, bdquadtest=false)
debug_disc(:coarse; N=12, bvectest=true, nJac=10)
debug_disc(:fine; N=12)

debug_ellipse(:fine; N=12, bvectest=true) #precompstest=true

debug_annulus(:coarse; N=12, slptest=true)

debug_ellipseNh(:coarse; N=10, invtest=true, dlptest=true, slptest=true)

debug_peanut(:fine; N=10, bvectest=true)

debug_squircle(:fine; N=10, dlptest=true)

debug_kite(:fine; N=10, dlptest=true)

debug_bean(:fine; N=10, precompstest=true)

debug_star(:petal4; N=10, invtest=true, dlptest=true)

debug_star(:petal6; N=10, invtest=true, dlptest=true)