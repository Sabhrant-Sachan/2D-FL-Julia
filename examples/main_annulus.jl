using Revise, FractionalLaplace2D

using FractionalLaplace2D.FLdata

# =============================================
# Run this once and save the exact solution
# =============================================
if ~isfile("annulus_ex_Eg1.bin")
    s, p = 0.75, 8

    δ, δclsbd = 0.1, 0.01

    domf = FractionalLaplace2D.annulus(b=[4, 5, 9, 9, 5, 4, 5, 5, 5, 5, 4, 4],
        a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 3, 3])

    FractionalLaplace2D.refine!(domf, 1, 2, [9, 12, 13, 16, 181, 184, 185, 188])

    f!, _, fv = makeannulusfuex(5, s)

    nₚᵣ, Nf = 256, 12

    prob = Problem(; N=Nf, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=5, (f!)=f!, uex=nothing, dom=domf)

    opts = Options(; save_Uapp_pre="annulus_ex_Eg1.bin", solver=:direct, 
    cond_num=false, benchmark=false)

    coref = solveFL(prob; opts=opts)

    println(SolveView(prob, opts, coref))

    display("==Reference for Kite finished==")
end
# =============================================
