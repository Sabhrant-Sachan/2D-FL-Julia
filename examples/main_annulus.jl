using Revise, Dates, FL2D

using FL2D.FLdata

dobenchmark, docondnum = false, false

s, p = 0.75, 4

δ, δ_near, δ_intp = 0.1, 0.15, 5e-3;

dom = FL2D.annulus(b=[2, 2, 2, 2, 2, 2, 2, 2])

#----------- Checking if the SLP is correct -----------
dp = FL2D.domprop(12, 0.1, 0.15, 5e-3, dom)

IV = FL2D.compress_vars(dom, dp, s, p; matrix_form=true);

β, βexact, bZ, exact, βerr, slperr = FL2D.SLPtest(dom, dp, IV);

FL2D.plotfunc(dp, dom, vec(bZ))
#-------------------------------------------------------

Anₚᵣ = [32, 32, 32, 32, 64, 64, 64, 128, 128];

AN = [3, 4, 5, 6, 7, 8, 9, 10, 12];

f!, uex, fv = makeannulusfuex(5, s);

for i in 1:9

    nₚᵣ, N = Anₚᵣ[i], AN[i]

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δ_near=δ_near, δ_intp=δ_intp, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

    core_res = solveFL(prob; opts=opts)

    println(SolveView(prob, opts, core_res))

end

# =============================================
# Run this once and save the exact solution
# =============================================
if ~isfile("annulus_ex_Eg1.bin")
    s, p = 0.75, 8

    δ, δclsbd = 0.1, 0.01

    domf = FL2D.annulus(b=[4, 5, 9, 9, 5, 4, 5, 5, 5, 5, 4, 4],
        a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 3, 3])

    FL2D.refine!(domf, 1, 2, [9, 12, 13, 16, 181, 184, 185, 188])

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
