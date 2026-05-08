using Revise, FL2D, MAT

using FL2D.FLdata

# =============================================
# Run this once and save the exact solution
# =============================================
if ~isfile("kite_ex_Eg1.bin")
    s, p = 0.75, 8

    δ, δclsbd = 0.1, 0.01

    nₚᵣ, Nf = 256, 10

    domf = FL2D.kite(b=[4, 4, 5, 5, 4, 4, 4, 3, 3, 4, 4, 4],
        a=[4, 2, 4, 4, 2, 4, 2, 5, 5, 2, 2, 2])

    f!, _, fv = makekitefuex(0, s)

    prob = Problem(; N=Nf, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=3, (f!)=f!, uex=nothing, dom=domf)

    #save_Uapp_pre="kite_ex_Eg1.bin"

    opts = Options(; solver=:direct, cond_num=false, benchmark=false)

    coref = solveFL(prob; opts=opts)

    println(SolveView(prob, opts, coref))

    display("==Reference for Kite finished==")

    FL2D.save_uapp_bin("kite_ex_Eg1.bin", coref.Uapp)
end
# ======================0=======================

s, p = 0.75, 8;

δ, δclsbd = 0.1, 0.01;

Uexf = FL2D.load_vec_bin("kite_ex_Eg1.bin");

domf = FL2D.kite(b=[4, 4, 5, 5, 4, 4, 4, 3, 3, 4, 4, 4],
        a=[4, 2, 4, 4, 2, 4, 2, 5, 5, 2, 2, 2]);

Nf = 10;

f!, _, fv = makekitefuex(0, s);

# ======================1=======================

dom = FL2D.kite(b = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);

dobenchmark, docondnum = false, false;

nₚᵣ, N = 64, 10

Uexc = FL2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=1, (f!)=f!, uex=Uexc, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# =======================2======================

dom = FL2D.kite(b=[2, 1, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1],
    a=[1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1]);

nₚᵣ, N = 128, 10;

Uex = FL2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=1, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# ======================3=======================

dom = FL2D.kite(b=[3, 2, 4, 4, 2, 3, 2, 2, 2, 2, 2, 2],
    a=[2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2]);

nₚᵣ, N = 128, 10;

Uex = FL2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=2, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# ======================4=======================

dom = FL2D.kite(b=[3, 2, 4, 4, 2, 3, 2, 2, 2, 2, 2, 2],
    a=[2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2]);

FL2D.refine!(dom, 1, 2, [4, 6, 34, 36]);

nₚᵣ, N = 128, 10;

Uex = FL2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=2, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# ======================5=======================

dom = FL2D.kite(b=[4, 3, 5, 5, 3, 4, 3, 3, 3, 3, 3, 3],
    a=[3, 2, 3, 3, 2, 3, 2, 4, 4, 2, 2, 2]);

nₚᵣ, N = 256, 10;

Uex = FL2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=3, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# =============================================

display("==End of Kite==")