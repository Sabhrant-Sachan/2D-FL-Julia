using Revise, FractionalLaplace2D

using FractionalLaplace2D.FLdata

# =============================================
# Run this once and save the exact solution
# =============================================
if ~isfile("kite_ex_Eg2.bin")
    s, p = 0.75, 8

    δ, δclsbd = 0.1, 0.01

    domf = FractionalLaplace2D.kite(b=[4, 5, 9, 9, 5, 4, 5, 5, 5, 5, 4, 4],
        a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 3, 3])

    FractionalLaplace2D.refine!(domf, 1, 2, [9, 12, 13, 16, 181, 184, 185, 188])

    f!, _, fv = makekitefuex(5, s)

    nₚᵣ, Nf = 256, 12

    prob = Problem(; N=Nf, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=5, (f!)=f!, uex=nothing, dom=domf)

    opts = Options(; save_Uapp_pre="kite_ex_Eg2.bin", solver=:direct, 
    cond_num=false, benchmark=false)

    coref = solveFL(prob; opts=opts)

    println(SolveView(prob, opts, coref))

    display("==Reference for Kite finished==")
end
# =============================================

s, p = 0.75, 8;

δ, δclsbd = 0.1, 0.01;

Uexf = FractionalLaplace2D.load_vec_bin("kite_ex_Eg2.bin");

domf = FractionalLaplace2D.kite(b=[4, 5, 9, 9, 5, 4, 5, 5, 5, 5, 4, 4],
        a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 3, 3]);

FractionalLaplace2D.refine!(domf, 1, 2, [9, 12, 13, 16, 181, 184, 185, 188]);

Nf = 12;

f!, _, fv = makekitefuex(5, s);

dobenchmark, docondnum = false, false;

# ======================1=======================

dom = FractionalLaplace2D.kite(b = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);

Anₚᵣ = [64, 64];

AN = [10, 12];

for i in 1:2

    nₚᵣ, N = Anₚᵣ[i], AN[i]

    Uex = FractionalLaplace2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=Uex, dom=dom);

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark )

    core_res = solveFL(prob; opts=opts);

    println(SolveView(prob, opts, core_res))

end

# =======================2======================

dom = FractionalLaplace2D.kite(b=[2, 1, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1],
    a=[1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1]);

nₚᵣ, N = 64, 12;

Uex = FractionalLaplace2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=1, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# ======================3=======================

dom = FractionalLaplace2D.kite(b=[3, 2, 4, 4, 2, 3, 2, 2, 2, 2, 2, 2],
    a=[2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2]);

nₚᵣ, N = 128, 12;

Uex = FractionalLaplace2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=2, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# ======================4=======================

dom = FractionalLaplace2D.kite(b=[3, 2, 5, 5, 2, 3, 2, 2, 2, 2, 2, 2],
    a=[2, 2, 3, 3, 2, 2, 2, 3, 3, 2, 2, 2]);

nₚᵣ, N = 128, 12;

Uex = FractionalLaplace2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=2, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# =======================5======================

dom = FractionalLaplace2D.kite(b=[4, 3, 7, 7, 3, 4, 3, 4, 4, 3, 3, 3],
    a=[3, 3, 4, 4, 3, 3, 3, 4, 4, 3, 2, 2]);

# nₚᵣ, N = 128, 12;
nₚᵣ, N = 256, 12;

Uex = FractionalLaplace2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=3, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))
# =====================6========================

dom = FractionalLaplace2D.kite(b=[4, 4, 8, 8, 4, 4, 4, 4, 4, 4, 4, 4],
    a=[3, 3, 7, 7, 3, 3, 3, 5, 5, 3, 2, 2]);

FractionalLaplace2D.refine!(dom, 1, 2, [5, 8, 9, 12, 153, 156, 157, 160]);

nₚᵣ, N = 256, 12;

Uex = FractionalLaplace2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=4, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# ======================7=======================

dom = FractionalLaplace2D.kite(b=[4, 4, 8, 8, 4, 4, 4, 4, 4, 4, 3, 3],
    a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 3, 3]);

FractionalLaplace2D.refine!(dom, 1, 2, [13, 16, 165, 168]);

nₚᵣ, N = 256, 12;

Uex = FractionalLaplace2D.build_ref_coarse(Uexf, domf, Nf, dom, N);

prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
    dₙₕ=5, (f!)=f!, uex=Uex, dom=dom);

opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark);

core_res = solveFL(prob; opts=opts);

println(SolveView(prob, opts, core_res))

# =============================================

display("==End of Kite==")
