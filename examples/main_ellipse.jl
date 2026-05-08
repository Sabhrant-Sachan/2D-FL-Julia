using Revise, FL2D

using FL2D.FLdata

# =============================================
s, p = 0.75, 4;

δ, δclsbd = 0.1, 0.01;

dom = FL2D.ellipse(b=[1, 1, 1, 1, 1], R1=1.0, R2=2.0, L1=2.0, L2=0.8);

Anₚᵣ = [32, 32, 32, 32, 64, 64, 64, 128, 128];

AN = [3, 4, 5, 6, 7, 8, 9, 10, 12];

f!, uex, fv = makeellipsefuex(2, s);

dobenchmark, docondnum = true, true;

for i in 1:9

    nₚᵣ, N = Anₚᵣ[i], AN[i]

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom);

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark )

    core_res = solveFL(prob; opts=opts);

    println(SolveView(prob, opts, core_res))

end

# =============================================
dom =  FL2D.ellipse(b=[2, 2, 2, 2, 2], a=[1, 1, 1, 1, 1], R1=1.0, R2=2.0, L1=2.0, L2=0.8);

Anₚᵣ = [128, 128];

AN = [10, 12];

for i in 1:2

    nₚᵣ, N = Anₚᵣ[i], AN[i]

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom);

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark )

    core_res = solveFL(prob; opts=opts);

    println(SolveView(prob, opts, core_res))

end

# =============================================
dom =  FL2D.ellipse(b=[4, 3, 4, 3, 3], a=[2, 2, 2, 2, 2], R1=1.0, R2=2.0, L1=2.0, L2=0.8);

Anₚᵣ = [128, 128];

AN = [10, 12];

for i in 1:2

    nₚᵣ, N = Anₚᵣ[i], AN[i]

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=2, (f!)=f!, uex=uex, dom=dom);

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark )

    core_res = solveFL(prob; opts=opts);

    println(SolveView(prob, opts, core_res))

end

# =============================================
dom =  FL2D.ellipse(b=[5, 4, 5, 4, 5],a=[3, 3, 3, 3, 4], R1=1.0, R2=2.0, L1=2.0, L2=0.8);

Anₚᵣ = [128, 256];

AN = [10, 12];

for i in 1:2

    nₚᵣ, N = Anₚᵣ[i], AN[i]

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=2, (f!)=f!, uex=uex, dom=dom);

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark )

    core_res = solveFL(prob; opts=opts);

    println(SolveView(prob, opts, core_res))

end

# =============================================
dom =  FL2D.ellipse(b=[6, 6, 6, 6, 6],a=[3, 3, 3, 3, 4], R1=1.0, R2=2.0, L1=2.0, L2=0.8);

Anₚᵣ = [128, 256];

AN = [10, 12];

for i in 1:2

    nₚᵣ, N = Anₚᵣ[i], AN[i]

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=3, (f!)=f!, uex=uex, dom=dom);

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark )

    core_res = solveFL(prob; opts=opts);

    println(SolveView(prob, opts, core_res))

end

display("==End of Ellipse==")

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