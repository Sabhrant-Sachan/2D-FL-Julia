using Revise, FL2D

import FL2D.FLdata as FLdata

#----------------------------

d = FL2D.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5])

dp = FL2D.domprop(24, 0.1, 0.1, d)

FL2D.plotns(dp, d, 1)

FL2D.plotnsbd(dp, d, 14)

FL2D.plotprojbd(dp, d, 18)

d = FL2D.annulus(b = [4, 5, 4, 5, 4, 5, 4, 5])

dp = FL2D.domprop(12, 0.1, 0.1, d)

f!, uex, fv = FLdata.makediscfuex(5, 0.1);

FL2D.plotu(dp, d, uex)

FL2D.plotfunc(dp, d, f!)

IntSbdex = FL2D.precompsDLP(d, dp, 2; n=2048);

for n in (64, 128, 256, 512, 1024)
    IntSbd = FL2D.precompsDLP(d, dp, 2; n=n);
    Err = maximum(abs.(IntSbd .- IntSbdex))
    FL2D.@printf("Err = %.2e\n", Err)
end


IntSbd = FL2D.precompsDLP(d, dp, 2; n=256);

FL2D.testDLPint(d, dp, IntSbd; nreg=64, plot_err = true)

#----------------------------

d = FL2D.kite(b=[4, 5, 9, 9, 5, 4, 5, 5, 5, 5, 4, 4],
a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 3, 3])

FL2D.refine!(d, 1, 2, [9, 12, 13, 16, 181, 184, 185, 188])

# d = FL2D.kite(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
# a=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

FL2D.chk_map(d)

b = [4, 5, 4, 5, 4, 5, 4, 5];
a = [1, 2, 1, 2, 1, 2, 1, 2];

FL2D.@btime FL2D.annulus(b = $b, a=$a, R11=1/2,R21=1/4,R22=1/2);

d = FL2D.annulus(b = [4, 5, 4, 5, 4, 5, 4, 5], a=[1, 2, 1, 2, 1, 2, 1, 2],
              R11=1/2,R21=1/4,R22=1/2,tht2=pi/3)

            
FL2D.@btime FL2D.annulus(b = $b, a=$a,  R11=1/2,R21=1/4,R22=1/2,tht2=pi/3);

FL2D.draw(d,1)


d = FL2D.annulus(b = [3, 3, 3, 3, 3, 3, 3, 3])

FL2D.drawbd(d)

s, p = 0.9, 5;

δ, δclsbd = 0.1, 0.01;

N = 12; Np = N*N;

#d = FLdata.domainbuild()

dp = FL2D.domprop(N, δ, δclsbd, d)
 
FL2D.plotbm(dp,d)

FL2D.chkinvpts(dp,d)
#FL2D.draw(d)

#FL2D.plotdp(dp,d;label=:int)
#IntS = (s >= 0.5) ? FL2D.precompsH(d, dp, s, p) :
#       FL2D.precompsL(d, dp, s, p)
#Now performing convergence analysis!

size(dp.prepts,2)

LL = 6

M = dp.pthgo[7]

Iapp = Vector{Float64}(undef, Np);

ErrI = Matrix{Float64}(undef, LL, M);

for PI in 1:M

    Iex = FL2D.precompsH(d, dp, s, p, PI, 2^(LL + 5))

    for i in 1:LL
        Iapp = FL2D.precompsH(d, dp, s, p, PI, 2^(i + 4))

        E = maximum(abs.(Iex .- Iapp))

        ErrI[i, PI] = E
    end

    #println("Error for PI $(Printf.@sprintf("%d", PI))")
    # for i in 1:LL
    #     E = maximum(abs.(Iex .- Iapp[:, i]))
    #     #println("$(Printf.@sprintf("%d : %.2e", 2^(i + 4), E))")
    #     ErrI[i, PI] = E
    # end

    # if ErrI[3] > 5 * (1e-14)
    #     println("Error for PI $(Printf.@sprintf("%d", PI))")
    #     display(ErrI)
    # end
end

tol = 1e-12

dp.pthgo[6]

idx1 = findall(i -> ErrI[4, i] > tol, axes(ErrI, 2))

idx2 = findall(i -> ErrI[3, i] > tol && ErrI[4, i] > tol && ErrI[5, i] > tol, axes(ErrI, 2))

idx3 = findall(i -> ErrI[3, i] > tol && ErrI[4, i] > tol, axes(ErrI, 2))

ErrI[:, 2]

PI = 1395


s, p = 0.9999, 5000;

Iex = FL2D.precompsH(d, dp, s, p, PI, 2^(LL + 5));

for i in 1:LL
    Iapp = FL2D.precompsH(d, dp, s, p, PI, 2^(i + 4))

    E = maximum(abs.(Iex .- Iapp))

    println("$(Printf.@sprintf("%d : %.2e", 2^(i + 4), E))")
end


Iex = FL2D.precompsH8(d, dp, s, p, PI, 2^(LL + 5));

Err8 = Vector{Float64}(undef, 8);

println("n:     I1,1     I1,2     I2,1     I2,2     I3,1     I3,2     I4,1     I4,2")

for i in 1:LL
    Iapp = FL2D.precompsH8(d, dp, s, p, PI, 2^(i + 4))

    for j in 1:8
        Err8[j] = maximum(abs.(Iex[:,j] .- Iapp[:,j]))
    end
    println("$(Printf.@sprintf("%d : %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e", 2^(i + 4), 
    Err8[1],Err8[2],Err8[3],Err8[4],Err8[5],Err8[6],Err8[7],Err8[8]))")
end

#Iapp = FL2D.precompsH8(d, dp, 0.999, 50, PI, 2^(3 + 4));

#----------------------------


IntS = (s >= 0.5) ? FL2D.precompsH(d, dp, s, p) :
                    FL2D.precompsL(d, dp, s, p)

@inline function f!(F, x, y)
        fill!(F,1.0)
        return nothing
end

b = FL2D.bvec(d, dp, s, f!)

FL2D.plotfunc(dp,d,b)

function FLsoldebug(N, δ, δclsbd, dₙₕ, s, p, f!, domainbuild)
    d  = domainbuild()
    dp = FL2D.domprop(N, δ, δclsbd, d)

    IntS = (s >= 0.5) ? FL2D.precompsH(d, dp, s, p) :
                        FL2D.precompsL(d, dp, s, p)

    b = FL2D.bvec(d, dp, s, f!)

    δeff = dp.delclsbd
    IV = FL2D.compress_vars(d, dp.N, s, p, dₙₕ, δeff)
    (; N, Np, M, Mbd) = IV.IV1

    Lpn = M*Np
    Lp  = Lpn + Mbd*N
    A = zeros(Float64, Lp, Lp)

    for k in 1:M
        @views v = A[:, (1 + Np*(k-1)) : (Np*k)]
        if k in d.kd
            FL2D.Axbdpth!(v, k, IntS, d, dp, s, IV)
        else
            FL2D.Axintpth!(v, k, IntS, d, dp, s, IV)
        end
    end

    for k in 1:Mbd
        @views v = A[:, Lpn + N*(k-1) + 1 : Lpn + N*k]
        FL2D.Axbdop!(v, k, d, dp, s, IV)
    end

    Uapp = copy(b)
    Uapp, ch = gmres!(Uapp, A, b; reltol=3e-15, abstol=0.0, restart=min(450, Lp), log=true)

    return dp, d, Uapp, ch, A, b, IntS
end


s, p = 0.75, 4;

δ, δclsbd = 0.1, 0.01;

N, dₙₕ = 10, 2;

f!, uex, fv = FLdata.makediscfuex(2, s);

dp, d, Uapp, ch, A, b, IntS = FLsoldebug(N, δ, δclsbd, dₙₕ, s, p, f!, FLdata.domainbuild);

#Uapp = A\b;

#@benchmark  FLsol($N, $δ, $δclsbd, $s, $p, $f!, $domainbuild) evals=1 samples=10 seconds=1000

Np = N^2;
M = d.Npat;
Lₚ = M * Np + length(d.kd) * N;

uappv= Vector{Float64}(undef, M * Np);
uexv = Vector{Float64}(undef, M * Np);

for i in 1:M*Np

    #Compute uexv
    x = dp.tgtpts[1, i]
    y = dp.tgtpts[2, i]

    uexv[i] = uex(x, y)

    #Compute uappv
    ℓ = cld(i, Np)

    j = i - (ℓ - 1) * Np

    _, r = divrem(j - 1, N)

    j1 = r + 1

    uappv[i] = Uapp[i] * FL2D.dfunc(d, ℓ, 2 * sinpi((2j1 - 1) / (4N))^2, s)

end

#Compute errors:

err_u  = uappv .- uexv;

maxerr = maximum(abs.(err_u));
relmax = maxerr / maximum(abs.(uexv));
l2err  = norm(err_u, 2);
rmse   = sqrt(FL2D.mean(err_u.^2));

println("Vars = N^2 *M + N*Mbd = $(Printf.@sprintf("%.d", Lₚ))")
println("Iters = $(Printf.@sprintf("%.d", ch.iters))")
println("\n-----------Error-----------")
println("Max Error   : $(Printf.@sprintf("%.2e", maxerr))")
println("Max Rel err : $(Printf.@sprintf("%.2e", relmax))")
println("L^2 error   : $(Printf.@sprintf("%.2e", l2err))")
println("Root MSE err: $(Printf.@sprintf("%.2e", rmse))")

plotfunc(dp,d,uappv)
plotfunc(dp,d,uexv)
plotfunc(dp,d,abs.(err_u))
