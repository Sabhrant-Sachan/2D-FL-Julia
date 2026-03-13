module FLdata

import ..disc, ..ellipse, ..kite, ..abstractdomain, ..refine!

#This file supplies the function f for which the FL problem will be solved
#and a exact solution of the problem, if one exists.

export domainbuild
export makediscfuex
export makeellipsefuex
export makekitefuex

using SpecialFunctions: gamma
using HypergeometricFunctions: pFq 

function makediscfuex(n::Int, s::Float64)
    nsgn = isodd(n) ? -1.0 : 1.0

    C = (2.0^(2s) * nsgn * gamma(1 + s + n)^2) / gamma(n+1)^2

    # Precompute coeffs for P_n^(s,0)(ξ) as [a0, a1, ..., an]
    coeffs = if n == 0
        (1.0,)

    elseif n == 1
        (0.5 * s, 0.5 * (s + 2.0))

    elseif n == 2
        a2 = (s + 3.0) * (s + 4.0) * 0.125
        a1 = (2.0 * s * (s + 3.0)) * 0.125
        a0 = (s * s - s - 4.0) * 0.125
        (a0, a1, a2)

    elseif n == 3
        a3 = (s^3 + 15.0 * s^2 + 74.0 * s + 120.0) / 48.0
        a2 = (s^3 + 9.0 * s^2 + 20.0 * s) / 16.0
        a1 = (s^3 + 3.0 * s^2 - 10.0 * s - 24.0) / 16.0
        a0 = (s^3 - 3.0 * s^2 - 16.0 * s) / 48.0
        (a0, a1, a2, a3)

    elseif n == 4
        a4 = (s^4 + 26.0 * s^3 + 251.0 * s^2 + 1066.0 * s + 1680.0) / 384.0
        a3 = s * (s^3 + 18.0 * s^2 + 107.0 * s + 210.0) / 96.0
        a2 = (s^4 + 10.0 * s^3 + 11.0 * s^2 - 118.0 * s - 240.0) / 64.0
        a1 = s * (s^3 + 2.0 * s^2 - 37.0 * s - 110.0) / 96.0
        a0 = (s^4 - 6.0 * s^3 - 37.0 * s^2 + 42.0 * s + 144.0) / 384.0
        (a0, a1, a2, a3, a4)

    elseif n == 5
        a5 = (s^5 + 40.0 * s^4 + 635.0 * s^3 + 5000.0 * s^2 + 19524.0 * s + 30240.0) / 3840.0
        a4 = s * (s^4 + 30.0 * s^3 + 335.0 * s^2 + 1650.0 * s + 3024.0) / 768.0
        a3 = (s^5 + 20.0 * s^4 + 115.0 * s^3 - 20.0 * s^2 - 1796.0 * s - 3360.0) / 384.0
        a2 = s * (s^4 + 10.0 * s^3 - 25.0 * s^2 - 490.0 * s - 1176.0) / 384.0
        a1 = (s^5 - 85.0 * s^3 - 240.0 * s^2 + 564.0 * s + 1440.0) / 768.0
        a0 = s * (s^4 - 10.0 * s^3 - 65.0 * s^2 + 250.0 * s + 1024.0) / 3840.0
        (a0, a1, a2, a3, a4, a5)

    else
        error("implemented n = 0..5")
    end

    @inline function polJac(ξ::Float64)
        # Horner: coeffs are (a0, a1, ..., an)
        p = coeffs[n+1]
        @inbounds for k in n:-1:1
            p = muladd(p, ξ, coeffs[k])   # p = p*ξ + a_{k-1}
        end
        return p
    end

    @inline function f!(F, x, y)
        @inbounds for i in eachindex(F, x, y)
            r2 = x[i] * x[i] + y[i] * y[i]
            F[i] = C * polJac(2.0 * r2 - 1.0)
        end
        return nothing
    end

    @inline function fv(x::Float64, y::Float64)::Float64
        r2 = x * x + y * y
        return C * polJac(2 * r2 - 1)
    end

    @inline function uex(x::Float64, y::Float64)::Float64
        r2 = x*x + y*y
        return (1.0 - r2)^s * nsgn * polJac(2.0*r2 - 1.0)
    end

    return f!, uex, fv
end

function makeellipsefuex(n::Int, s::Float64)

    #1st example
    # a, b = 1.5, 1.0

    # C₁ = a^(-(2 * s + 1)) * b * pFq((s + 1, 1 / 2), (1,), 1 - (b / a)^2)

    # C₂ = 2.0^(2s) * gamma(1 + s)^2 * C₁

    # @inline function f!(F, x, y)
    #     fill!(F, C₂)
    #     return nothing
    # end

    # @inline function fv(x::Float64, y::Float64)::Float64
    #     return C₂
    # end

    # @inline function uex(x::Float64, y::Float64)::Float64
    #     r2 = (x / a) * (x / a) + (y / b) * (y / b)
    #     return (1 - r2)^s
    # end

    # return f!, uex, fv

    #2nd example
    a, b = 1.0, 2.0

    jo = 4 * pi * pFq((s + 1, 1 / 2), (1,), -3)

    j1 = 2 * pi * pFq((s + 1, 1 / 2), (2,), -3)

    j2 = (3*pi/2) * pFq((s + 1, 1/2), (3,), -3)

    c1 = s * (jo - j1 + 2 * (s - 1) * (j1 - j2))

    c2 = (2 * s + 1) * (s + 1) * jo - s * (4 * s + 1) * j1 + 2 * s * (s - 1) * j2

    C = 2^(2 * s - 1) * gamma(1 + s)^2 / π

    @inline function f!(F, x, y)
        @inbounds for i in eachindex(F, x, y)
            F[i] = C * (c1 * x[i]^2 - s * (jo - j1) + c2 * y[i]^2 / 4)
        end
        return nothing
    end

    @inline function fv(x::Float64, y::Float64)::Float64
        return C * (c1 * x^2 - s * (jo - j1) + c2 * y^2 / 4)
    end

    @inline function uex(x::Float64, y::Float64)::Float64
        r2 = (x / a) * (x / a) + (y / b) * (y / b)
        return y * y * (1 - r2)^s / 4
    end

    return f!, uex, fv  

end

function makekitefuex(n::Int, s::Float64)

   nsgn = isodd(n) ? -1.0 : 1.0

    C = (2.0^(2s) * nsgn * gamma(1 + s + n)^2) / gamma(n+1)^2

    # Precompute coeffs for P_n^(s,0)(ξ) as [a0, a1, ..., an]
    coeffs = if n == 0
        (1.0,)

    elseif n == 1
        (0.5 * s, 0.5 * (s + 2.0))

    elseif n == 2
        a2 = (s + 3.0) * (s + 4.0) * 0.125
        a1 = (2.0 * s * (s + 3.0)) * 0.125
        a0 = (s * s - s - 4.0) * 0.125
        (a0, a1, a2)

    elseif n == 3
        a3 = (s^3 + 15.0 * s^2 + 74.0 * s + 120.0) / 48.0
        a2 = (s^3 + 9.0 * s^2 + 20.0 * s) / 16.0
        a1 = (s^3 + 3.0 * s^2 - 10.0 * s - 24.0) / 16.0
        a0 = (s^3 - 3.0 * s^2 - 16.0 * s) / 48.0
        (a0, a1, a2, a3)

    elseif n == 4
        a4 = (s^4 + 26.0 * s^3 + 251.0 * s^2 + 1066.0 * s + 1680.0) / 384.0
        a3 = s * (s^3 + 18.0 * s^2 + 107.0 * s + 210.0) / 96.0
        a2 = (s^4 + 10.0 * s^3 + 11.0 * s^2 - 118.0 * s - 240.0) / 64.0
        a1 = s * (s^3 + 2.0 * s^2 - 37.0 * s - 110.0) / 96.0
        a0 = (s^4 - 6.0 * s^3 - 37.0 * s^2 + 42.0 * s + 144.0) / 384.0
        (a0, a1, a2, a3, a4)

    elseif n == 5
        a5 = (s^5 + 40.0 * s^4 + 635.0 * s^3 + 5000.0 * s^2 + 19524.0 * s + 30240.0) / 3840.0
        a4 = s * (s^4 + 30.0 * s^3 + 335.0 * s^2 + 1650.0 * s + 3024.0) / 768.0
        a3 = (s^5 + 20.0 * s^4 + 115.0 * s^3 - 20.0 * s^2 - 1796.0 * s - 3360.0) / 384.0
        a2 = s * (s^4 + 10.0 * s^3 - 25.0 * s^2 - 490.0 * s - 1176.0) / 384.0
        a1 = (s^5 - 85.0 * s^3 - 240.0 * s^2 + 564.0 * s + 1440.0) / 768.0
        a0 = s * (s^4 - 10.0 * s^3 - 65.0 * s^2 + 250.0 * s + 1024.0) / 3840.0
        (a0, a1, a2, a3, a4, a5)

    else
        error("implemented n = 0..5")
    end

    @inline function polJac(ξ::Float64)
        # Horner: coeffs are (a0, a1, ..., an)
        p = coeffs[n+1]
        @inbounds for k in n:-1:1
            p = muladd(p, ξ, coeffs[k])   # p = p*ξ + a_{k-1}
        end
        return p
    end

    @inline function f!(F, x, y)
        @inbounds for i in eachindex(F, x, y)
            yy = (y[i] / 1.5)^2
            r2 = (x[i] + 1.3 * yy)^2 + yy
            F[i] = C * polJac(2.0 * r2 - 1.0)
        end
        return nothing
    end

    @inline function fv(x::Float64, y::Float64)::Float64
        yy = (y / 1.5)^2
        r2 = (x + 1.3 * yy)^2 + yy
        return C * polJac(2 * r2 - 1)
    end

    uex = nothing

    return f!, uex, fv

end

end 

#A closure is a function bundled together with the variables it needs 
#from the scope where it was created, stored as fields in a callable object.
#Conceptually, Julia is doing something very close to:

# struct AnonymousF!
#     C::Float64
#     coeffs::NTuple{N,Float64}
# end

# function (f::AnonymousF!)(F, x, y)
#     # body of f!, using f.C and f.coeffs
# end

#When makediscfuex is called, Julia computes C, coeffs and then
#allocates a small immutable struct holding C and coeffs and
#returns a callable instance of that struct. 
#Closures are fast if:
#    1.captured variables are local
#    2.captured variables are type-stable
#    3.closures are created once, used many times
#Closures hurt only when:
#    1.they capture globals
#    2.they capture "Any" type 
#    3.they are created repeatedly in hot loops