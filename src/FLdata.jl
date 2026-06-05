module FLdata

import ..disc, ..ellipse, ..kite, ..annulus, ..squircle

import ..peanut, ..bean, ..star, ..ellipseNh, ..abstractdomain, ..refine!

#This file supplies the function f for which the FL problem will be solved
#and a exact solution of the problem, if one exists.

export domainbuild
export makediscfuex
export makeellipsefuex
export makekitefuex
export makeannulusfuex
export makesquirclefuex
export makepeanutfuex
export makebeanfuex
export makestarfuex
export makeellipseNhfuex

using SpecialFunctions: gamma, loggamma
using HypergeometricFunctions: pFq

function jacobi(x::Float64, n::Int, a::Float64, b::Float64)

   if n == 0
      return 1.0
   elseif n == 1
      return 0.5 * (a - b + (a + b + 2) * x)
   end

   p0 = 1.0
   p1 = 0.5 * (a - b + (a + b + 2) * x)
   p2 = 0.0

   for i = 1:(n-1)
      a1 = 2 * (i + 1) * (i + a + b + 1) * (2 * i + a + b)
      a2 = (2 * i + a + b + 1) * (a * a - b * b)
      a3 = (2 * i + a + b) * (2 * i + a + b + 1) * (2 * i + a + b + 2)
      a4 = 2 * (i + a) * (i + b) * (2 * i + a + b + 2)
      p2 = 1.0 / a1 * ((a2 + a3 * x) * p1 - a4 * p0)

      p0 = p1
      p1 = p2
   end

   return p2
end

@inline function jacobiC(n::Int, s::Float64)
   return exp(2s * log(2.0) + 2 * (loggamma(1.0 + s + n) - loggamma(n + 1.0)))
   #return (2.0^(2s) * gamma(1.0 + s + n)^2) / gamma(n + 1)^2
end

function makediscfuex(n::Int, s::Float64)

   C = jacobiC(n, s)

   @inline function f!(F, x, y)
      @inbounds for i in eachindex(F, x, y)
         r2 = x[i] * x[i] + y[i] * y[i]
         F[i] = C * jacobi(1.0 - 2.0 * r2, n, 0.0, s)
      end
      return nothing
   end

   @inline function uex(x::Float64, y::Float64)::Float64
      r2 = x * x + y * y
      return (1.0 - r2)^s * jacobi(1.0 - 2.0 * r2, n, 0.0, s)
   end

   return f!, uex
end

function makeellipsefuex(n::Int, s::Float64)

   a, b = 1.0, 2.0

   jo = 4 * pi * pFq((s + 1, 1 / 2), (1,), -3)

   j1 = 2 * pi * pFq((s + 1, 1 / 2), (2,), -3)

   j2 = (3 * pi / 2) * pFq((s + 1, 1 / 2), (3,), -3)

   c1 = s * (jo - j1 + 2 * (s - 1) * (j1 - j2))

   c2 = (2 * s + 1) * (s + 1) * jo - s * (4 * s + 1) * j1 + 2 * s * (s - 1) * j2

   C = 2^(2 * s - 1) * gamma(1 + s)^2 / π

   @inline function f!(F, x, y)
      @inbounds for i in eachindex(F, x, y)
         F[i] = C * (c1 * x[i]^2 - s * (jo - j1) + c2 * y[i]^2 / 4)
      end
      return nothing
   end

   @inline function uex(x::Float64, y::Float64)::Float64
      r2 = (x / a) * (x / a) + (y / b) * (y / b)
      return y * y * (1 - r2)^s / 4
   end

   return f!, uex

end

function makekitefuex(n::Int, s::Float64)

   C = jacobiC(n, s)

   @inline function f!(F, x, y)
      @inbounds for i in eachindex(F, x, y)
         yy = (y[i] / 1.5)^2
         r2 = (x[i] + 1.3 * yy)^2 + yy
         F[i] = C * jacobi(1.0 - 2.0 * r2, n, 0.0, s)
      end
      return nothing
   end

   uex = nothing

   return f!, uex
end

function makeannulusfuex(n::Int, s::Float64)

   C = jacobiC(n, s)

   @inline function annulus_ξ(x::Float64, y::Float64)::Float64
      r2 = x * x + y * y
      return (8.0 * r2 - 5.0) / 3.0
   end

   @inline function f!(F, x, y)
      @inbounds for i in eachindex(F, x, y)
         ξ = annulus_ξ(x[i], y[i])
         F[i] = C * jacobi(-ξ, n, 0.0, s)
      end
      return nothing
   end

   uex = nothing

   return f!, uex
end

function makesquirclefuex(n::Int, s::Float64, P::Int)

   C = jacobiC(n, s)

   @inline function f!(F, x, y)
      @inbounds for i in eachindex(F, x, y)
         r2 = x[i]^P + y[i]^P
         F[i] = C * jacobi(1.0 - 2.0 * r2, n, 0.0, s)
      end
      return nothing
   end

   uex = nothing

   return f!, uex
end

function makepeanutfuex(n::Int, s::Float64)

   C = jacobiC(n, s)

   @inline function f!(F, x, y)
      @inbounds for i in eachindex(F, x, y)
         y2 = y[i]^2
         ξ = ((x[i] - 1.0)^2 + y2) * ((x[i] + 1.0)^2 + y2)
         F[i] = C * jacobi(1.0 - ξ, n, 0.0, s)
      end
      return nothing
   end

   uex = nothing

   return f!, uex
end

function makebeanfuex(n::Int, s::Float64)

   C = jacobiC(n, s)

   @inline function beaneta(x::Float64, y::Float64)::Float64
      r2 = x * x + y * y
      return x^3 + y^3 - r2 * r2
   end

   function f!(F, x, y)
      @inbounds for i in eachindex(F, x, y)
         η = beaneta(x[i], y[i])
         F[i] = C * jacobi(1.0 - 10.0 * η, n, 0.0, s)
      end
      return nothing
   end

   uex = nothing

   return f!, uex
end

function makestarfuex(n::Int, s::Float64)

   R = 1.0
   Ps = 1.0

   C = jacobiC(n, s)

   # Smooth radial test coordinate.
   # h(θ) = 1 + Ps/2 + Ps*cos(Pθ)/2,
   # so hmax = 1 + Ps and Rmax = R*(1 + Ps).
   Rmax = R * (1.0 + Ps)
   invR = 1.0 / (Rmax * Rmax)

   @inline function f!(F, x, y)
      @inbounds for i in eachindex(F, x, y)
         r2 = (x[i] * x[i] + y[i] * y[i]) * invR
         F[i] = C * jacobi(1.0 - 2.0 * r2, n, 0.0, s)
      end
      return nothing
   end

   uex = nothing

   return f!, uex
end

function makeellipseNhfuex(n::Int, s::Float64)

   C = jacobiC(n, s)

   # Default ellipseNh outer ellipse:
   # center = (0,0), R1 = 1, R2 = 1, tht = 0.
   @inline function f!(F, x, y)
      @inbounds for i in eachindex(F, x, y)
         r2 = x[i] * x[i] + y[i] * y[i]
         F[i] = C * jacobi(1.0 - 2.0 * r2, n, 0.0, s)
      end
      return nothing
   end

   uex = nothing

   return f!, uex

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