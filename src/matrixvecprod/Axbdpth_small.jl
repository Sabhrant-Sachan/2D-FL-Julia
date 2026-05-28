function Axbdpth_small!(v::SubArray{Float64}, kM::Int, IntS::Matrix{Float64},
   d::D, dp::domprop, s::Float64, IV::IVT) where {D<:abstractdomain,IVT}

   # --------- Unpack IV (NamedTuple) ---------
   (; IV1, IVr, IVbdth, IVbt1, IVbt2) = IV
   (; N, Np, Cs, M, Nb, Ni) = IV1
   (; nr, fwr) = IVr
   (; zx2, zy2, coeffs, nbd) = IVbdth
   (; Zx₂, Zy₂, DJ₂, KIbd) = IVbt1
   (; mfw, B1bd, B2bd, Gbdtmp, Gbd, Wbd) = IVbt2

   VOL_FAR = UInt8(0)

   # v is M*Np + Mbd*N by Np matrix.
   # size(v) == (M*Np + Mbd*N, Np)

   # ------------------------------------------------------------------
   # Boundary-patch quadrature grid for source patch kM
   # ------------------------------------------------------------------
   mapxy_Dmap!(Zx₂, Zy₂, DJ₂, d, zx2, zy2, kM)

   hc = d.pths[kM].ck1 - d.pths[kM].ck0
   Cs₂ = Cs * hc^s

   @inbounds for j in 1:nr
      @inbounds for i in 1:nbd
         Wbd[i, j] = mfw[i] * DJ₂[i, j] * fwr[j]
      end
   end

   # ------------------------------------------------------------------
   # Interior target rows
   # ------------------------------------------------------------------
   colns = dp.pthgo[kM] + Np

   @inbounds for row in 1:Ni

      x1 = dp.tgtpts[1, row]
      x2 = dp.tgtpts[2, row]

      l = cld(row, Np)
      jloc = row - (l - 1) * Np

      if l == kM
         # ---------------- Singular/self block ----------------
         col = dp.pthgo[kM] + jloc - 1

         @views SI = IntS[:, col]
         @views IM = reshape(SI, N, N)

         @inbounds for jj in 1:Np
            q, r = divrem(jj - 1, N)
            j1 = r + 1
            j2 = q + 1

            @views c1 = coeffs[:, j1]
            @views c2 = coeffs[:, j2]

            v[row, jj] = Cs * dot(c1, IM, c2)
         end

      elseif dp.volmode[row, kM] != VOL_FAR
         # ---------------- Near/close block ----------------
         @views NSI = IntS[:, colns]
         @views IM = reshape(NSI, N, N)

         @inbounds for jj in 1:Np
            q, r = divrem(jj - 1, N)
            j1 = r + 1
            j2 = q + 1

            @views c1 = coeffs[:, j1]
            @views c2 = coeffs[:, j2]

            v[row, jj] = Cs * dot(c1, IM, c2)
         end

         colns += 1

      else
         # ---------------- Regular/far block ----------------
         @inbounds for jj in 1:nr
            @inbounds for ii in 1:nbd
               dx = x1 - Zx₂[ii, jj]
               dy = x2 - Zy₂[ii, jj]
               r2 = dx * dx + dy * dy

               KIbd[ii, jj] = Wbd[ii, jj] * expm1(-s * log(r2))
            end
         end

         mul!(Gbdtmp, B1bd', KIbd)
         mul!(Gbd, Gbdtmp, B2bd)

         @inbounds for jj in 1:Np
            v[row, jj] = Cs₂ * Gbd[jj]
         end
      end
   end

   # ------------------------------------------------------------------
   # Boundary target rows
   # ------------------------------------------------------------------
   Nt = Ni + Nb
   
   bcol = dp.pthgo[M+1] + Nb + (dp.invgo[M+kM] - dp.invgo[M+1])

   @inbounds for row in (Ni + 1):Nt

      x1 = dp.tgtpts[1, row]
      x2 = dp.tgtpts[2, row]

      k₀ = cld(row - Ni, N)

      l = d.kd[k₀]
      jloc = row - Ni - (k₀ - 1) * N

      if l == kM
         # ---------------- Boundary singular block ----------------
         col = dp.pthgo[M+1] + jloc - 1 + N * (k₀ - 1)

         @views SI = IntS[:, col]
         @views IM = reshape(SI, N, N)

         @inbounds for jj in 1:Np
            q, r = divrem(jj - 1, N)
            j1 = r + 1
            j2 = q + 1

            @views c1 = coeffs[:, j1]
            @views c2 = coeffs[:, j2]

            v[row, jj] = Cs * dot(c1, IM, c2)
         end

      elseif dp.volmode[row, kM] != VOL_FAR
         # ---------------- Boundary near/close block ----------------
         @views NSI = IntS[:, bcol]
         @views IM = reshape(NSI, N, N)

         @inbounds for jj in 1:Np
            q, r = divrem(jj - 1, N)
            j1 = r + 1
            j2 = q + 1

            @views c1 = coeffs[:, j1]
            @views c2 = coeffs[:, j2]

            v[row, jj] = Cs * dot(c1, IM, c2)
         end

         bcol += 1

      else
         # ---------------- Boundary regular/far block ----------------
         @inbounds for jj in 1:nr
            @inbounds for ii in 1:nbd
               dx = x1 - Zx₂[ii, jj]
               dy = x2 - Zy₂[ii, jj]
               r2 = dx * dx + dy * dy

               KIbd[ii, jj] = Wbd[ii, jj] * expm1(-s * log(r2))
            end
         end

         mul!(Gbdtmp, B1bd', KIbd)
         mul!(Gbd, Gbdtmp, B2bd)

         @inbounds for jj in 1:Np
            v[row, jj] = Cs₂ * Gbd[jj]
         end
      end
   end

   return nothing
end