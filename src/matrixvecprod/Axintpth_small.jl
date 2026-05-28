function Axintpth_small!(v::SubArray{Float64}, kM::Int, IntS::Matrix{Float64},
   d::D, dp::domprop, s::Float64, IV::IVT) where {D<:abstractdomain,IVT}

   # --------- Unpack IV (NamedTuple) ---------
   (; IV1, IVr, IVbdth, IVt) = IV
   (; N, Np, Cs, M, Nb, Ni) = IV1
   (; nr, fwr, zx, zt, zy, Df, idctrg) = IVr
   (; Zx, Zy, DJ, KIr, Wrg, Grgtmp, Grg) = IVt
   (; coeffs) = IVbdth

   VOL_FAR = UInt8(0)

   # v is M*Np + Mbd*N by Np matrix.
   # size(v) == (M*Np + Mbd*N, Np)

   # ------------------------------------------------------------------
   # Regular volume quadrature grid for source patch kM
   # ------------------------------------------------------------------
   mapxy_Dmap!(Zx, Zy, DJ, d, zx, zy, kM)

   dfunc!(Df, d, kM, zt, s)

   @inbounds for j in 1:nr
      @inbounds for i in 1:nr
         Wrg[i, j] = fwr[i] * DJ[i, j] * Df[i, j] * fwr[j]
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
         @inbounds for j in 1:nr
            @inbounds for i in 1:nr
               dx = x1 - Zx[i, j]
               dy = x2 - Zy[i, j]
               r2 = dx * dx + dy * dy

               KIr[i, j] = Wrg[i, j] * expm1(-s * log(r2))
            end
         end

         mul!(Grgtmp, idctrg', KIr)
         mul!(Grg, Grgtmp, idctrg)

         @inbounds for jj in 1:Np
            v[row, jj] = Cs * Grg[jj]
         end
      end
   end

   # ------------------------------------------------------------------
   # Boundary target rows
   # ------------------------------------------------------------------
   Nt = Ni + Nb

   bcol = dp.pthgo[M+1] + Nb + (dp.invgo[M+kM] - dp.invgo[M+1])

   @inbounds for row in (Ni+1):Nt

      x1 = dp.tgtpts[1, row]
      x2 = dp.tgtpts[2, row]

      ib = row - Ni
      ibp = cld(ib, N)
      kt = d.kd[ibp]
      ip = ib - (ibp - 1) * N

      if kt == kM
         # ---------------- Boundary singular block ----------------
         col = dp.pthgo[M+1] + (ibp - 1) * N + ip - 1

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
         @inbounds for j in 1:nr
            @inbounds for i in 1:nr
               dx = x1 - Zx[i, j]
               dy = x2 - Zy[i, j]
               r2 = dx * dx + dy * dy

               KIr[i, j] = Wrg[i, j] * expm1(-s * log(r2))
            end
         end

         mul!(Grgtmp, idctrg', KIr)
         mul!(Grg, Grgtmp, idctrg)

         @inbounds for jj in 1:Np
            v[row, jj] = Cs * Grg[jj]
         end
      end
   end

   return nothing
end