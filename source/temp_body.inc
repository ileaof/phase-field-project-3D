! ============================================================================
! temp_body.inc -- per-cell body of the heat-equation update (temperatura)
! ----------------------------------------------------------------------------
! Single source of truth, INCLUDEd verbatim from both the CPU worksharing loop
! and the GPU offload loop. TE..TB, PP, h_linha, FONTE, TNEW1 are listed PRIVATE
! on both enclosing directives; TP/cp/D/DELTA_T/delta_H/dx/dy/dz are read-only
! scalars (implicitly firstprivate on the device). TOLD is read-only, TNEW is
! written at the owned (i,j,k) only -> no race.
! ============================================================================

TE= Told(i+1,j,k)

TW= Told(i-1,j,k)

TN= Told(i,j+1,k)

TS= Told(i,j-1,k)

TT= Told(i,j,k+1)

TB= Told(i,j,k-1)


	   PP = phase_field_old(I,J,K)


!=======================================TERMO FONTE======================================


 h_linha = 30.d0*PP*PP*(1.d0-PP)*(1.d0-PP)           ! (ileao mesmo do fi) quintico; anterior: 2.d0*PP*(3.d0-2.d0*PP)-2.d0*PP*PP
 FONTE   = (delta_H/Cp)*h_linha*DPP(I,J,K)             ! (DELTAH/cp)*(DX*DY/Dtime)*Hlinha*(Phase_Field(I,J)-Phase_Field_Old(I,J))


!========================================================================================

TNEW1= D*DELTA_T*(((TE-2.d0*TP+TW)/dx**2.d0)+((TN-2.d0*TP+TS)/dy**2.d0)+((TT-2.d0*TP+TB)/dz**2.d0))+FONTE+TP



	 TNEW(I,J,K) = TNEW1
