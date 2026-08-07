! ============================================================================
! diag_body.inc -- per-cell body of the verification diagnostics
! ----------------------------------------------------------------------------
! Single source of truth, INCLUDEd from both the CPU and GPU reduction loops.
! phi_p, grad2, w_term are PRIVATE; L2_phi/solute_total/free_energy/phi_max/
! phi_min are reduction variables on both directives.
! ============================================================================
      phi_p = phase_field_new(i,j,k)
      L2_phi = L2_phi + phi_p*phi_p
      if (phi_p>phi_max) phi_max = phi_p
      if (phi_p<phi_min) phi_min = phi_p
      solute_total = solute_total + conc_soluto_new(i,j,k)

      w_term = w0*phi_p*phi_p*(1.d0-phi_p)*(1.d0-phi_p)
      if (i>1 .and. i<Nx .and. j>1 .and. j<Ny .and. k>1 .and. k<Nz) then
        grad2 = ((phase_field_new(i+1,j,k)-phase_field_new(i-1,j,k))*0.5d0/dx)**2 &
              + ((phase_field_new(i,j+1,k)-phase_field_new(i,j-1,k))*0.5d0/dy)**2 &
              + ((phase_field_new(i,j,k+1)-phase_field_new(i,j,k-1))*0.5d0/dz)**2
      else
        grad2 = 0.d0
      end if
      free_energy = free_energy + (0.5d0*E0*E0*grad2 + w_term)*dx*dy*dz
