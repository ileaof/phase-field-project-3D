Program  principal

use variaveis
use omp_lib

implicit none

double precision :: L2_phi, phi_max, phi_min, solute_total, free_energy

! ---- lightweight wall-clock profiling state (see internal sub print_timing) ----
! omp_get_wtime() returns real(8) wall-clock seconds and works identically for the
! CPU-OpenMP and GPU-offload builds (it measures the HOST wall clock, which is what
! "how long is my run taking" means in both cases). cpu_time() is deliberately NOT
! used: it reports summed CPU-thread time, which is meaningless once work runs on
! many threads or on the GPU. All timing state below is host-local scalars, updated
! only at output events (every WriteFileD_timeCounter steps) -> overhead negligible.
real(8) :: t_start          ! wall-clock at simulation start
real(8) :: t_prev_output    ! wall-clock at the previous output event
integer :: n_outputs        ! how many output events have been written so far
integer :: it_prev_output   ! iteration index at the previous output event

! ---- separate performance / results log file ------------------------------
! Every output event appends its full record (the ESCREVENDO RESULTADOS header,
! the iteration/time/stoptime lines, the L2/energy diagnostics, and the timing
! block) to this plain-text file, IN ADDITION to printing to the console. The
! file is opened in append mode and closed at each event so every record is
! flushed to disk immediately -- nothing is lost if the run is interrupted.
integer,          parameter :: iu_log       = 77
character(len=*), parameter :: perf_log_file = 'performance_log.txt'

call entrada
call W0_E0

! ---------------------------------------------------------------------------
! Runtime device selection. omp_get_num_devices() returns the number of usable
! OpenMP target (offload) devices. If >0 (an NVIDIA/Intel/AMD accelerator is
! present AND this binary was built with offload enabled), the kernels take
! their GPU path; otherwise use_gpu stays .false. and the existing CPU OpenMP
! implementation runs. A single executable serves both cases -- on a
! non-offload build the target directives simply fall back to host execution
! and omp_get_num_devices() returns 0, so we never enter the GPU branches.
! ---------------------------------------------------------------------------
if (omp_get_num_devices() > 0) then
    use_gpu = .true.
    write(*,'(/,A,I0,A,/)') '### GPU offload ENABLED -- ', omp_get_num_devices(), &
                            ' OpenMP target device(s) detected ###'
else
    use_gpu = .false.
    write(*,'(/,A,/)') '### No OpenMP target device -- running CPU OpenMP path ###'
end if

! Make the seven solver fields resident on the device for the WHOLE time loop.
! They are initialised on the host in entrada.f90, so map(to:) uploads those
! initial values once here; every kernel then reads/writes the device copies and
! nothing else crosses the bus except the tiny diagnostic samples and the field
! snapshots taken just before file output. This is the single most important
! performance decision -- it turns thousands of per-iteration array transfers
! into one upload at start and one release at the end.
if (use_gpu) then
!$omp target enter data map(to: phase_field_old, phase_field_new, &
!$omp                          conc_soluto_old,  conc_soluto_new,  &
!$omp                          TOLD, TNEW, DPP)
end if


time=0.d0

! Start the wall clock. Taken once, on the host, just before the first time step
! so it excludes the one-off entrada/W0_E0/target-enter-data setup cost.
t_start        = omp_get_wtime()
t_prev_output  = t_start
n_outputs      = 0
it_prev_output = 0

! Create (or overwrite) the performance log for this run and write a one-time
! header. Subsequent output events append to it (see print_timing).
open(iu_log, file=perf_log_file, status='replace', action='write')
write(iu_log,'(A)') '###########################################################'
write(iu_log,'(A)') '###  Phase-field 3D -- performance / results log        ###'
write(iu_log,'(A)') '###########################################################'
write(iu_log,'(A,I0,A,I0,A,I0)') 'Grid  Nx x Ny x Nz     : ', nx, ' x ', ny, ' x ', nz
if (use_gpu) then
   write(iu_log,'(A,I0,A)') 'Execution              : GPU offload (', &
                            omp_get_num_devices(), ' device(s))'
else
   write(iu_log,'(A,I0,A)') 'Execution              : CPU OpenMP (', &
                            omp_get_max_threads(), ' threads)'
end if
write(iu_log,'(A)') ''
close(iu_log)

Do it=1, it_max


            time=delta_T+time

            call temperatura
            !!call Calc_Temperature  ! Estava comentada
            !!TNEW = T_

            ! TP (temperatura isotermica do banho) agora e lido do entrada.dat em entrada.f90
            ! (antes: TP = 919.0d0 fixo aqui). Historico: !T0 - 19.84d0*DELTA_T*IT

            !========================================================================================
            !DIF_CONC=ABS(C_NEW-C0)
            !nyf=ny
            !IF(DIF_CONC>=1.D-02) THEN
            !NYF=NY+1
            !NY=NYF
            !END IF

            !WRITE(*,*)NYF,DIF_CONC

            !========================================================================================


            !IF (DELTA_T*it>=16.2D-02) THEN     !  7.920D-02) THEN
            !WRITE(*,*)'FIM DOS CALCULOS'
            !STOP
            !END IF


            call equacao_campo_fase

            call concentracao_soluto
            
            !IF(1.0d+13>=IT>=1) THEN
            call auxiliar
            !END IF


            if(E_R >= WriteFileD_timeCounter) then

            E_R = 0

                      WRITE(*,2000)
                      WRITE(*,2010)
                      WRITE(*,2000)
                      WRITE(*,2000)
                      WRITE(*,2020) IT
                      WRITE(*,2030) time
                      WRITE(*,2040) stoptime
                      WRITE(*,2050) WriteFileD_timeCounter*Delta_T
 !                     WRITE(*,*)'ESCREVENDO OS RESULTADOS ','  ','TETA_0====',ran_0,'  ','TETA_1====',ran_1
 !                     WRITE(*,*)'ESCREVENDO OS RESULTADOS ','  ','TETA_2====',ran_2,'  ','TETA_3====',ran_3
                      WRITE(*,2000)
            call diagnosticos(L2_phi, phi_max, phi_min, solute_total, free_energy)
            WRITE(*,2060) L2_phi, phi_max, phi_min, solute_total, free_energy
            ! saida.f90 writes the full phase_field_new/conc_soluto_new fields to
            ! disk from the HOST, so refresh the host copies from the device first.
            ! This runs only at write checkpoints (rare), not every iteration.
            if (use_gpu) then
!$omp target update from(phase_field_new, conc_soluto_new)
            end if
            call saida
            call print_timing        ! performance summary for this output event
            end if
            e_r = e_r + 1



            ! CRITICAL: this was "TNEW = TOLD" in the original -- backwards for advancing
            ! time. It only worked because temperatura.f90 used to write TOLD(I,J,K)=
            ! TNEW(I,J,K) in-place during its sweep, which made TOLD already equal to TNEW
            ! by this point (a harmless no-op copy). Now that temperatura.f90 has been
            ! converted to double-buffering (required for race-free OpenMP), TOLD is no
            ! longer advanced anywhere else, so this line must be "TOLD = TNEW" or the
            ! temperature field would never advance in time.
            ! End-of-iteration double-buffer swap: advance old<=new for all three
            ! fields. Runs entirely on whichever device owns the data (no transfer).
            if (use_gpu) then
            !$omp target teams distribute parallel do collapse(3) private(i,j,k)
            do k=1,nz
              do j=1,ny
                do i=1,nx
                  phase_field_old(i,j,k) = phase_field_new(i,j,k)
                  TOLD(i,j,k)            = TNEW(i,j,k)
                  conc_soluto_old(i,j,k) = conc_soluto_new(i,j,k)
                end do
              end do
            end do
            !$omp end target teams distribute parallel do
            else
            !$omp parallel default(none) shared(phase_field_old,phase_field_new,TOLD,TNEW,conc_soluto_old,conc_soluto_new,nx,ny,nz) private(i,j,k)
            !$omp do collapse(2) schedule(static)
            do k=1,nz
              do j=1,ny
                do i=1,nx
                  phase_field_old(i,j,k) = phase_field_new(i,j,k)
                  TOLD(i,j,k)            = TNEW(i,j,k)
                  conc_soluto_old(i,j,k) = conc_soluto_new(i,j,k)
                end do
              end do
            end do
            !$omp end do
            !$omp end parallel
            end if


            !if(IT>= 8000) then
            !go to 100
            !end if

            if ((time .GE. (stoptime)).AND.(stoptime .NE. 0.d0)) then
                      write(*,1000) stoptime
                      WRITE(*,2000)
                      WRITE(*,2010)
                      WRITE(*,2000)
                      WRITE(*,2000)
                      WRITE(*,2020) IT
                      WRITE(*,2030) time
                      WRITE(*,2040) stoptime
                      WRITE(*,2050) WriteFileD_timeCounter*Delta_T
 !                     WRITE(*,*)'ESCREVENDO OS RESULTADOS ','  ','TETA_0====',ran_0,'  ','TETA_1====',ran_1 
 !                     WRITE(*,*)'ESCREVENDO OS RESULTADOS ','  ','TETA_2====',ran_2,'  ','TETA_3====',ran_3 
                      WRITE(*,2000)
                      if(E_R < WriteFileD_timeCounter ) then
                        if (use_gpu) then
!$omp target update from(phase_field_new, conc_soluto_new)
                        end if
                        call saida
                      end if
                      call print_timing    ! final performance summary
                      exit
             endif
                
            


end do

!! 100 call saida   ! �ltima saida repetida by ileao ! must be included a Delta_T on the stoptime eq. above

!============== liberando os vetores da mem�ria alocada =====================

! Release the device-resident fields before freeing the host allocations.
! map(delete:) drops them from the device data environment (their final values
! were already pulled back to the host by the target update before the last saida).
if (use_gpu) then
!$omp target exit data map(delete: phase_field_old, phase_field_new, &
!$omp                             conc_soluto_old,  conc_soluto_new,  &
!$omp                             TOLD, TNEW, DPP)
end if

deallocate(phase_field_old,phase_field_new)
deallocate(xc,yc,zc,DPP,TOLD,TNEW)
deallocate(conc_soluto_new, conc_soluto_old)
deallocate( T_ )

!=============================================================================
1000  format (//,5X,'Program terminated by stop time condition = ',f14.7,' s',/)
2000  format ('=====================================================')
2010  format ('=============  ESCREVENDO RESULTADOS  ===============')
2015  format ('Parada por criterio de tempo total < tempo de escrita')
2020  format ('IT==> ', I8,' iteracoes')
2030  format ('Tempo total simulado ==> ',d12.4,' s')
2040  format ('Tempo de parada (stoptime) ==> ', d12.4,' s')
2060  format ('L2(phi)==>',e14.6,'  phi_max==>',e14.6,'  phi_min==>',e14.6,'  solute_total==>',e14.6,'  free_energy==>',e14.6)
2050  format ('Intervalo de escrita em arquivo  ==> ', d12.4,' s')

contains

! ===========================================================================
! print_timing -- concise wall-clock performance summary, emitted exactly once
! per output event (called immediately after each `call saida`). It is an
! INTERNAL subroutine, so it reads it/time/it_max/stoptime/use_gpu (module
! globals from `variaveis`) and the timing state (t_start, n_outputs, ...) by
! host association -- no argument plumbing needed.
!
! Runs entirely on the host, in serial code, OUTSIDE every parallel/target
! region, so it cannot interfere with parallel execution or introduce races.
! Cost is one omp_get_wtime() plus a few prints, once every
! WriteFileD_timeCounter steps -> completely negligible.
!
! CPU vs GPU detection: use_gpu was set once at startup from
! omp_get_num_devices() (>0 => GPU path taken). We re-query the OpenMP runtime
! here for the live device count and, on the GPU path, the default device id.
! ===========================================================================
subroutine print_timing
implicit none
real(8) :: t_now, elapsed_total, elapsed_interval
real(8) :: avg_step, avg_output, eta, pct
integer :: nthreads, ndev, devid, steps_int
integer :: units(2), nunits, ku, u, ios
logical :: log_ok

t_now            = omp_get_wtime()
elapsed_total    = t_now - t_start
elapsed_interval = t_now - t_prev_output
n_outputs        = n_outputs + 1
steps_int        = it - it_prev_output

nthreads = omp_get_max_threads()   ! configured CPU-thread count (valid on both paths)
ndev     = omp_get_num_devices()

avg_step   = elapsed_total / dble(max(it,1))
avg_output = elapsed_total / dble(max(n_outputs,1))

! Progress and ETA are driven by stoptime (the real termination criterion here;
! it_max is left huge). Fall back to it/it_max if stoptime is disabled (==0).
if (stoptime > 0.d0) then
   pct = 100.d0 * time / stoptime
else
   pct = 100.d0 * dble(it) / dble(max(it_max,1))
end if
if (pct > 0.d0) then
   eta = elapsed_total * (100.d0 - pct) / pct
else
   eta = 0.d0
end if

! Open the separate log file in APPEND mode. If it cannot be opened (log_ok=.false.)
! we still print to the console -- logging never blocks the simulation.
open(iu_log, file=perf_log_file, status='old', position='append', &
     action='write', iostat=ios)
log_ok = (ios == 0)

! Units to write the timing block to: always the console (6); plus the log file
! when it is open. Writing through this list avoids duplicating the format lines.
nunits = 1 ; units(1) = 6
if (log_ok) then ; nunits = 2 ; units(2) = iu_log ; end if

! ---- header + diagnostics: written ONLY to the log (the console already got
!      these from the main loop's WRITE(*,20xx) / diagnosticos calls) ----------
if (log_ok) then
   write(iu_log,'(A)') '====================================================='
   write(iu_log,'(A)') '=============  ESCREVENDO RESULTADOS  ==============='
   write(iu_log,'(A)') '====================================================='
   write(iu_log,'(A)') '====================================================='
   write(iu_log,'(A,I8,A)')    'IT==> ', it, ' iteracoes'
   write(iu_log,'(A,D12.4,A)') 'Tempo total simulado ==> ', time, ' s'
   write(iu_log,'(A,D12.4,A)') 'Tempo de parada (stoptime) ==> ', stoptime, ' s'
   write(iu_log,'(A,D12.4,A)') 'Intervalo de escrita em arquivo  ==> ', &
                               WriteFileD_timeCounter*Delta_T, ' s'
   write(iu_log,'(A)') '====================================================='
   write(iu_log,'(A,E14.6,A,E14.6,A,E14.6,A,E14.6,A,E14.6)') &
        'L2(phi)==>',       L2_phi,       '  phi_max==>',      phi_max, &
        '  phi_min==>',     phi_min,      '  solute_total==>', solute_total, &
        '  free_energy==>', free_energy
end if

! ---- timing block: written to the console AND the log, once each -------------
do ku = 1, nunits
   u = units(ku)
   write(u,'(A)')            '--------------------------------------------------------'
   write(u,'(A,I0)')         'Output Number          : ', n_outputs
   write(u,'(A,F0.4)')       'Simulation Time        : ', time
   write(u,'(A,I0)')         'Iteration              : ', it
   write(u,'(A)')            ''
   if (use_gpu) then
      devid = omp_get_default_device()
      write(u,'(A)')         'Execution Device       : GPU (OpenMP target offload)'
      write(u,'(A,I0)')      'GPU Device ID          : ', devid
      write(u,'(A,I0)')      'Available Devices      : ', ndev
   else
      write(u,'(A)')         'Execution Device       : CPU (OpenMP threads)'
      write(u,'(A,I0)')      'Available Devices      : ', ndev
   end if
   write(u,'(A,I0)')         'OpenMP Threads         : ', nthreads
   write(u,'(A)')            ''
   write(u,'(A,F0.2,A)')     'Elapsed Time           : ', elapsed_total,    ' s'
   write(u,'(A,F0.2,A)')     'Time Since Last Output : ', elapsed_interval, ' s'
   write(u,'(A)')            ''
   write(u,'(A,F0.5,A)')     'Average Time/Step      : ', avg_step,   ' s'
   write(u,'(A,F0.2,A)')     'Average Time/Output    : ', avg_output, ' s'
   write(u,'(A,I0)')         'Steps This Interval    : ', steps_int
   write(u,'(A,I0)')         'Total Outputs Written  : ', n_outputs
   write(u,'(A,F0.2,A)')     'Completed              : ', pct, ' %'
   write(u,'(A,F0.1,A)')     'Estimated Remaining    : ', eta, ' s'
   write(u,'(A)')            '--------------------------------------------------------'
end do

if (log_ok) close(iu_log)

! advance the interval markers for the next output event
t_prev_output  = t_now
it_prev_output = it

end subroutine print_timing

end program  principal