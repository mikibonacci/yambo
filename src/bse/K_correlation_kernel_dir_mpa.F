!
! License-Identifier: GPL
!
! Copyright (C) 2022 The Yambo Team
!
! Authors (see AUTHORS file for details): MB
!
#include<dev_defs.h>
! 
function K_correlation_kernel_dir_mpa(i_block,i_p,i_pmq, i_k, i_kmq, Ken, omega_s, eta,                             &
&                                 i_k_s,i_p_s,i_n_k,i_n_p,i_kmq_s,i_pmq_s,i_m_k,i_m_p,  &
&                                 i_kmq_t,i_pmq_t,i_k_sp_pol_n,i_p_sp_pol_n,i_k_sp_pol_m,i_p_sp_pol_m, &
&                                 iq_W,iq_W_s,ig_W,i_k_s_m1,iq_W_s_mq,ig_W_mq,i_kmq_s_m1,&
&                                 BS_n_g_W,O1,O2,O_times_W)
 !
 use pars,         ONLY:SP,cZERO,pi
 use BS,           ONLY:BS_blk
 use mpa_m,        ONLY:MPA_E_store,MPA_R_store,MPA_V_store
 use R_lattice,    ONLY:DEV_VAR(G_m_G),DEV_VAR(g_rot)
 use wrapper,      ONLY:Vstar_dot_V_gpu,M_by_V_gpu
 use electrons, ONLY:levels
 !
 implicit none
 !
 integer,     intent(in) :: i_block,i_p,i_pmq,                                            &
&                                 i_k_s,i_p_s,i_n_k,i_n_p,i_kmq_s,i_pmq_s,i_m_k,i_m_p, &
&                                 i_kmq_t,i_pmq_t,i_k_sp_pol_n,i_p_sp_pol_n,i_k_sp_pol_m,i_p_sp_pol_m, &
&                                 iq_W,iq_W_s,ig_W,i_k_s_m1,iq_W_s_mq,ig_W_mq,i_kmq_s_m1
 integer,     intent(in) :: BS_n_g_W,i_k,i_kmq
 real(SP),    intent(in) :: omega_s, eta
 complex(SP), intent(inout) DEV_ATTR :: O1(BS_n_g_W),O2(BS_n_g_W),O_times_W(BS_n_g_W)
 !
 complex(SP)  :: K_correlation_kernel_dir_mpa
 !
 integer      :: iO1,iO2,i_g1,i_g2,i_g3,i_np,iq_mpa,n_poles,nibz_mpa
 integer      :: ng1,ng2
 real(SP)     :: d_cv_prime,d_c_prime_v
 complex(SP)  :: z_cv_prime,z_c_prime_v,pole_factor
 complex(SP)  :: E_eff,R_eff
 real(SP)     :: pole_sign
 complex(SP), allocatable DEV_ATTR :: W_pole(:,:)
 complex(SP), pointer DEV_ATTR :: O_c_d_iO1(:)
 complex(SP), pointer DEV_ATTR :: O_c_d_iO2(:)
 logical      :: l_time_reversal_block
 !
 type(levels), intent(in) :: Ken

 K_correlation_kernel_dir_mpa=cZERO
 !
 if (i_k_sp_pol_n/=i_p_sp_pol_n .or. i_k_sp_pol_m/=i_p_sp_pol_m) then
   K_correlation_kernel_dir_mpa=cZERO
   return
 endif
 !
 ! The MPA database is indexed by the q sections written by io_MPA.  The
 ! symmetry and reciprocal-vector transformations needed by the collision
 ! contraction are handled below through iq_W_s/iq_W_s_mq and ig_W/ig_W_mq,
 ! exactly as in K_correlation_kernel_dir.  The MPA database does not
 ! write a second q block; if K_kernel selects it, we reconstruct it by
 ! transposing the stored IBZ pole matrices below.
 if (.not.allocated(MPA_E_store).or..not.allocated(MPA_R_store)) then
   call error('[BSE MPA] Pole coefficients are not loaded')
   K_correlation_kernel_dir_mpa=cZERO
   return
 endif
 ! The MPA database stores only the q points in the IBZ.  For the
 ! time-reversal sector, K_kernel uses the second BS_W block and shifts
 ! iq_W by q%nibz.  BS_W contains that second block explicitly, whereas
 ! the MPA arrays do not.  Map the shifted index back to the IBZ and
 ! transpose G,G' below, reproducing the BS_W construction.
 nibz_mpa=size(MPA_E_store)
 l_time_reversal_block=iq_W>nibz_mpa
 if (l_time_reversal_block) then
   iq_mpa=iq_W-nibz_mpa
 else
   iq_mpa=iq_W
 endif
 if (iq_mpa<1.or.iq_mpa>nibz_mpa.or. &
&    .not.allocated(MPA_E_store(iq_mpa)%blc).or. &
&    .not.allocated(MPA_R_store(iq_mpa)%blc)) then
   call error('[BSE MPA] Invalid q point in the pole database')
   K_correlation_kernel_dir_mpa=cZERO
   return
 endif
 n_poles=size(MPA_E_store(iq_mpa)%blc,3)
 if (n_poles<1) then
   call error('[BSE MPA] The MPA representation contains no poles')
   K_correlation_kernel_dir_mpa=cZERO
   return
 endif
 if (eta<=0._SP) call error('[BSE MPA] BSEDynEta must be positive')
 !
 ! The two denominators are crossed: the first combines c(k) with v'(p-q),
 ! while the second combines c'(p) with v(k-q).  Ken uses IBZ k indices;
 ! the symmetry labels were already used when the collision vectors O1/O2
 ! were generated.
 d_cv_prime=Ken%E(i_n_k,i_k,i_k_sp_pol_n)-Ken%E(i_m_p,i_pmq,i_p_sp_pol_m)
 d_c_prime_v=Ken%E(i_n_p,i_p,i_p_sp_pol_n)-Ken%E(i_m_k,i_kmq,i_k_sp_pol_m)
 z_cv_prime=cmplx(omega_s-d_cv_prime,eta,kind=SP)
 z_c_prime_v=cmplx(omega_s-d_c_prime_v,eta,kind=SP)
 !
 ! W_pole is a workspace.  The pole matrix is rebuilt for every pole,
 ! because its denominator depends on G,G' through E_n(G,G').
 allocate(W_pole(BS_n_g_W,BS_n_g_W))
 W_pole=cZERO
 !
 iO1 = BS_blk(i_block)%O_table(i_k_s,i_p_s,1,1,i_n_k,i_n_p,i_k_sp_pol_n)
 iO2 = BS_blk(i_block)%O_table(i_kmq_s,i_pmq_s,i_kmq_t,i_pmq_t,i_m_k,i_m_p,i_k_sp_pol_m)
 !
 O_c_d_iO1 => DEV_VAR(BS_blk(i_block)%O_c)(:,iO1)
 O_c_d_iO2 => DEV_VAR(BS_blk(i_block)%O_c)(:,iO2)
 !
 !DEV_ACC_DEBUG data present(O_c_d_iO1,O_c_d_iO2,G_m_G,g_rot,O1,O2)
 !DEV_ACC parallel loop
 !DEV_CUF kernel do(1) <<<*,*>>>
 !DEV_OMP parallel do default(shared), private(i_g1,i_g2,i_g3)
 do i_g1=1,BS_n_g_W
   i_g2=DEV_VAR(G_m_G)(DEV_VAR(g_rot)(i_g1,iq_W_s)   ,ig_W   )
   i_g3=DEV_VAR(G_m_G)(DEV_VAR(g_rot)(i_g1,iq_W_s_mq),ig_W_mq)
   O1(i_g1)=O_c_d_iO1(i_g2)
   O2(i_g1)=O_c_d_iO2(i_g3)
 enddo
 !DEV_OMP end parallel do
 !DEV_ACC_DEBUG end data
 ! MPA_V_store contains exactly the gamp produced by scatter_Gamp_gpu for
 ! this q; it is not read from the MPA database.  The bare identity*gamp
 ! term is added to the first pole matrix below, avoiding a second
 ! matrix-vector contraction.
 if (.not.allocated(MPA_V_store).or.iq_mpa>size(MPA_V_store).or. &
&    .not.allocated(MPA_V_store(iq_mpa)%blc)) then
   call error('[BSE MPA] Bare Coulomb prefactor is not retained')
   deallocate(W_pole)
   K_correlation_kernel_dir_mpa=cZERO
   return
 endif
 if (size(MPA_V_store(iq_mpa)%blc,1)<BS_n_g_W.or. &
&    size(MPA_V_store(iq_mpa)%blc,2)<BS_n_g_W) then
   call error('[BSE MPA] Retained gamp matrix is smaller than BS_n_g_W')
   deallocate(W_pole)
   K_correlation_kernel_dir_mpa=cZERO
   return
 endif
 do i_np=1,n_poles
   ! Build the screened pole contribution element by element.  MPA_R already
   ! contains the residue of the fitted X*gamp quantity, so no second gamp
   ! multiplication belongs here.
   ! The first matrix also receives delta_GG'*gamp, corresponding to the
   ! identity term in (X+I)*gamp.
   ! Each G,G' element is independent.  Keep this loop parallel in CPU
   ! builds; the DEV_ACC/DEV_CUF/DEV_OMPGPU markers are translated by the
   ! usual Yambo preprocessing for accelerator builds.
   !DEV_ACC_DEBUG data present(W_pole,MPA_E_store,MPA_R_store)
   !DEV_ACC parallel loop collapse(2)
   !DEV_CUF kernel do(2) <<<*,*>>>
   !$omp parallel do collapse(2), default(shared), &
   !$omp private(i_g1,i_g2,ng1,ng2,pole_sign,E_eff,R_eff,pole_factor)
   do i_g2=1,BS_n_g_W
     do i_g1=1,BS_n_g_W
       if (l_time_reversal_block) then
         ! BS_W(:,:,nibz+iq) is built from X(G',G,q), so the
         ! corresponding MPA pole matrix uses the exchanged indices.
         ng1=i_g2
         ng2=i_g1
       else
         ng1=i_g1
         ng2=i_g2
       endif
       ! MPA_R is the residue of X*gamp, because the fit input was
       ! X_par*isc%gamp.  The Coulomb prefactor must not be applied again.
       ! NB: here we do use the same pole convention as QP_mpa.F,
       ! positive Re and negative Im are enforced in mod_MPA, 
       pole_sign=-sign(1._SP,aimag(MPA_E_store(iq_mpa)%blc(ng1,ng2,i_np)))
       E_eff=pole_sign*MPA_E_store(iq_mpa)%blc(ng1,ng2,i_np)
       R_eff=pole_sign*MPA_R_store(iq_mpa)%blc(ng1,ng2,i_np)
       pole_factor=R_eff*( &
&        1._SP/(z_cv_prime-E_eff)+ &
&        1._SP/(z_c_prime_v-E_eff))
       W_pole(i_g1,i_g2)=pole_factor
       if (i_np==1.and.i_g1==i_g2) then
         W_pole(i_g1,i_g2)=W_pole(i_g1,i_g2)+ &
&          MPA_V_store(iq_mpa)%blc(i_g1,i_g1,1)
       endif
     enddo
   enddo
   !$omp end parallel do
   !DEV_ACC_DEBUG end data
   !
   call M_by_V_gpu('T',BS_n_g_W,W_pole,O1,O_times_W)
   K_correlation_kernel_dir_mpa=K_correlation_kernel_dir_mpa+ &
&       Vstar_dot_V_gpu(BS_n_g_W,O2,O_times_W)
 enddo
 deallocate(W_pole)
 K_correlation_kernel_dir_mpa=K_correlation_kernel_dir_mpa*4._SP*pi
 !
end function K_correlation_kernel_dir_mpa

!
! Prepare host MPA coefficients before entering the transition/OpenMP loop.
! All MPI ranks must call this helper. io_MPA uses global single-q buffers;
! it must not be called lazily from individual kernel-element evaluations.
! TODO:check carefully this routine.
!
!
subroutine K_mpa_load_if_needed(X,Xw,q,l_loaded)
 use X_m,          ONLY:X_t
 use frequency,    ONLY:w_samp,W_reset
 use R_lattice,    ONLY:bz_samp
 use mpa_m,        ONLY:MPA_E_store,MPA_R_store,MPA_V_store,MPA_store_free,MPA_E_par,MPA_R_par
 use matrix,       ONLY:PAR_matrix,MATRIX_reset,MATRIX_copy
 use IO_int,       ONLY:io_control
 use IO_m,         ONLY:OP_RD_CL,DUMP,REP,NONE,io_MULTIPOLE,IO_NO_DATABASE,deliver_IO_error_message
 use parallel_int, ONLY:PP_redux_wait
 use parallel_m,   ONLY:ncpu
#include<y_memory.h>
 type(X_t), intent(in) :: X
 type(w_samp), intent(in) :: Xw
 type(bz_samp), intent(in) :: q
 logical, intent(out), optional :: l_loaded
 type(X_t) :: X_read
 type(w_samp) :: W_read
 type(PAR_matrix), allocatable :: E_saved(:),R_saved(:)
 integer :: iq,ID,io_err,n_ready
 integer, external :: io_MPA
 logical :: ready,multipole_saved,load_ok

 if (present(l_loaded)) l_loaded=.false.
 !
 ! Allocation alone is not sufficient: every q must have complete matrices.
 ! For this first consumer, store(:) contains IBZ q points only; the caller
 ! handles transposed screening indices separately.
 !
 ready=.false.
 if (allocated(MPA_E_store).and.allocated(MPA_R_store)) then
   ! E/R are the objects loaded from mpa_ER.  gamp is deliberately checked
   ! separately: it is not part of that database and is regenerated by the
   ! caller with scatter_Gamp_gpu when needed.
   ready=size(MPA_E_store)==q%nibz.and.size(MPA_R_store)==q%nibz
   if (ready) then
     do iq=1,q%nibz
       if (.not.allocated(MPA_E_store(iq)%blc)) then
         ready=.false.
         exit
       endif
       if (.not.allocated(MPA_R_store(iq)%blc)) then
         ready=.false.
         exit
       endif
       if (any(lbound(MPA_E_store(iq)%blc)/=1).or. &
&          any(lbound(MPA_R_store(iq)%blc)/=1)) ready=.false.
       if (any(shape(MPA_E_store(iq)%blc)/=(/X%ng,X%ng,X%mpa_npoles/)).or. &
&          any(shape(MPA_R_store(iq)%blc)/=(/X%ng,X%ng,X%mpa_npoles/))) ready=.false.
     enddo
   endif
 endif
 n_ready=0
 if (ready) n_ready=1
 call PP_redux_wait(n_ready)
 if (n_ready==ncpu) then
   if (present(l_loaded)) l_loaded=.true.
   return
 endif
 !
 ! Keep the caller's screening/frequency metadata and scratch buffers intact.
 ! io_MPA reallocates Xw%p on reads, although the ER database has no samples.
 !
 X_read=X
 call W_reset(W_read)
 W_read%n_freqs=Xw%n_freqs
 multipole_saved=io_MULTIPOLE
 io_MULTIPOLE=.true.
 call MPA_store_free()
 call move_alloc(MPA_E_par,E_saved)
 call move_alloc(MPA_R_par,R_saved)
 allocate(MPA_E_par(1),MPA_R_par(1))
 !
 ! AAA: do we need a control first? 
 call io_control(ACTION=OP_RD_CL,MODE=DUMP,COM=REP,SEC=(/1/),ID=ID)
 io_err=io_MPA(X_read,W_read,ID)
 ! A missing mpa_ER database is a normal signal to the caller: it should
 ! perform a fresh fit.  Other I/O errors remain fatal because silently
 ! fitting from a corrupt or incompatible database would hide the problem.
 if (io_err==IO_NO_DATABASE) goto 100
 call deliver_IO_error_message(io_err,'BSE MPA/header',STOP_it=.TRUE.)
 if (io_err/=0) goto 100
 if (X_read%mpa_npoles/=X%mpa_npoles.or.X_read%mpa_npoles<1) then
   call error('[BSE MPA] Database pole count differs from the requested screening')
   goto 100
 endif
 !
 ! Replicate the complete matrices: transition blocks on a rank may need
 ! any screening q and any G pair, independently of response parallelism.
 !
 call MATRIX_init('SER',MPA_E_par(1),X%ng,X%mpa_npoles)
 call MATRIX_init('SER',MPA_R_par(1),X%ng,X%mpa_npoles)
 YAMBO_ALLOC(MPA_E_par(1)%blc,(X%ng,X%ng,X%mpa_npoles))
 YAMBO_ALLOC(MPA_R_par(1)%blc,(X%ng,X%ng,X%mpa_npoles))
 allocate(MPA_E_store(q%nibz),MPA_R_store(q%nibz))
 load_ok=.true.
 do iq=1,q%nibz
   call io_control(ACTION=OP_RD_CL,MODE=DUMP,COM=NONE,SEC=(/2*iq,2*iq+1/),ID=ID)
   io_err=io_MPA(X_read,W_read,ID)
   call deliver_IO_error_message(io_err,'BSE MPA/poles and residues',STOP_it=.TRUE.)
   if (io_err/=0) then
     call MPA_store_free()
     load_ok=.false.
     exit
   endif
   call MATRIX_copy(MPA_E_par(1),MPA_E_store(iq),.FALSE.)
   call MATRIX_copy(MPA_R_par(1),MPA_R_store(iq),.FALSE.)
   YAMBO_ALLOC(MPA_E_store(iq)%blc,(X%ng,X%ng,X%mpa_npoles))
   YAMBO_ALLOC(MPA_R_store(iq)%blc,(X%ng,X%ng,X%mpa_npoles))
   MPA_E_store(iq)%blc=MPA_E_par(1)%blc
   MPA_R_store(iq)%blc=MPA_R_par(1)%blc
 enddo
 if (present(l_loaded)) l_loaded=load_ok
 !
100 continue
 call MATRIX_reset(MPA_E_par(1))
 call MATRIX_reset(MPA_R_par(1))
 deallocate(MPA_E_par,MPA_R_par)
 call move_alloc(E_saved,MPA_E_par)
 call move_alloc(R_saved,MPA_R_par)
 call W_reset(W_read)
 io_MULTIPOLE=multipole_saved
end subroutine K_mpa_load_if_needed
