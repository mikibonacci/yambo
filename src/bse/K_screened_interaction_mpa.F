!
! License-Identifier: GPL
!
! Copyright (C) 2018 The Yambo Team
!
! Authors (see AUTHORS file for details): MB
!
#include<dev_defs.h>
!
subroutine K_screened_interaction_mpa(X,Xw,q,io_X_err)
 !
 ! MPA poles/residues for BSE-screened interaction.
 !
 ! Focus of this routine (for now): fit X(G,G',q,w) -> (E,R) per iq_W
 ! and write poles to the mpa_ER database. The dynamical BSE kernel
 ! (K_correlation_kernel_mpa) is NOT built here.
 !
 ! Two modes:
 !   streaming (default, low memory): single MPA_E/R_par(1) working buffer,
 !     fitted and written iq-by-iq via io_MPA (OP_APP_CL).
 !   retention (optional): additionally keep all iq in memory via
 !     MPA_E/R_store(iq) (array of PAR_matrix).
 !
 use pars,         ONLY:SP
 use stderr,       ONLY:intc
 use com,          ONLY:exp_user,msg
 use frequency,    ONLY:w_samp
 use drivers,      ONLY:l_td_hf,l_tddft,l_rim_w
 use collision_el, ONLY:elemental_collision,elemental_collision_free
  use R_lattice,    ONLY:bz_samp,RIM_W,RIM_W_ng,RIM_W_is_diagonal,RIM_qpg,RIM_anisotropy
 use D_lattice,    ONLY:i_time_rev,i_space_inv
 use ALLOC,        ONLY:X_ALLOC_elemental
 use X_m,          ONLY:X_t,X_par
 use BS,           ONLY:BS_n_g_W,BS_W_is_diagonal,BS_W
 use IO_int,       ONLY:io_control,IO_and_Messaging_switch
 use IO_m,         ONLY:manage_action,OP_RD_CL,REP,VERIFY,NONE,&
&                       OP_RD,RD_CL_IF_END,deliver_IO_error_message,&
&                       io_RESPONSE,io_MULTIPOLE,OP_WR_CL,OP_APP_CL
 use devxlib,      ONLY:devxlib_memcpy_d2h
 use gpu_m,        ONLY:have_gpu
 use mpa_m,        ONLY:mpa_RE_solver,MPA_freqs,MPA_E_par,MPA_R_par,MPA_free,&
&                       MPA_E_store,MPA_R_store,MPA_V_store,MPA_store_free
 use timing_m,     ONLY:timing
 use parser_m,     ONLY:parser
 use parallel_m,   ONLY:PAR_COM_X_WORLD,PAR_COM_WORLD,PAR_COM_RL_INDEX,PAR_COM_Q_INDEX,&
&                       PP_indexes,PP_indexes_reset
 use parallel_int, ONLY:PP_wait,PP_redux_wait,PARALLEL_index
 use matrix,       ONLY:MATRIX_reset,MATRIX_copy
 !
#include<y_memory.h>
 !
 type(bz_samp) ::q
 type(X_t)     ::X
 type(w_samp)  ::Xw
 !
 integer,     intent(out) :: io_X_err
 external :: K_mpa_load_if_needed
 !
 type(elemental_collision), target ::isc
 integer    :: nibz
 !
 integer    :: BS_W_size2,i_g1,i_g1_mat,i_g2,i_g2_mat,iq_W
 logical    :: W_is_unscreened,l_RIM_W_g
 !
 ! I/O
 !
 integer           :: XID,ID2,IO_ACT,io_err
 integer, external :: io_X,io_MPA
 logical           :: condition
 !
 ! MPA variables
 !
 integer           :: X_rows(2),X_cols(2),ig_first,ig_last
 logical           :: MPred,l_para_G_MPA,l_mpa_loaded,l_user_keep_MPA
 real(SP)          :: MPred_rate,PPcond_rate,PPcond_Qrate,&
&                     MP_err,MP_Qerr,cond_num(2),cond_numQ(2)
 complex(SP), allocatable :: MPA_Xo(:)
 type(PP_indexes)  :: px
 !
 ! The dynamical BSE kernel consumes E/R and gamp after this routine returns,
 ! so the fit path retains all q blocks in memory.  The legacy
 ! MPA_on_the_fly parser is still read for compatibility and reported when
 ! it requests the old streaming-only behavior.
 !
 logical :: l_keep_MPA_in_memory
 !
 ! 2021/02/11 to be modified for functionals with screening
 !
 W_is_unscreened=l_td_hf.or.l_tddft
 !
 if (W_is_unscreened) then
   io_X_err=0
   return
 endif
 !
 if (X%ng<BS_n_g_W) then
   if (     exp_user) call warning("Using screening with less G than Weh interaction")
   if (.not.exp_user) call   error("Using screening with less G than Weh interaction!")
 endif
 !
 call elemental_collision_free(isc)
 !
 ! First try to reuse a compatible mpa_ER database.  This must happen before
 ! reading the response database with io_X: if E/R are already available,
 ! there is no reason to repeat the expensive X -> MPA fit.
 !
 call K_mpa_load_if_needed(X,Xw,q,l_mpa_loaded)
 if (l_mpa_loaded) then
   ! The database contains E/R, but never gamp.  Rebuild only gamp here,
   ! once per IBZ q point, using the same collision setup as the fit path.
   call X_ALLOC_elemental("X")
   if (BS_W_is_diagonal) call error("[MPA] BS_W_is_diagonal=.TRUE. not implemented with MPA")
   BS_W_size2=BS_n_g_W
   if (BS_W_is_diagonal) BS_W_size2=1
   if (allocated(MPA_V_store)) then
     do iq_W=1,size(MPA_V_store)
       call MATRIX_reset(MPA_V_store(iq_W))
     enddo
     deallocate(MPA_V_store)
   endif
   allocate(MPA_V_store(q%nibz))
   YAMBO_ALLOC(isc%gamp,(BS_n_g_W,BS_W_size2))
   if (have_gpu) then
     YAMBO_ALLOC_GPU(DEV_VAR(isc%gamp),(BS_n_g_W,BS_W_size2))
   endif
   do iq_W=1,q%nibz
     isc%qs(2)=iq_W
     call scatter_Gamp_gpu(isc,'c')
     if (have_gpu) call devxlib_memcpy_d2h(isc%gamp,DEV_VAR(isc%gamp))
     call MATRIX_init("SER",MPA_V_store(iq_W),BS_n_g_W,1)
     YAMBO_ALLOC(MPA_V_store(iq_W)%blc,(BS_n_g_W,BS_W_size2,1))
     MPA_V_store(iq_W)%blc(:,:,1)=isc%gamp(:BS_n_g_W,:BS_W_size2)
   enddo
   ! Match the static screened-interaction head correction.  The RIM part
   ! of the dynamical fit is already represented by RIM_W/2 below; only the
   ! frequency-independent identity term needs this correction.
   if (RIM_anisotropy/=0._SP.and.allocated(RIM_qpg)) then
     MPA_V_store(1)%blc(1,1,1)=MPA_V_store(1)%blc(1,1,1)+ &
&                                cmplx((RIM_qpg(1,1,1)-RIM_anisotropy)/2._SP, &
&                                      0._SP,kind=SP)
   endif
   call elemental_collision_free(isc)
   io_X_err=0
   return
 endif
 !
 ! No compatible mpa_ER database was found.  Only now verify/open the
 ! response database needed by the fallback X -> MPA fit.
 if(.not.l_rim_w) call section('+','Screened interaction header I/O')
 if(     l_rim_w) call section('=','Screened interaction header I/O')
 !===================================================
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2/),MODE=VERIFY,ID=XID)
 io_X_err=io_X(X,Xw,XID)
 if(io_X_err<0) return
 !
 ! Allocation of X_par (single-q working buffer, cf. io_X.F + QP_mpa.F).
 ! X_par has size 1; io_X fills X_par(1) for the requested iq.
 !
  call X_ALLOC_elemental("X")
  if (.not.io_RESPONSE) call error("[MPA] .not. io_RESPONSE untested with MPA")
  if (BS_W_is_diagonal) call error("[MPA] BS_W_is_diagonal=.TRUE. not implemented with MPA")
 if (.not.allocated(X_par).and.io_RESPONSE) allocate(X_par(1))
 if (io_RESPONSE) call X_ALLOC_parallel(X_par(1),X%ng,2*X%mpa_npoles,"X")
 !
 X_rows=X_par(1)%rows
 X_cols=X_par(1)%cols
 !
 ! Single-buffer MPA workspace (IO-compatible: io_MPA reads MPA_E/R_par(1)).
 !
 if (.not.allocated(MPA_E_par)) allocate(MPA_E_par(1))
 if (.not.allocated(MPA_R_par)) allocate(MPA_R_par(1))
 call MPA_ALLOC_parallel(MPA_E_par(1),MPA_R_par(1),X%ng,X%mpa_npoles,"CPU_ONLY")
 !
 ! Optional retention storage
 !
 ! The dynamical kernel consumes every q point after this routine returns.
 ! Therefore E/R and gamp must be retained even when the historical
 ! MPA_on_the_fly option is disabled; otherwise the later transition loop
 ! would have no complete in-memory representation of W(omega).
 l_keep_MPA_in_memory=.true.
 call parser('MPA_on_the_fly',l_user_keep_MPA)
 if (.not.l_user_keep_MPA) call msg('r',' :: [MPA] retaining E/R/gamp for the dynamical BSE kernel')
 ! The fit below replaces the screening; discard any previously retained fit.
 call MPA_store_free()
 if (l_keep_MPA_in_memory) then
   allocate(MPA_E_store(q%nibz))
   allocate(MPA_R_store(q%nibz))
   allocate(MPA_V_store(q%nibz))
 endif
 !
 BS_W_size2=BS_n_g_W
 if (BS_W_is_diagonal) BS_W_size2=1
 !
 YAMBO_ALLOC(isc%gamp,(BS_n_g_W,BS_W_size2))
 if (have_gpu) then
   YAMBO_ALLOC_GPU(DEV_VAR(isc%gamp),(BS_n_g_W,BS_W_size2))
 endif
 !
 nibz=q%nibz
 !
 call io_control(ACTION=OP_RD,COM=NONE,SEC=(/1/),ID=XID)
 io_X_err=io_X(X,Xw,XID)
 if (io_X_err<0) then
   call warning('Reading of io_X header failed despite previous check.')
   return
 endif
 !
 ! MPA DB: header + per-iq appends (QP_mpa.F:296-301,441-442 pattern).
 !
 ! Parallel-I/O note: the mpa_ER file MUST be opened collectively in
 ! parallel mode (nf90_create_par). Hence a REAL communicator covering all
 ! ranks executing this replicated q-loop is passed (PAR_COM_WORLD, which
 ! coincides with PAR_COM_X_WORLD here), with DO_IT=.TRUE. on every rank:
 ! all ranks define (hence the 'independent' par-access in io_MPA is legal)
 ! and all ranks write the identical post-redux buffers. A NULL comm would
 ! fall back to a serial open, on which nf90_var_par_access fails with
 ! "Parallel operation on file opened for non-parallel access".
 ! `condition` below only gates messaging, not writing.
 !
 condition=PAR_COM_X_WORLD%CPU_id==0
 call parser('mpERdb',io_MULTIPOLE)
 !
 if (io_MULTIPOLE) then
   call io_control(ACTION=OP_WR_CL,COM=REP,SEC=(/1/),ID=ID2,COMM=PAR_COM_WORLD,DO_IT=.TRUE.)
   io_err=io_MPA(X,Xw,ID2)
   call IO_and_Messaging_switch("SAVE")
   call IO_and_Messaging_switch("+io_out",CONDITION=condition)
 endif
 !
 call PP_wait()
 !
 ! Optional extra G-distribution of the fit (as in QP_mpa.F:262-273).
 ! Incompatible with RL- or Q-parallelism of the self-energy; here the
 ! BSE q-loop is serial so we only guard on RL parallelism.
 ! Here below we are basically slicing 
 !
 call PP_indexes_reset(px)
 l_para_G_MPA=.not.(PAR_COM_RL_INDEX%n_cpu>1.or.PAR_COM_Q_INDEX%n_cpu>1)
  if (l_para_G_MPA) then
    call PARALLEL_index(px,uplim=(/BS_n_g_W/), low_range=(/1/),&
&                      COMM=PAR_COM_X_WORLD,CONSECUTIVE=.true.)
    ig_first=px%first_of_1D(PAR_COM_X_WORLD%cpu_id+1)
    ig_last =px%last_of_1D(PAR_COM_X_WORLD%cpu_id+1)
  else
    ig_first=1
    ig_last =BS_n_g_W
  endif
 !
 do iq_W=1,q%nibz
   !
   isc%qs(2)=iq_W
   call scatter_Gamp_gpu(isc,'c')
   if (have_gpu) call devxlib_memcpy_d2h(isc%gamp,DEV_VAR(isc%gamp))
   ! Keep the bare Coulomb prefactor beside the fitted coefficients.  It is
   ! an in-memory object only: io_MPA writes E/R, while gamp depends on the
   ! current Coulomb cutoff and q-grid and is regenerated here on every run.
   if (l_keep_MPA_in_memory) then
     call MATRIX_init("SER",MPA_V_store(iq_W),BS_n_g_W,1)
     YAMBO_ALLOC(MPA_V_store(iq_W)%blc,(BS_n_g_W,BS_W_size2,1))
     MPA_V_store(iq_W)%blc(:,:,1)=isc%gamp(:BS_n_g_W,:BS_W_size2)
     if (iq_W==1.and.RIM_anisotropy/=0._SP.and.allocated(RIM_qpg)) then
       ! The RIM pole fit already supplies RIM_W/2.  Add only the same
       ! static bare-head correction applied in K_screened_interaction.F.
       MPA_V_store(iq_W)%blc(1,1,1)=MPA_V_store(iq_W)%blc(1,1,1)+ &
&                                    cmplx((RIM_qpg(1,1,1)-RIM_anisotropy)/2._SP, &
&                                          0._SP,kind=SP)
     endif
   endif
   !
   IO_ACT=manage_action(RD_CL_IF_END,iq_W,1,q%nibz)
   call io_control(ACTION=IO_ACT,COM=NONE,SEC=(/2*iq_W,2*iq_W+1/),ID=XID)
   io_X_err=io_X(X,Xw,XID)
   call deliver_IO_error_message(io_X_err,'MPA/Xm',STOP_it=.TRUE.)
   !
   ! Sampling frequencies for this run (auto-allocates on first iq,
   ! reallocated only if n_freqs changes)
   !
   MPA_freqs=Xw%p
   !
   ! Poles and Residues
   !
   call timing('BSE(mpa)_poles',OPR='start')
   !
   MPred_rate=0._SP
   PPcond_Qrate=0._SP
   MP_Qerr=0._SP
   cond_numQ(:)=0._SP
   !
   MPA_R_par(1)%blc=0.0_SP
   MPA_E_par(1)%blc=0.0_SP
   !
   ! Thread-private MPA_Xo: PRIVATE allocatables start unallocated on each
   ! thread, so each thread allocates its own copy once with plain allocate
   ! (avoids YAMBO counter races and the unallocated-PRIVATE bug).
   !
   !$omp parallel default(shared),private(i_g1,i_g2,i_g1_mat,i_g2_mat, &
   !$omp & MPred,PPcond_rate,MP_err,cond_num,MPA_Xo,l_rim_w_g), &
   !$omp &        reduction(+:MPred_rate,PPcond_Qrate,MP_Qerr,cond_numQ)
   allocate(MPA_Xo(Xw%n_freqs))
   !$omp do collapse(2)
   do i_g2= ig_first, ig_last
     i_g2_mat=min(i_g2,X%ng)
     do i_g1=1, BS_n_g_W
       i_g1_mat=min(i_g1,X%ng)
       !
        l_RIM_W_g=(l_rim_w.and.i_g1<=RIM_W_ng.and.i_g2<=RIM_W_ng)
        if (RIM_W_is_diagonal) l_RIM_W_g=(l_RIM_W_g.and.(i_g1==i_g2))
        !
       if (l_RIM_W_g) then
         MPA_Xo(1:Xw%n_freqs)=RIM_W(1:Xw%n_freqs,iq_W,i_g1_mat,i_g2_mat)/2._SP
       else
         ! X_par(1): single-q buffer (io_RESPONSE). gamp on host (d2h above).
         MPA_Xo(1:Xw%n_freqs)=X_par(1)%blc(i_g1_mat,i_g2_mat,1:Xw%n_freqs)*isc%gamp(i_g1_mat,i_g2_mat)
       endif
       !
       ! mpa_RE_solver leaves MPred/cond_num untouched for np==1: init here.
       !
       MPred=.false.
       PPcond_rate=0._SP
       MP_err=0._SP
       cond_num(:)=0._SP
       !
       call mpa_RE_solver(X%mpa_npoles,MPA_freqs,MPA_Xo, &
&                         MPA_R_par(1)%blc(i_g1_mat,i_g2_mat,:),MPA_E_par(1)%blc(i_g1_mat,i_g2_mat,:), &
&                         Xw%mpa_solver,MPred,PPcond_rate,MP_err,cond_num)
       !
       !
       if (io_MULTIPOLE) then
         if (MPred) MPred_rate=MPred_rate+1._SP
         PPcond_Qrate=PPcond_Qrate+PPcond_rate
         MP_Qerr=MP_Qerr+MP_err
         cond_numQ(:)=cond_numQ(:)+cond_num(:)
       endif
       !
     enddo
   enddo
   !$omp end do
   deallocate(MPA_Xo)
   !$omp end parallel
   !
   ! Collect distributed G-columns before write (QP_mpa.F:422-425)
   !
    if (l_para_G_MPA) then
      call PP_redux_wait(MPA_R_par(1)%blc,COMM=PAR_COM_X_WORLD%COMM)
      call PP_redux_wait(MPA_E_par(1)%blc,COMM=PAR_COM_X_WORLD%COMM)
      if (io_MULTIPOLE) then
        ! MB: this reduction is needed for the I/O MPA part
        call PP_redux_wait(MPred_rate,COMM=PAR_COM_X_WORLD%COMM)
        call PP_redux_wait(PPcond_Qrate,COMM=PAR_COM_X_WORLD%COMM)
        call PP_redux_wait(MP_Qerr,COMM=PAR_COM_X_WORLD%COMM)
        call PP_redux_wait(cond_numQ,COMM=PAR_COM_X_WORLD%COMM)
      endif
    endif
   !
   call timing('BSE(mpa)_poles',OPR='stop')
   call timing('BSE(mpa)_poles_io',OPR='start')
   !
   ! I/O MPA: append poles for this iq (io_MPA reads MPA_E/R_par(1))
   !
   if (io_MULTIPOLE) then
     !
      MPred_rate=X%mpa_npoles-MPred_rate/real(BS_n_g_W,SP)**2
      PPcond_Qrate=PPcond_Qrate/real(BS_n_g_W,SP)**2
      MP_Qerr=MP_Qerr/real(BS_n_g_W,SP)**2
      cond_numQ(:)=cond_numQ(:)/real(BS_n_g_W,SP)**2
     !
     call io_control(ACTION=OP_APP_CL,COM=REP,SEC=(/2*iq_W,2*iq_W+1/),ID=ID2,&
&                    COMM=PAR_COM_WORLD,DO_IT=.TRUE.)
     io_err=io_MPA(X,Xw,ID2)
     call msg('r',' :: [MPA@Q'//trim(intc(iq_W))//'] PP cond fix/tot   ',PPcond_Qrate)
     call msg('r',' :: [MPA@Q'//trim(intc(iq_W))//'] Mean np reduction ',MPred_rate)
     call msg('r',' :: [MPA@Q'//trim(intc(iq_W))//'] Mean Xm rel dev   ',MP_Qerr)
     if (Xw%mpa_solver=='LA') call msg('r',' :: [MPA@Q'//trim(intc(iq_W))//'] Condition numbers ',cond_numQ)
     !
   endif
   !
   ! Optional retention: copy working buffer into per-q storage
   !
    if (l_keep_MPA_in_memory) then
      call MATRIX_copy(MPA_E_par(1),MPA_E_store(iq_W),.TRUE.)
      MPA_E_store(iq_W)%blc = MPA_E_par(1)%blc

      call MATRIX_copy(MPA_R_par(1),MPA_R_store(iq_W),.TRUE.)
      MPA_R_store(iq_W)%blc = MPA_R_par(1)%blc
    endif
    !
    call timing('BSE(mpa)_poles_io',OPR='stop')
   !
 enddo
 !
 if (io_MULTIPOLE) call IO_and_Messaging_switch("RESTORE")
 !
 call elemental_collision_free(isc)
 !
 call MPA_free()
 !
 call MATRIX_reset(X_par(1))
 deallocate(X_par)
 !
 ! MPA workspace cleanup (working buffer only; store(:) survives if retained)
 !
 call MATRIX_reset(MPA_E_par(1))
 call MATRIX_reset(MPA_R_par(1))
 if (allocated(MPA_E_par)) deallocate(MPA_E_par)
 if (allocated(MPA_R_par)) deallocate(MPA_R_par)
 if (l_para_G_MPA) call PP_indexes_reset(px)
 !
end subroutine
