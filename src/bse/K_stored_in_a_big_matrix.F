!
! License-Identifier: GPL
!
! Copyright (C) 2018 The Yambo Team
!
! Authors (see AUTHORS file for details): AM MG DS
!
subroutine K_stored_in_a_big_matrix(i_BS_mat,BS_mat_dim,BS_mat,iq,what)
 !
 ! Here I fill the kernel (coupling included) in an entire BIG matrix to be
 ! used for diagonalization and/or invertion
 !
 use pars,           ONLY:SP,cI,cZERO,cONE
 use parallel_int,   ONLY:PP_redux_wait
 use BS,             ONLY:BS_K_coupling,n_BS_blks,BS_res_ares_n_mat,&
&                         l_BS_ares_from_res,BS_dyn_kernel_full,BS_dyn_kernel_hermitianize,BS_dump_full_H
 use BS_solvers,     ONLY:BS_blk,BS_H_dim,BS_K_dim,&
&                         BSS_perturbative_width,run_inversion,run_Haydock,&
&                         BSS_eh_E,BSS_eh_W
 use interfaces,      ONLY:K_dump_full_matrix
 use timing_m,       ONLY:timing
 use com,           ONLY:msg
 !
#include <y_memory.h>
 !
 integer,      intent(in) :: i_BS_mat,BS_mat_dim,iq
 complex(SP), allocatable, target :: BS_mat(:,:) 
 character(*), intent(in) :: what
 !
 ! Work Space
 !
 complex(SP) :: Mij
 integer     :: i_Tk,i_Tp,i_r,i_c,H_pos(2),H_shift(2),i_B,i,j
 !
 if (allocated(BS_mat)) return
 !
 call timing('BSE matrix filling',OPR='start')
 !
 YAMBO_ALLOC(BS_mat,(BS_mat_dim,BS_mat_dim))
 BS_mat    = cZERO
 !
 do i_B=1,n_BS_blks
   !
   i_Tk  =BS_blk(i_B)%iT_k
   i_Tp  =BS_blk(i_B)%iT_p
   !
   if (i_BS_mat/=BS_blk(i_B)%ira_k .and. BS_res_ares_n_mat==2) cycle
   !
   H_shift=0
   if(BS_blk(i_B)%mode=="C") H_shift(2)=BS_K_dim(1)
   if(BS_blk(i_B)%mode=="A") H_shift(:)=BS_K_dim(1)
   !
   do i_r=1,BS_blk(i_B)%size(1)
     !
     H_pos(1)=BS_blk(i_B)%coordinate(1)+i_r-1
     !
     do i_c=1,BS_blk(i_B)%size(2)
       !
       H_pos(2)=BS_blk(i_B)%coordinate(2)+i_c-1
       !
       ! Static BSE storage contains only the upper triangle and rebuilds
       ! the other half from known symmetries.  Dynamic MPA stores every
       ! ordered element, so no element may be skipped here.
       if (.not.BS_dyn_kernel_full) then
         if (H_pos(1)+H_shift(1)>H_pos(2)+H_shift(2)) cycle
         if (l_BS_ares_from_res.and.H_pos(1)>H_pos(2)) cycle
       endif
       !
       if(     allocated(BS_blk(i_B)%mat)) Mij=BS_blk(i_B)%mat(i_r,i_c)
       if(.not.allocated(BS_blk(i_B)%mat)) Mij=cZERO
       !
       ! Add energies to the diagonal
       !
       if(H_pos(1)+H_shift(1)==H_pos(2)+H_shift(2)) then
         if (BS_dyn_kernel_full) then
           ! Preserve the raw ordered diagonal until the global projection.
           Mij=Mij+BSS_eh_E(H_pos(1)+H_shift(1))*cONE
         else
           Mij=real(Mij)+BSS_eh_E(H_pos(1)+H_shift(1))*cONE
         endif
         if (allocated(BSS_eh_W).and..not.BSS_perturbative_width) then
           Mij=Mij+cI*BSS_eh_W(H_pos(1)+H_shift(1))
         endif
       endif
       !
       if (BS_dyn_kernel_full) then
         ! Dynamic MPA currently uses the resonant block.  Store its ordered
         ! value directly; do not manufacture the opposite element here.
         BS_mat(H_pos(1)+H_shift(1),H_pos(2)+H_shift(2))=Mij
       else
       select case(BS_blk(i_B)%mode)
       case("R")
           BS_mat(H_pos(1)            ,H_pos(2)            )=       Mij
           ! The resonant block is hermitial
           BS_mat(H_pos(2)            ,H_pos(1)            )= conjg(Mij)
           if (l_BS_ares_from_res.and.BS_K_coupling) then
             ! The anti-resonant block is A=-R*
             BS_mat(H_pos(1)+BS_K_dim(1),H_pos(2)+BS_K_dim(1))=-conjg(Mij)
             ! The anti-resonant block is hermitian
             BS_mat(H_pos(2)+BS_K_dim(1),H_pos(1)+BS_K_dim(1))=      -Mij
           endif
       case("C")
           BS_mat(H_pos(1)            ,H_pos(2)+BS_K_dim(1))=       Mij
           ! Anti-coupling from coupling: the whole BSE matrix is Pseudo-HErmitian
           BS_mat(H_pos(2)+BS_K_dim(1),H_pos(1)            )=-conjg(Mij)
           if (l_BS_ares_from_res) then
             ! The coupling block and the anti-coupling block are symmetric
             BS_mat(H_pos(2)            ,H_pos(1)+BS_K_dim(1))=       Mij
             BS_mat(H_pos(1)+BS_K_dim(1),H_pos(2)            )=-conjg(Mij)
           endif
       case("A")
         ! The anti-resonant block is hermitial
         if(BS_res_ares_n_mat==1) then
           BS_mat(H_pos(1)+BS_K_dim(1),H_pos(2)+BS_K_dim(1))=       Mij
           BS_mat(H_pos(2)+BS_K_dim(1),H_pos(1)+BS_K_dim(1))= conjg(Mij)
         else
           BS_mat(H_pos(1)            ,H_pos(2)            )=       Mij
           BS_mat(H_pos(2)            ,H_pos(1)            )= conjg(Mij)
         endif
       case("Q")
           BS_mat(H_pos(1)+BS_K_dim(1),H_pos(2)            )=       Mij
           ! Coupling from anti-coupling: the whole BSE matrix is Pseudo-HErmitian
           BS_mat(H_pos(2)            ,H_pos(1)+BS_K_dim(1))=-conjg(Mij)
       end select
       endif
       !
     enddo
     !
   enddo
   !
 enddo
 !
 call timing('BSE matrix filling',OPR='stop')
 !
 ! ALL 2 ALL
 !===========
 !
 ! BS_mat lives in the entire world.
 !
 call timing('BSE matrix (REDUX)',OPR='start')
 call PP_redux_wait(BS_mat)
 call timing('BSE matrix (REDUX)',OPR='stop')
 !
 if (BS_dyn_kernel_full.and.BS_dump_full_H) then
   ! The first dump is deliberately placed before the projection: it keeps
   ! the independently computed ordered matrix available for external
   ! diagonalization and for comparing MPI/rank-dependent kernel elements.
   call msg('s','[BSE] Dynamic MPA: Dumping the full matrix') ! Miki Debug
   call K_dump_full_matrix(iq,BS_mat,.FALSE.)
 endif
 !
 if (BS_dyn_kernel_full.and.BS_dyn_kernel_hermitianize) then
   ! Hermitianize only after every ordered block has reached the global
   ! dense matrix.  Doing this earlier would mix MPI-local partial data.
   do i=1,BS_mat_dim
     BS_mat(i,i)=real(BS_mat(i,i),kind=SP)
     do j=i+1,BS_mat_dim
       Mij=0.5_SP*(BS_mat(i,j)+conjg(BS_mat(j,i)))
       BS_mat(i,j)=Mij
       BS_mat(j,i)=conjg(Mij)
     enddo
  enddo
 endif
 !
 if (BS_dyn_kernel_full.and.BS_dump_full_H.and.BS_dyn_kernel_hermitianize) then
   ! Store the projected object in a separate database.  Keeping two files
   ! makes it impossible to confuse the raw ordered matrix with its effective
   ! Hermitian approximation during post-processing.
   call msg('s','[BSE] Dynamic MPA: Dumping the Hermitianized matrix') ! Miki Debug
   call K_dump_full_matrix(iq,BS_mat,.TRUE.)
 endif
 !
 ! After the folding the %mat's are (in some cases) not needed anymore 
 !
 if (run_Haydock) return
 if (run_inversion.and.what=="DIAGO") return
 !
 call BS_blks_free(i_BS_mat)
 !
end subroutine K_stored_in_a_big_matrix

subroutine K_dump_full_matrix(iq,BS_mat,hermitianized)
 ! Write one globally reduced dense BSE matrix and the transition metadata
 ! needed to interpret its rows and columns outside Yambo.
 use pars,        ONLY:SP,schlen,IP_YIO,LP_YIO
 use stderr,      ONLY:intc
 use BS,          ONLY:BS_dyn_omega,BS_dyn_eta,BS_dyn_kernel_full,BSqpts
 use BS_solvers,  ONLY:BSS_eh_table,BSS_eh_E,BSS_eh_f,BSS_eh_f_RES,BSS_eh_W
 use IO_m,        ONLY:io_status,REP,OP_WR_CL
 use IO_int,      ONLY:io_control,io_connect,io_header,io_disconnect,&
&                       def_variable_elemental,io_variable_elemental,&
&                       def_variable_bulk,io_variable_bulk
 use timing_m,    ONLY:timing
 implicit none
 integer, intent(in) :: iq
 complex(SP), target :: BS_mat(:,:)
 logical, intent(in) :: hermitianized
 integer :: ID,io_err,N,var_id
 character(schlen) :: db_name
 logical :: marker
 logical :: connected
 !
 N=size(BS_mat,1)
 marker=hermitianized
 if (hermitianized) then
   db_name='BS_full_H_hermitian_Q'//trim(intc(iq))
 else
   db_name='BS_full_H_raw_Q'//trim(intc(iq))
 endif
 !
 call timing('io_'//db_name,OPR='start')
 !
 ! This routine creates a new dump database.  OP_RD_CL would make
 ! io_connect look for an existing database and consequently fail for
 ! the first dump (and for every new q-point name).
 call io_control(ACTION=OP_WR_CL,COM=REP,SEC=(/1,2/),ID=ID)
 connected=.false.
 io_err=io_connect(desc=trim(db_name),type=2,ID=ID)
 if (io_err/=0) goto 1
 connected=.true.
 io_err=io_header(ID,QPTS=.true.,R_LATT=.true.,WF=.true.,IMPOSE_SN=.true.,&
&                XC_KIND="K_WF force Xs",CUTOFF=.true.,&
&                PARALLEL_STRUCTURE="Response_T_space",FRAG=.false.)
 if (io_err/=0) goto 1
 ! Scalar metadata describes the matrix and the approximation applied to it.
 call def_variable_elemental(ID,'BSE_FULL_H_DIM',1,IP_YIO,0)
 call io_variable_elemental(ID,VAR='Dense BSE matrix dimension',I0=N)
 call def_variable_elemental(ID,'BSE_FULL_H_HERMITIANIZED',1,LP_YIO,0)
 call io_variable_elemental(ID,VAR='Dense BSE matrix hermitianized',L0=marker)
 call def_variable_elemental(ID,'BSE_FULL_H_DYNAMIC',1,LP_YIO,0)
 call io_variable_elemental(ID,VAR='Dense BSE matrix dynamic',L0=BS_dyn_kernel_full)
 call def_variable_elemental(ID,'BSE_FULL_H_Q_INDEX',1,IP_YIO,0)
 call io_variable_elemental(ID,VAR='Dense BSE matrix q index',I0=iq)
 if (allocated(BSqpts)) then
   if (size(BSqpts,2)>=iq) then
     call def_variable_elemental(ID,'BSE_FULL_H_Q_POINT',3,SP,0)
     call io_variable_elemental(ID,VAR='Dense BSE matrix q point',R1=BSqpts(:,iq))
   endif
 endif
 call def_variable_elemental(ID,'BSE_FULL_H_OMEGA',1,SP,0)
 call io_variable_elemental(ID,VAR='Dense BSE matrix omega',R0=BS_dyn_omega)
 call def_variable_elemental(ID,'BSE_FULL_H_ETA',1,SP,0)
 call io_variable_elemental(ID,VAR='Dense BSE matrix eta',R0=BS_dyn_eta)
 call def_variable_elemental(ID,'',0,0,1)
 ! Dense matrix and transition-basis metadata.
 var_id=1
 call def_variable_bulk(ID,'BSE_FULL_H',var_id,(/2,N,N/),SP)
 call io_variable_bulk(ID,var_id,C2=BS_mat)
 if (allocated(BSS_eh_table)) then
   var_id=var_id+1
   call def_variable_bulk(ID,'BS_TABLE',var_id,(/size(BSS_eh_table,1),size(BSS_eh_table,2)/),IP_YIO)
   call io_variable_bulk(ID,var_id,I2=BSS_eh_table)
 endif
 if (allocated(BSS_eh_E)) then
   var_id=var_id+1
   call def_variable_bulk(ID,'BS_EH_ENERGIES',var_id,(/size(BSS_eh_E)/),SP)
   call io_variable_bulk(ID,var_id,R1=BSS_eh_E)
 endif
 if (allocated(BSS_eh_f)) then
   var_id=var_id+1
   call def_variable_bulk(ID,'BS_EH_OSCILLATORS',var_id,(/size(BSS_eh_f)/),SP)
   call io_variable_bulk(ID,var_id,R1=BSS_eh_f)
 endif
 if (allocated(BSS_eh_f_RES)) then
   var_id=var_id+1
   call def_variable_bulk(ID,'BS_EH_OSCILLATORS_RES',var_id,(/size(BSS_eh_f_RES)/),SP)
   call io_variable_bulk(ID,var_id,R1=BSS_eh_f_RES)
 endif
 if (allocated(BSS_eh_W)) then
   var_id=var_id+1
   call def_variable_bulk(ID,'BS_EH_WIDTHS',var_id,(/size(BSS_eh_W)/),SP)
   call io_variable_bulk(ID,var_id,R1=BSS_eh_W)
 endif
 io_err=io_status(ID)
1 if (connected) call io_disconnect(ID)
 !
 call timing('io_'//db_name,OPR='stop')
end subroutine K_dump_full_matrix
