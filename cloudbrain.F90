#define BRAINDEBUG
module cloudbrain
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use cam_history,         only: outfld, addfld, add_default
use physconst,       only: gravit,cpair,latvap,latice
use spmd_utils, only: masterproc
use camsrfexch,       only: cam_out_t, cam_in_t
use constituents,    only: cnst_get_ind
use physics_types,    only: physics_state,  physics_ptend
use cam_logfile,       only: iulog

! -------- NEURAL-FORTRAN --------
! imports
use mod_kinds, only: ik, rk
use mod_network , only: network_type
use mod_ensemble, only: ensemble_type
! --------------------------------

  implicit none
  save 

  private
  ! Define variables for this entire module
  integer, parameter :: inputlength = 108 ! 26*4 + 4 scalars
  integer, parameter :: outputlength = 112 ! 26*4 + 8 scalars

  type(network_type) :: cloudbrain_net

  real :: inp_sub(inputlength)
  real :: inp_div(inputlength)
  real :: out_weight(outputlength)

  public neural_net, init_neural_net
  contains

  subroutine neural_net (state,nn_solin,cam_in,ztodt,ptend,cam_out)
 ! note state is meant to have the "BP" state saved earlier. 
   type(physics_state), intent(in)    :: state
    real(r8), intent(in)               :: nn_solin(pcols) 
    type(cam_in_t),intent(in)          :: cam_in
    real(r8), intent(in) :: ztodt 
    type(physics_ptend),intent(inout)     :: ptend            ! indivdual parameterization tendencies
    type(cam_out_t),     intent(inout) :: cam_out

    ! local variables
    real :: input(pcols,inputlength)
    real :: output(pcols,outputlength)
    integer :: i,k,ncol,ixcldice,ixcldliq
    real (r8) :: s_bctend(pcols,pver), q_bctend(pcols,pver), qc_bctend(pcols,pver), qi_bctend(pcols,pver), qafter, safter
    
    ncol  = state%ncol
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
  
  ! Ankitesh says on Slack that ['QBP','TBP','CLDLIQBP','CLDICEBP','PS', 'SOLIN', 'SHFLX', 'LHFLX']
    input(:ncol,1:pver) = state%q(:ncol,:pver,1)
    input(:ncol,(pver+1):(2*pver)) = state%t(:ncol,:pver)
    input(:ncol,(2*pver+1):(3*pver)) = state%q(:ncol,:pver,ixcldliq)
    input(:ncol,(3*pver+1):(4*pver)) = state%q(:ncol,:pver,ixcldice)
    input(:ncol,(4*pver+1)) = state%ps(:ncol)
    input(:ncol,(4*pver+2)) = nn_solin(:ncol) ! WARNING this is being lazily mined from part of SP solution... should be avoidable in future when bypassing SP totally but will take work.
    input(:ncol,(4*pver+3)) = cam_in%shf(:ncol)
    input(:ncol,(4*pver+4)) = cam_in%lhf(:ncol) 
 

#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG input pre norm=',input(1,:)
      endif
#endif

    ! 2. Normalize input
    do k=1,inputlength
      input(:ncol,k) = (input(:ncol,k) - inp_sub(k))/inp_div(k)
    end do
#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG input post norm=',input(1,:)
      endif
#endif

    do i=1,ncol
      output(i,:) = cloudbrain_net % output(input(i,:))
    end do

#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG output = ',output(1,:)
      endif
#endif
   ! output normalization (un-weighting, really).
   do i=1,ncol
     do k=1,outputlength
      output(i,k) = output(i,k) / out_weight(k)
     end do
   end do

#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG out post scale = ',output(1,:)
      endif
#endif

! From Tom's Slack the first 120 elements of the output vector contain atmos forcing, ordered (but need to check with Ankitesh and vs. output norms) as:
! TBCTEND=(TBC-TBP)/DT,QBCTEND=(QBC-QBP)/DT,CLDLIQBCTEND=(CLDLIQBC-CLDLIQBP)/DT,CLDICEBCTEND=(CLDICEBC-CLDICEBP)/DT,
! Then the land variables (there should be 8) are as:
! NN2L_FLWDS,NN2L_NETSW(misleading name as it's downwelling shortwave),NN2L_PRECC,NN2L_PRECSC,NN2L_SOLL,NN2L_SOLLD,NN2L_SOLS,NN2L_SOLSD)
! Thus total NN output vector length is:
! 4*pver + 8 = 4*26 +8 = 112

! ---------- 1. NN output to atmosphere forcing --------

   s_bctend = cpair*output(:ncol,1:pver) ! K/s --> J/kg/s
   q_bctend = output(:ncol,(pver+1):(2*pver)) ! kg/kg/s 
   qc_bctend = output(:ncol,(2*pver+1):(3*pver)) ! kg/kg/s 
   qi_bctend = output(:ncol,(3*pver+1):(4*pver)) ! kg/kg/s 

! -- atmos positivity constraints ---- 
   do i=1,ncol
     do k=1,pver
! energy positivity:
       safter = state%s(i,k) + s_bctend(i,k)*ztodt ! predicted DSE after NN tendency
       if (safter < 0.) then ! can only happen when bctend < 0...
         s_bctend(i,k) = s_bctend(i,k) + abs(safter)/ztodt ! in which case reduce cooling rate
         write (iulog,*) 'HEY CBRAIN made a negative absolute temperature, corrected but BEWARE!!!'
       endif

 ! vapor positivity:
       qafter = state%q(i,k,1) + q_bctend(i,k)*ztodt ! predicted vapor after NN tendency
       if (qafter < 0.) then ! can only happen when qbctend < 0...
         q_bctend(i,k) = q_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
       endif
 ! liquid positivity:
       qafter = state%q(i,k,ixcldliq) + qc_bctend(i,k)*ztodt ! predicted liquid after NN tendency
       if (qafter < 0.) then ! can only happen when qbctend < 0...
         qc_bctend(i,k) = qc_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
       endif
! ice positivity:
       qafter = state%q(i,k,ixcldice) + qi_bctend(i,k)*ztodt ! predicted ice after NN tendency
       if (qafter < 0.) then ! can only happen when qbctend < 0...
         qi_bctend(i,k) = qi_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
       endif
     end do
   end do

! Wire to ptend:
   ptend%ls = .true.
   ptend%s(:ncol,:pver) = s_bctend(:ncol,:pver) ! physics_tend expects J/kg/s for heating 
   ptend%lq = .true.
   ptend%q(:ncol,:pver,1) = q_bctend(:ncol,:pver)
   ptend%q(:ncol,:pver,ixcldliq) = qc_bctend(:ncol,:pver)
   ptend%q(:ncol,:pver,ixcldice) = qi_bctend(:ncol,:pver)
   ptend%lu = .false.

! ------------- 2. NN output to land forcing ---------

end subroutine neural_net

  subroutine init_neural_net()

    call cloudbrain_net % load('/scratch/07064/tg863631/fortran_models/BF_RG_config.txt')
    write (iulog,*) '------- FKB: loaded network from txt file -------'
    
    open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/inp_sub.txt',status='old',action='read')
    read(555,*) inp_sub(:)
    
    open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/inp_div.txt',status='old',action='read')
    read(555,*) inp_div(:)
    
    open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/scale_dict_output.txt',status='old',action='read')
    read(555,*) out_weight(:)
 #ifdef BRAINDEBUG
    if (masterproc) then
       write (iulog,*) 'BRAINDEBUG read input norm inp_sub=', inp_sub(:)
       write (iulog,*) 'BRAINDEBUG read input norm inp_div=', inp_div(:)       
       write (iulog,*) 'BRAINDEBUG read output norm out_weight=', out_weight(:)       
    endif
 #endif

  end subroutine init_neural_net
    
end module cloudbrain
