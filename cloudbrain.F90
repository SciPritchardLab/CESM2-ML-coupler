#define BRAINDEBUG
module cloudbrain
use constituents,    only: pcnst
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use cam_history,         only: outfld, addfld, add_default
use physconst,       only: gravit,cpair,latvap,latice
use spmd_utils, only: masterproc,iam
use camsrfexch,       only: cam_out_t, cam_in_t
use constituents,    only: cnst_get_ind
use physics_types,    only: physics_state,  physics_ptend, physics_ptend_init
use physics_buffer,     only: physics_buffer_desc
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
  integer, parameter :: outputlength_atm = 104 ! 26*4
  integer, parameter :: outputlength_lnd = 8 !  8 scalars

  type(network_type) :: cloudbrain_net_atm
  type(network_type) :: cloudbrain_net_lnd

  real :: inp_sub(inputlength)
  real :: inp_div(inputlength)
  real :: out_weight_atm(outputlength_atm)
  real :: out_weight_lnd(outputlength_lnd)
 
  public neural_net, init_neural_net
  contains

  subroutine neural_net (state,pbuf,nn_solin,cam_in,ztodt,ptend,cam_out)
 ! note state is meant to have the "BP" state saved earlier. 

   implicit none

   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in)               :: nn_solin(pcols) 
   type(cam_in_t),intent(in)          :: cam_in
   real(r8), intent(in) :: ztodt 
   type(physics_ptend),intent(out) :: ptend            ! indivdual parameterization tendencies
   type(cam_out_t),     intent(inout) :: cam_out

    ! local variables
   real :: input(pcols,inputlength)
   real :: output_atm(pcols,outputlength_atm)
   real :: output_lnd(pcols,outputlength_lnd)
   integer :: i,k,ncol,ixcldice,ixcldliq,ii,kk
   real (r8) :: s_bctend(pcols,pver), q_bctend(pcols,pver), qc_bctend(pcols,pver), qi_bctend(pcols,pver), qafter, safter
   logical :: doconstraints
   logical ::  lq(pcnst)
   ncol  = state%ncol
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   lq(:)        = .FALSE.
   lq(1) = .TRUE.
   lq(ixcldliq) = .TRUE.
   lq(ixcldice) = .TRUE.
   call physics_ptend_init(ptend, state%psetcols, 'neural-net', ls=.true.,lq=lq)   ! Initialize local physics_ptend object

   doconstraints = .true.
   
   s_bctend(:,:) = 0.
   q_bctend(:,:) = 0.
   qc_bctend(:,:) = 0.
   qi_bctend(:,:) = 0.
 
  
  ! Input variable order hardwired here:
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
      output_atm(i,:) = cloudbrain_net_atm % output(input(i,:))
    end do

#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG output_atm = ',output_atm(1,:)
      endif
#endif
   ! output_atm normalization (un-weighting, really).
   do i=1,ncol
     do k=1,outputlength_atm
      output_atm(i,k) = output_atm(i,k) / out_weight_atm(k)
     end do
   end do

#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG out post scale = ',output_atm(1,:)
      endif
#endif


! ---------- 1. NN output_atm to atmosphere forcing --------

   q_bctend(:ncol,:pver) = real(output_atm(:ncol,1:pver),r8) ! kg/kg/s 
   s_bctend(:ncol,1:pver) = cpair*real(output_atm(:ncol,(pver+1):(2*pver)),r8) ! K/s --> J/kg/s (ptend expects that)
   qc_bctend(:ncol,:pver) = real(output_atm(:ncol,(2*pver+1):(3*pver)),r8) ! kg/kg/s 
   qi_bctend(:ncol,:pver) = real(output_atm(:ncol,(3*pver+1):(4*pver)),r8) ! kg/kg/s 

! -- atmos positivity constraints ---- 
   if (doconstraints) then
   do i=1,ncol
     do k=1,pver
! energy positivity:
       safter = state%s(i,k) + s_bctend(i,k)*ztodt ! predicted DSE after NN tendency
       if (safter .lt. 0.) then ! can only happen when bctend < 0...
         s_bctend(i,k) = s_bctend(i,k) + abs(safter)/ztodt ! in which case reduce cooling rate
         write (iulog,*) 'HEY CBRAIN made a negative absolute temperature, corrected but BEWARE!!!'
       endif

 ! vapor positivity:
       qafter = state%q(i,k,1) + q_bctend(i,k)*ztodt ! predicted vapor after NN tendency
       if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
         q_bctend(i,k) = q_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
       endif
 ! liquid positivity:
       qafter = state%q(i,k,ixcldliq) + qc_bctend(i,k)*ztodt ! predicted liquid after NN tendency
       if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
         qc_bctend(i,k) = qc_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
       endif
! ice positivity:
       qafter = state%q(i,k,ixcldice) + qi_bctend(i,k)*ztodt ! predicted ice after NN tendency
       if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
         qi_bctend(i,k) = qi_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
       endif
     end do
   end do
   endif
! Wire to ptend:
    ptend%s(:ncol,:pver) = s_bctend(:ncol,:pver)
    ptend%q(:ncol,:pver,1) = q_bctend(:ncol,:pver)
    ptend%q(:ncol,:pver,ixcldliq) = qc_bctend(:ncol,:pver)
    ptend%q(:ncol,:pver,ixcldice) = qi_bctend(:ncol,:pver)

! ------------- 2. NN output to land forcing ---------

    do i=1,ncol
      output_lnd(i,:) = cloudbrain_net_lnd % output(input(i,:))
    end do

#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG output_lnd = ',output_lnd(1,:)
      endif
#endif
   ! output_lnd normalization (un-weighting, really).
   do i=1,ncol
     do k=1,outputlength_lnd
      output_lnd(i,k) = output_lnd(i,k) / out_weight_lnd(k)
     end do
   end do

#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG output_lnd post scale = ',output_lnd(1,:)
      endif
#endif

! populate the all-important cam_out structure that packages up info needed by
! CLM:

   ! Phase 1: Auto-repopulate most properties based on the lowest-model-level state given that it has been updated by NN tendencies upstream during physics_update call::
   call cam_export (state,cam_out,pbuf) 
   ! Phase 2: overwrite what remains to be overwritten with the 8 items learned,
   ! mimicking Liran's NN2L output variables:

   ! assuming the following order.
   ! 'NN2L_FLWDS', 'NN2L_PRECC', 'NN2L_PRECSC', 'NN2L_SOLL', 'NN2L_SOLLD',
   ! 'NN2L_SOLS', 'NN2L_SOLSD', 'NN2L_NETSW']
  cam_out%flwds(:ncol) = output_lnd(:ncol,1)
  cam_out%precc(:ncol) = output_lnd(:ncol,2)
  cam_out%precsc(:ncol) = output_lnd(:ncol,3)
  cam_out%cam_out%soll(:ncol) = output_lnd(:ncol,4)
  cam_out%cam_out%solld(:ncol) = output_lnd(:ncol,5)
  cam_out%cam_out%sols(:ncol) = output_lnd(:ncol,6)
  cam_out%cam_out%solsd(:ncol) = output_lnd(:ncol,7)
  cam_out%netsw(:ncol) = output_lnd(:ncol,8)
  cam_out%precl(:ncol) = 0. ! probably redundant (when SP diagnosed)... won't be otherwise 

! TODO: Is it important to wire things into the PBUF arteries (need to ask
! Walter...)

! Pasting Slack notes from Jul 16 when we thought this through a bit more:
! Notes:
!srfrad = fsns+flwds
!netsw = srfrad-flwds
!So netsw is a shortcut to get the downwelling shortwave flux felt by the land
!model
!LHFLX and SHFLX are saved from the previous timestep (shf and lhf in the code)
!Check for SOLIN , whose variable is solin in the radiation.F90, stored in
!solin_m in the local sw routine.
! ***  NEED to calculate psl upstream of cam_export somewhere?
!psl is calculated using cpslec in diag_phys_writeout so it will not be emulated
!because it is a direct FORTRAN diagnostic calculation
!Assumptions made on the way:
!(1) Momentum tendency from physics phase 1 before coupling = 0
!TODO for @Tom Beucler :
!Has to be done before training
!Check whether the solar radiation sub-components sum up to the net shortwave
! ****
!Add consistent diagnostic calculation of PSL so that we don't have to predict it
!using the NN
! **** (think this is done via cam_export call above?)
!Add diagnostic calculation for QBOT, TBOT, THBOT , UBOT , VBOT
!cam_out%thbot(i) = state%t(i,pver) * state%exner(i,pver) 
! *** NEED TO LOOK INTO:
!4. Figure out how geopotential height and pressure are repopulated after
!heating/moistening tendencies are applied. This would give us unambiguous values
!for ZBOT and PBOT
!5. Same for RHO
!cam_out%rho(i)  = cam_out%pbot(i)/(rair*cam_out%tbot(i)) 
!6. Same for SRFRAD=NETSW+FLWDS
end subroutine neural_net

  subroutine init_neural_net()

    implicit none

    call cloudbrain_net_atm % load('/scratch/07064/tg863631/fortran_models/BF_RG_config.txt')
    write (iulog,*) '------- FKB: loaded atm network from txt file -------'

    call cloudbrain_net_lnd % load('INSERT/PATH/TO/LND/NN.txt')
    write (iulog,*) '------- FKB: loaded lnd network from txt file -------'
    
    open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/inp_sub.txt',status='old',action='read')
    read(555,*) inp_sub(:)
    
    open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/inp_div.txt',status='old',action='read')
    read(555,*) inp_div(:)
    
    open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/scale_dict_output.txt',status='old',action='read')
    read(555,*) out_weight_atm(:)

    open (unit=555,file='/INSERT/PATH/TO/LND/OUTPUT_scale.txt',status='old',action='read')
    read(555,*) out_weight_lnd(:)
#ifdef BRAINDEBUG
    if (masterproc) then
       write (iulog,*) 'BRAINDEBUG read input norm inp_sub=', inp_sub(:)
       write (iulog,*) 'BRAINDEBUG read input norm inp_div=', inp_div(:)       
       write (iulog,*) 'BRAINDEBUG read output norm out_weight_atm=', out_weight_atm(:)       
       write (iulog,*) 'BRAINDEBUG read output norm out_weight_lnd=', out_weight_lnd(:)       
    endif
#endif

  end subroutine init_neural_net
    
end module cloudbrain
