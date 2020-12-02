#define CBRAIN
#ifdef CBRAIN
#define BRAINDEBUG
#define RHDEBUG
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
!  integer, parameter :: outputlength = 112 ! 26*4 + 8 scalars
  integer, parameter :: outputlength = 111 ! 26*4 + 7 scalars (error, Ankitesh forgot one of the NN2L outputs)
  logical, parameter :: input_rh = .true. ! toggle to switch from q --> RH input

  type(network_type) :: cloudbrain_net

  real :: inp_sub(inputlength)
  real :: inp_div(inputlength)
  real :: out_weight(outputlength)

  public neural_net, init_neural_net
  contains

  subroutine neural_net (state,nn_solin,cam_in,ztodt,ptend,cam_out)
 ! note state is meant to have the "BP" state saved earlier. 

   implicit none

   type(physics_state), intent(in)    :: state
   real(r8), intent(in)               :: nn_solin(pcols) 
   type(cam_in_t),intent(in)          :: cam_in
   real(r8), intent(in) :: ztodt 
   type(physics_ptend),intent(out) :: ptend            ! indivdual parameterization tendencies
   type(cam_out_t),     intent(out) :: cam_out

    ! local variables
   real :: input(pcols,inputlength)
   real :: output(pcols,outputlength)
   integer :: i,k,ncol,ixcldice,ixcldliq,ii,kk
   real (r8) :: s_bctend(pcols,pver), q_bctend(pcols,pver), qc_bctend(pcols,pver), qi_bctend(pcols,pver), qafter, safter
   logical :: doconstraints
   logical ::  lq(pcnst)
   real :: rh_loc
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
 
  ! Ankitesh says on Slack that ['QBP','TBP','CLDLIQBP','CLDICEBP','PS', 'SOLIN', 'SHFLX', 'LHFLX']
    if (input_rh) then
       do i = 1,ncol
         do k=1,pver
           ! Port of tom's RH =  Rv*p*qv/(R*esat(T))
           rh_loc = 461.*state%pmid(i,k)*state%q(i,k,1)/(287.*tom_esat(real(state%t(i,k)))) ! note function tom_esat below refercing SAM's sat.F90
#ifdef RHDEBUG
           if (masterproc) then
             write (iulog,*) 'RHDEBUG:p,q,T,RH=',state%pmid(i,k),state%q(i,k,1),state%t(i,k),rh_loc
           endif
#endif
           input(i,k) = rh_loc
         end do
       end do
    else
      input(:ncol,1:pver) = state%q(:ncol,:pver,1) ! specific humidity input
    endif
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

! ['QBCTEND','TBCTEND','CLDLIQBCTEND', 'CLDICEBCTEND', 'NN2L_FLWDS',
! 'NN2L_PRECC', 'NN2L_PRECSC', 'NN2L_SOLL', 'NN2L_SOLLD', 'NN2L_SOLS',
! 'NN2L_SOLSD', 'NN2L_NETSW'] 

! Thus total NN output vector length is:
! 4*pver + 8 = 4*26 +8 = 112

! ---------- 1. NN output to atmosphere forcing --------

   q_bctend(:ncol,:pver) = real(output(:ncol,1:pver),r8) ! kg/kg/s 
   s_bctend(:ncol,1:pver) = cpair*real(output(:ncol,(pver+1):(2*pver)),r8) ! K/s --> J/kg/s (ptend expects that)
   qc_bctend(:ncol,:pver) = real(output(:ncol,(2*pver+1):(3*pver)),r8) ! kg/kg/s 
   qi_bctend(:ncol,:pver) = real(output(:ncol,(3*pver+1):(4*pver)),r8) ! kg/kg/s 

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
! TODO:

end subroutine neural_net

  subroutine init_neural_net()

    implicit none
    
    if (input_rh) then
      call cloudbrain_net % load('/scratch/07064/tg863631/fortran_models/RH_RGV1_config.txt')
    else
      call cloudbrain_net % load('/scratch/07064/tg863631/fortran_models/BF_RG_config.txt')
    end if
    write (iulog,*) '------- FKB: loaded network from txt file -------'
    
    if (input_rh) then
      open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/inp_sub_RH.txt',status='old',action='read')
    else
      open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/inp_sub.txt',status='old',action='read')
    end if
    read(555,*) inp_sub(:)
    
    if (input_rh) then
      open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/inp_div_RH.txt',status='old',action='read')
    else
      open (unit=555,file='/scratch/07064/tg863631/frontera_data/data/inp_div.txt',status='old',action='read')
    end if
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

  real function tom_esat(T) 
  ! For consistency with the python port of Tom's RH-calculator, this is how it
  ! is formed from SAM's esatw(T) and esati(T) at the time of training.
    implicit none
    real T
    real, parameter :: T0 = 273.16
    real, parameter :: T00 = 253.16
    real, external :: esatw_crm,esati_crm ! register functions from crm source.
    real :: omtmp,omega
    omtmp = (T-T00)/(T0-T00)
    omega = max(0.,min(1.,omtmp))
!tf.where(T>T0,eliq(T),tf.where(T<T00,eice(T),(omega*eliq(T)+(1-omega)*eice(T))))
    if (T .gt. T0) then
      tom_esat = esatw_crm(T)
    elseif (T .lt. T00) then
      tom_esat = esati_crm(T)
    else
      tom_esat = omega*esatw_crm(T) + (1-omega)*esati_crm(T)
    endif
  end

end module cloudbrain
#endif
