#define CBRAIN
#ifdef CBRAIN
#define BRAINDEBUG
!#define RHDEBUG
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
   integer :: i,k,ncol,ixcldice,ixcldliq,ii,kk,klev_crmtop
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

! deny any moisture activity in the stratosphere:
   do i=1,ncol
     call detect_tropopause(state%t(i,:),state%exner(i,:),state%zm(i,:),state%pmid(i,:),klev_crmtop)
     q_bctend(i,1:klev_crmtop) = 0.
     qc_bctend(i,1:klev_crmtop) = 0.
     qi_bctend(i,1:klev_crmtop) = 0.
   end do
! -- atmos positivity constraints ---- 
   if (doconstraints) then
   do i=1,ncol
     do k=1,pver
! deny activity in the ice phase where it is above freezing.
       if (state%t(i,k) .gt. 273.16) then
          qi_bctend(i,k) = 0.
! deny activitiy in the water phase where it is below freezing.
       elseif (state%t(i,k) .lt. 253.16) then
          qc_bctend(i,k) = 0.
       end if
!eliminate all activity in the water phase on top 10 levels:

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
  ! was done in the training environment (Caution: could be porting bugs here)
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
      tom_esat = tom_eliq(T)
    elseif (T .lt. T00) then
      tom_esat = tom_eice(T)
    else
      tom_esat = omega*tom_eliq(T) + (1.-omega)*tom_eice(T)
    endif
  end

  real function tom_eliq(T)
    implicit none
    real T
    real, parameter :: T0 = 273.16
    real, parameter :: cliq = -80. 
    real a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
       6.11239921, 0.443987641, 0.142986287e-1, &
       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
    real :: dt
    dt = max(cliq,T-T0)
    tom_eliq = 100.*(a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))))  
  end 


  real function tom_eice(T)
    implicit none
    real T
    real, parameter :: T0 = 273.16
    real a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/       
    real cice(6)
    real dt
    dt = T-T0
    cice(1) = 273.15
    cice(2) = 185.
    cice(3) = -100.
    cice(4) = 0.00763685
    cice(5) = 0.000151069
    cice(6) = 7.48215e-07
    if (T .gt. cice(1)) then
      tom_eice = tom_eliq(T)
    else if (T .le. cice(2)) then
      tom_eice = 100.*(cice(4) + max(cice(2),dt)*(cice(5)+max(cice(3),dt)*cice(6))) 
    else
      tom_eice = 100.*(a0 +dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))))
    end if
  end      


 subroutine detect_tropopause (t,exner,zmid,pmid,klev_crmtop)
   real(r8), intent(in) :: t(pver),exner(pver),zmid(pver),pmid(pver)
   integer, intent(out) :: klev_crmtop
   integer :: k
   real (r8) :: theta(pver),dthetadz

   do k=1,pver
     theta(k) = t(k)*exner(k)
   end do

   klev_crmtop = 1

   do k=2,pver-1
     dthetadz = (theta(k-1)-theta(k+1))/(zmid(k-1)-zmid(k+1))*1000. ! K/km
     ! assume theta in K and pmid in Pa then
     if (pmid(k) .le. 40000. .and. dthetadz > 10.) then
       klev_crmtop = k
     endif
   end do
  end subroutine detect_tropopause

end module cloudbrain
#endif
