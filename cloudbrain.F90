#define CBRAIN
#ifdef CBRAIN
#define BRAINDEBUG
!#define RHDEBUG
!#define CBRAIN_OCN_ONLY

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
use physics_buffer,   only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
use cam_history_support, only: pflds, fieldname_lenp2
use cam_abortutils, only: endrun

! -------- NEURAL-FORTRAN --------
! imports
use mod_kinds, only: ik, rk
use mod_network , only: network_type
! --------------------------------

  implicit none
  save 

  private
  ! Define variables for this entire module
  ! These are nameliest variables.
  ! If not specified in atm_in, the following defaults values are used.
  integer :: inputlength  = 108     ! length of NN input vector
  integer :: outputlength = 112     ! length of NN output vector
  logical :: input_rh     = .false.  ! toggle to switch from q --> RH input
  logical :: cb_use_input_prectm1 = .false.  ! use previous timestep PRECT for input variable 
  character(len=256)    :: cb_fkb_model   ! absolute filepath for a fkb model txt file
  character(len=256)    :: cb_inp_sub     ! absolute filepath for a inp_sub.txt
  character(len=256)    :: cb_inp_div     ! absolute filepath for a inp_div.txt
  character(len=256)    :: cb_out_scale   ! absolute filepath for a out_scale.txt
  logical :: cb_partial_coupling  = .false.
  character(len=fieldname_lenp2) :: cb_partial_coupling_vars(pflds)

  type(network_type), allocatable :: cloudbrain_net(:)
  real(r8), allocatable :: inp_sub(:)
  real(r8), allocatable :: inp_div(:)
  real(r8), allocatable :: out_scale(:)

  logical :: cb_do_ensemble  = .false.
  integer :: cb_ens_size
  integer :: max_nn_ens = 100 ! Max. ensemble size is arbitrarily set to 100.
  character(len=256), allocatable :: cb_ens_fkb_model_list(:)

  integer :: cb_random_ens_size = 0

  ! local
  integer, allocatable :: ens_ind_shuffled(:)
  logical :: cb_do_random_ensemble = .false.

  public neural_net, init_neural_net, cbrain_readnl, &
         cb_partial_coupling, cb_partial_coupling_vars
  
contains

  subroutine neural_net (state,nn_solin,cam_in,ztodt,ptend,cam_out,pbuf)
 ! note state is meant to have the "BP" state saved earlier. 

   implicit none

   type(physics_state), intent(in)    :: state
   real(r8), intent(in)               :: nn_solin(pcols) 
   type(cam_in_t),intent(in)          :: cam_in
   real(r8), intent(in) :: ztodt 
   type(physics_ptend),intent(out)    :: ptend            ! indivdual parameterization tendencies
   type(cam_out_t),     intent(inout) :: cam_out  ! SY: changed to inout to use variables from a previous time step (e.g., PRECT)
   type(physics_buffer_desc), pointer  :: pbuf(:) ! SY: for precip variables.

   ! SY: for random ensemble averaging
   integer, external :: shuffled_1d 

    ! local variables
   real(r8) :: input(pcols,inputlength)
   real(r8) :: output(pcols,outputlength)
   integer :: i,k,ncol,ixcldice,ixcldliq,ii,kk,klev_crmtop,kens
   real (r8) :: s_bctend(pcols,pver), q_bctend(pcols,pver), qc_bctend(pcols,pver), qi_bctend(pcols,pver), qafter, safter
   logical :: doconstraints
   logical ::  lq(pcnst)
   real(r8) :: rh_loc
   integer :: pvert,ntrim
   integer :: prec_dp_idx, snow_dp_idx
   ! convective precipitation variables
   real(r8), pointer :: prec_dp(:)                ! total precip rate [m/s]
   real(r8), pointer :: snow_dp(:)                ! snow precip rate [m/s]

   ntrim = 0 !sungduk: cam_nz=26, crm_nz=24, However, NN input variable has nz=26.
   pvert = pver-ntrim ! after trimming.

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
  ! Gunnar's ANN (2022DEC)
  ! INPUT: ['QBP', 'TBP','PS', 'SOLIN', 'SHFLX','LHFLX','PRECTt-dt','CLDLIQBP','CLDICEBP']
  ! Gunnar's ANN (2023FEB)
  ! INPUT: ['QBP', 'TBP','PS', 'SOLIN', 'SHFLX', 'LHFLX','CLDLIQBP','CLDICEBP']
    if (input_rh) then
       do i = 1,ncol
         do k=ntrim+1,pver
           ! Port of tom's RH =  Rv*p*qv/(R*esat(T))
           rh_loc = 461.*state%pmid(i,k)*state%q(i,k,1)/(287.*tom_esat(state%t(i,k))) ! note function tom_esat below refercing SAM's sat.F90
#ifdef RHDEBUG
           if (masterproc) then
             write (iulog,*) 'RHDEBUG:p,q,T,RH=',state%pmid(i,k),state%q(i,k,1),state%t(i,k),rh_loc
           endif
#endif
           input(i,k-ntrim) = rh_loc
         end do
       end do
    else
      do i=1,ncol 
        do k=ntrim+1,pver
          input(i,k-ntrim) = state%q(i,k,1) ! specific humidity input
        end do
      end do
    endif
    input(:ncol,(pvert+1):(2*pvert))     = state%t(:ncol,(ntrim+1):pver) ! TBP
    input(:ncol,(2*pvert+1))             = state%ps(:ncol) ! PS
    input(:ncol,(2*pvert+2))             = nn_solin(:ncol) ! SOLIN / WARNING this is being lazily mined from part of SP solution... should be avoidable in future when bypassing SP totally but will take work.
    input(:ncol,(2*pvert+3))             = cam_in%shf(:ncol) ! SHFLX
    input(:ncol,(2*pvert+4))             = cam_in%lhf(:ncol) ! LHFLX
    input(:ncol,(2*pvert+5):(3*pvert+4)) = state%q(:ncol,(ntrim+1):pver,ixcldliq) ! CLDLIQBP
    input(:ncol,(3*pvert+5):(4*pvert+4)) = state%q(:ncol,(ntrim+1):pver,ixcldice) ! CLDICEBP
    if (cb_use_input_prectm1) then
       input(:ncol,(2*pvert+5))             = cam_out%precc(:ncol)  + cam_out%precl(:ncol) ! PRECTt-dt
       input(:ncol,(2*pvert+6):(3*pvert+5)) = state%q(:ncol,(ntrim+1):pver,ixcldliq) ! CLDLIQBP
       input(:ncol,(3*pvert+6):(4*pvert+5)) = state%q(:ncol,(ntrim+1):pver,ixcldice) ! CLDICEBP
    else
       input(:ncol,(2*pvert+5):(3*pvert+4)) = state%q(:ncol,(ntrim+1):pver,ixcldliq) ! CLDLIQBP
       input(:ncol,(3*pvert+5):(4*pvert+4)) = state%q(:ncol,(ntrim+1):pver,ixcldice) ! CLDICEBP
    endif


! Tue Jan 24 13:28:43 CST 2023
! Sungduk 
#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG TBP=',state%t(1,(ntrim+1):pver)
        write (iulog,*) 'BRAINDEBUG PS=',state%ps(1)
        write (iulog,*) 'BRAINDEBUG SOLIN=',nn_solin(1)
        write (iulog,*) 'BRAINDEBUG SHFLX=',cam_in%shf(1)
        write (iulog,*) 'BRAINDEBUG LHFLX=',cam_in%lhf(1)
        if (cb_use_input_prectm1) then
          write (iulog,*) 'BRAINDEBUG PRECTm1=',cam_out%precc(1) + cam_out%precl(1)
        endif
        write (iulog,*) 'BRAINDEBUG CLDLIQBP=', state%q(1,(ntrim+1):pver,ixcldliq)
        write (iulog,*) 'BRAINDEBUG CLDICEBP=', state%q(1,(ntrim+1):pver,ixcldice)
      endif
#endif# 

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
      if (cb_do_ensemble) then
        output(i,:) = 0.
        !! Random ensemble averaging
        if (cb_do_random_ensemble) then
          ens_ind_shuffled = shuffle_1d(ens_ind_shuffled) ! randomly shuffle ens indices
          do kens = 1,cb_ens_size
            if ( ens_ind_shuffled(kens) .le. cb_random_ens_size ) then
              output(i,:) = output(i,:) +  (1._r8/cb_random_ens_size) * cloudbrain_net(kens) % output(input(i,:))
            end if
          enddo
#ifdef BRAINDEBUG
          if (masterproc .and. i.eq.1) then
            write (iulog,*) 'BRAINDEBUG  random ensemble model IDs = ',output(1,1:cb_random_ens_size)
          endif
#endif
        !! All ensemble averaging
        else
          do kens = 1,cb_ens_size
            output(i,:) = output(i,:) +  (1._r8/cb_ens_size) * cloudbrain_net(kens) % output(input(i,:))
          enddo
        endif
      !! Using a single model
      else ! cb_do_ensemble
        output(i,:) = cloudbrain_net(1) % output(input(i,:))
      endif
    end do
#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG output = ',output(1,:)
      endif
#endif

   ! Manually applying ReLU activation for positive-definite variables
   do i=1,ncol
     do k=4*pvert+1,4*pvert+8
       output(i,k) = max(output(i,k), 0.)
     end do
     k=4*pvert+3
     output(i,k) = max(output(i,k), tiny(output(i,k))) ! flwds
                                                       ! preventing flwds==0 error
   end do
#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG output after ReLU = ',output(1,:)
      endif
#endif

   ! output normalization (un-weighting, really).
   do i=1,ncol
     do k=1,outputlength
      output(i,k) = output(i,k) / out_scale(k)
     end do
   end do
#ifdef BRAINDEBUG
      if (masterproc) then
        write (iulog,*) 'BRAINDEBUG output post scale = ',output(1,:)
      endif
#endif

! OUTPUT:
! ['QBCTEND','TBCTEND','CLDLIQBCTEND','CLDICEBCTEND','PREC_CRM_SNOW','PREC_CRM','NN2L_FLWDS','NN2L_DOWN_SW','NN2L_SOLL','NN2L_SOLLD','NN2L_SOLS','NN2L_SOLSD']

! ---------- 1. NN output to atmosphere forcing --------
! ['QBCTEND','TBCTEND','CLDLIQBCTEND','CLDICEBCTEND']
! ! don't use CRM tendencies from two crm top levels
   q_bctend (:ncol,ntrim+1:pver) = output(:ncol,1:pvert) ! kg/kg/s 
   s_bctend (:ncol,ntrim+1:pver) = output(:ncol,(pvert+1)  :(2*pvert))*cpair ! K/s --> J/kg/s (ptend expects that)
   qc_bctend(:ncol,ntrim+1:pver) = output(:ncol,(2*pvert+1):(3*pvert)) ! kg/kg/s 
   qi_bctend(:ncol,ntrim+1:pver) = output(:ncol,(3*pvert+1):(4*pvert)) ! kg/kg/s 

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
!!! Sungduk: these are original version.
!!!          It works, but I wrote it again to add 'ocean only coupling' option
!!!          (#CBRAIN_OCN_ONLY)
!!! ! ['PRECT','PREC_CRM_SNOW','PREC_CRM','NN2L_FLWDS','NN2L_DOWN_SW','NN2L_SOLL','NN2L_SOLLD','NN2L_SOLS','NN2L_SOLSD']
!!!    ! These are the cam_out members that are not assigned in cam_export
!!!    cam_out%flwds = output(:ncol,4*pvert+3)
!!!    cam_out%netsw = output(:ncol,4*pvert+4)
!!!    cam_out%soll  = output(:ncol,4*pvert+5)
!!!    cam_out%solld = output(:ncol,4*pvert+6)
!!!    cam_out%sols  = output(:ncol,4*pvert+7)
!!!    cam_out%solsd = output(:ncol,4*pvert+8)
!!! 
!!!    ! These are the cam_out members that are assigned in cam_export,
!!!    ! and so saved to pbuf, instead.
!!!    ! SY: Note that this uses SPCAM's pbuf register.
!!!    !     Once, SP is entirely removed, we still need to call crm_physics_register().
!!!    prec_dp_idx = pbuf_get_index('PREC_DP', errcode=i) ! Query physics buffer index
!!!    snow_dp_idx = pbuf_get_index('SNOW_DP', errcode=i)
!!!    call pbuf_get_field(pbuf, prec_dp_idx, prec_dp) ! Associate pointers withphysics buffer fields
!!!    call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
!!!    prec_dp(:ncol) = output(:ncol,4*pvert+2)   ! PREC_CRM
!!!    snow_dp(:ncol) = output(:ncol,4*pvert+1)   ! PREC_CRM_SNOW

   prec_dp_idx = pbuf_get_index('PREC_DP', errcode=i) ! Query physics buffer index
   snow_dp_idx = pbuf_get_index('SNOW_DP', errcode=i)
   call pbuf_get_field(pbuf, prec_dp_idx, prec_dp) ! Associate pointers withphysics buffer fields
   call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
   do i = 1,ncol
! SY: debugging
!     allowing surface coupling over ocean only
#ifdef CBRAIN_OCN_ONLY 
     if (cam_in%ocnfrac(i) .eq. 1.0_r8) then
#endif
       cam_out%flwds(i) = output(i,4*pvert+3)
       cam_out%netsw(i) = output(i,4*pvert+4)
       cam_out%soll(i)  = output(i,4*pvert+5)
       cam_out%solld(i) = output(i,4*pvert+6)
       cam_out%sols(i)  = output(i,4*pvert+7)
       cam_out%solsd(i) = output(i,4*pvert+8)

       prec_dp(i) = output(i,4*pvert+2)   ! PREC_CRM
       snow_dp(i) = output(i,4*pvert+1)   ! PREC_CRM_SNOW
#ifdef CBRAIN_OCN_ONLY
     end if
#endif
   end do 

end subroutine neural_net

  subroutine init_neural_net()

    implicit none

    integer :: i, k

    allocate(inp_sub (inputlength))
    allocate(inp_div (inputlength))
    allocate(out_scale (outputlength))
    
    ! ens-mean inference
    if (cb_do_ensemble) then
       write (iulog,*) 'CLOUDBRAIN: Ensemble is turned on with Ensemble size  ', cb_ens_size
       allocate(cloudbrain_net (cb_ens_size))
       do i = 1,cb_ens_size
          call cloudbrain_net(i) %load(cb_ens_fkb_model_list(i))
          write (iulog,*) 'CLOUDBRAIN: Ensemble fkb model (', i, ') : ', trim(cb_ens_fkb_model_list(i))
       enddo

       ! random ensemble
       if (cb_random_ens_size .ge. 1) then
          write (iulog,*) 'CLOUDBRAIN: Random ensemble averaging with N = ', cb_random_ens_size
          if (cb_random_ens_size .le. cb_ens_size) then
             allocate(ens_ind_shuffled(cb_ens_size))
             ens_ind_shuffled = (/ (k, k=1, cb_ens_size) /)
             cb_do_random_ensemble = .true.
          else
             call endrun("init_neural_net error: cb_random_ens_size should be less than or equal to cb_ens_size")
          endif

       endif

    ! single model inference
    else
       allocate(cloudbrain_net (1))
       call cloudbrain_net(1) %load(cb_fkb_model)
       write (iulog,*) 'CLOUDBRAIN: loaded network from txt file, ', trim(cb_fkb_model)
    endif

    open (unit=555,file=cb_inp_sub,status='old',action='read')
    read(555,*) inp_sub(:)
    close (555)
    if (masterproc) then
       write (iulog,*) 'CLOUDBRAIN: loaded inp_sub.txt, ', trim(cb_inp_sub)
    endif

    open (unit=555,file=cb_inp_div,status='old',action='read')
    read(555,*) inp_div(:)
    close (555)
    if (masterproc) then
       write (iulog,*) 'CLOUDBRAIN: loaded inp_div.txt, ', trim(cb_inp_div)
    endif

    open (unit=555,file=cb_out_scale,status='old',action='read')
    read(555,*) out_scale(:)
    close (555)
    if (masterproc) then
       write (iulog,*) 'CLOUDBRAIN: loaded out_scale.txt, ', trim(cb_out_scale)
    endif

#ifdef BRAINDEBUG
    if (masterproc) then
       write (iulog,*) 'BRAINDEBUG read input norm inp_sub=', inp_sub(:)
       write (iulog,*) 'BRAINDEBUG read input norm inp_div=', inp_div(:)       
       write (iulog,*) 'BRAINDEBUG read output norm out_scale=', out_scale(:)       
    endif
#endif

  end subroutine init_neural_net

  real(r8) function tom_esat(T) 
  ! For consistency with the python port of Tom's RH-calculator, this is how it
  ! was done in the training environment (Caution: could be porting bugs here)
    implicit none
    real(r8) T
    real(r8), parameter :: T0 = 273.16
    real(r8), parameter :: T00 = 253.16
    real(r8), external :: esatw_crm,esati_crm ! register functions from crm source.
    real(r8) :: omtmp,omega
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

  real(r8) function tom_eliq(T)
    implicit none
    real(r8) T
    real(r8), parameter :: T0 = 273.16
    real(r8), parameter :: cliq = -80. 
    real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
       6.11239921, 0.443987641, 0.142986287e-1, &
       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
    real(r8) :: dt
    dt = max(cliq,T-T0)
    tom_eliq = 100.*(a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))))  
  end 


  real(r8) function tom_eice(T)
    implicit none
    real(r8) T
    real(r8), parameter :: T0 = 273.16
    real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/       
    real(r8) cice(6)
    real(r8) dt
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

  ! Read namelist variables.
  subroutine cbrain_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr, f
      character(len=*), parameter :: subname = 'cbrain_readnl'
      
      namelist /cbrain_nl/ inputlength, outputlength, input_rh, &
                           cb_fkb_model, &
                           cb_inp_sub, cb_inp_div, cb_out_scale, &
                           cb_partial_coupling, cb_partial_coupling_vars,&
                           cb_use_input_prectm1, &
                           cb_do_ensemble, cb_ens_size, cb_ens_fkb_model_list, &
                           cb_random_ens_size

      ! Initialize 'cb_partial_coupling_vars'
      do f = 1, pflds
        cb_partial_coupling_vars(f) = ' '
      end do

      ! Initialize 'cb_ens_fkb_model_list'
      allocate(cb_ens_fkb_model_list(max_nn_ens))
      do f = 1, max_nn_ens
        cb_ens_fkb_model_list(f) = ' '
      end do

      ierr = 0
      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'cbrain_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, cbrain_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(inputlength,  1,                 mpiint,  0, mpicom)
      call mpibcast(outputlength, 1,                 mpiint,  0, mpicom)
      call mpibcast(input_rh,     1,                 mpilog,  0, mpicom)
      call mpibcast(cb_use_input_prectm1,1,          mpilog,  0, mpicom)
      call mpibcast(cb_fkb_model, len(cb_fkb_model), mpichar, 0, mpicom)
      call mpibcast(cb_inp_sub,   len(cb_inp_sub),   mpichar, 0, mpicom)
      call mpibcast(cb_inp_div,   len(cb_inp_div),   mpichar, 0, mpicom)
      call mpibcast(cb_out_scale, len(cb_out_scale), mpichar, 0, mpicom)
      call mpibcast(cb_partial_coupling, 1,          mpilog,  0, mpicom)
      call mpibcast(cb_partial_coupling_vars, len(cb_partial_coupling_vars(1))*pflds, mpichar, 0, mpicom, ierr)
      call mpibcast(cb_do_ensemble, 1,               mpilog,  0, mpicom)
      call mpibcast(cb_ens_size,    1,               mpiint,  0, mpicom)
      call mpibcast(cb_ens_fkb_model_list,    len(cb_ens_fkb_model_list(1))*max_nn_ens, mpichar, 0, mpicom, ierr)
      call mpibcast(cb_random_ens_size,    1,        mpiint,  0, mpicom)
      if (ierr /= 0) then
         call endrun(subname // ':: ERROR broadcasting namelist variable cb_partial_coupling_vars')
      end if
#endif

   end subroutine cbrain_readnl

  function shuffle_1d(array_1d) result(array_shuffled)
  ! Shuffling the entries of 1-d INTEGER array
  ! (using the Knuth shuffle algorithm: https://en.wikipedia.org/wiki/Fisher–Yates_shuffle)
    implicit none
    integer, intent(in) :: array_1d(:)
    integer :: array_shuffled(size(array_1d))
    integer :: j, k, tmp
    real :: u

    array_shuffled(:) = array_1d(:)

    j    = size(array_1d)
    do while (j > 1)
       call random_seed
       call random_number(u)
       k = 1 + FLOOR(j * u)
       tmp = array_shuffled(j)
       array_shuffled(j) = array_shuffled(k)
       array_shuffled(k) = tmp
       j = j -1
    end do
  end function shuffle_1d

end module cloudbrain
#endif
