#define BRAINDEBUG
module cloudbrain
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use cam_history,         only: outfld, addfld, add_default, phys_decomp
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

  subroutine neural_net (state,nn_solin,cam_in,ptend,cam_out)
    type(physics_state), intent(in)    :: state
    real(r8), intent(in)               :: nn_solin(pcols) 
    type(cam_in_t),intent(in)          :: cam_in
    type(physics_ptend),intent(in)     :: ptend            ! indivdual parameterization tendencies
    type(cam_out_t),     intent(inout) :: cam_out

    ! local variables
    real :: input(pcols,inputlength)
    real :: output(pcols,outputlength)
    integer :: i,k,ncol
    
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
      input(k) = (input(k) - inp_sub(k))/inp_div(k)
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

! INSERT wiring to ptend, cam_out
! From Tom's Slack the output vector is ordered (but need to check with Ankitesh and vs. output norms):
! TBCTEND=(TBC-TBP)/DT,QBCTEND=(QBC-QBP)/DT,CLDLIQBCTEND=(CLDLIQBC-CLDLIQBP)/DT,CLDICEBCTEND=(CLDICEBC-CLDICEBP)/DT,NN2L_FLWDS,
! NN2L_NETSW(misleading name as it's downwelling shortwave),NN2L_PRECC,NN2L_PRECSC,NN2L_SOLL,NN2L_SOLLD,NN2L_SOLS,NN2L_SOLSD)
! 4*pver + 8 = 4*26 +8 = 112
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
    if (masterproc)
       write (iulog,*) 'BRAINDEBUG read input norm inp_sub=', inp_sub(:)
       write (iulog,*) 'BRAINDEBUG read input norm inp_div=', inp_div(:)       
       write (iulog,*) 'BRAINDEBUG read output norm out_scale=', out_scale(:)       
    #endif

  end subroutine init_neural_net
    
end module cloudbrain
