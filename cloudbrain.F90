module cloudbrain
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use cam_history,         only: outfld, addfld, add_default, phys_decomp
use physconst,       only: gravit,cpair,latvap,latice
use spmd_utils, only: masterproc
use camsrfexch,       only: cam_out_t, cam_in_t
use constituents,    only: cnst_get_ind

!use runtime_opts, only: nn_nint, inputlength, outputlength, activation_type, width
use physics_types,    only: physics_state,  physics_ptend

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
  integer, parameter :: outputlength = 65

  type(network_type) :: cloudbrain_net
  integer(ik) :: fileunit, num_layers
  integer(ik) :: n

  real :: inp_sub(inputlength)
  real :: inp_div(inputlength)
  real :: out_sub(outputlength)
  real :: out_div(outputlength)

  public neural_net, init_keras_matrices, init_keras_norm
  contains

  subroutine neural_net (state,nn_solin,cam_in,ptend,cam_out)
    type(physics_state), intent(in)    :: state
    real(r8), intent(in)               :: nn_solin(pcols) 
    type(cam_in_t),intent(in)          :: cam_in
    type(physics_ptend),intent(in)     :: ptend            ! indivdual parameterization tendencies
    type(cam_out_t),     intent(inout) :: cam_out

    ! local variables
    real :: input(pcols,4*pver+4)
    integer :: ncol
    
    ncol  = state%ncol
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
  
    input(:ncol,1:pver) = state%t(:ncol,:pver)
    input(:ncol,(pver+1):(2*pver)) = state%q(:ncol,:pver,1)
    input(:ncol,(2*pver+1):(3*pver)) = state%q(:ncol,:pver,ixcldliq)
    input(:ncol,(3*pver+1):(4*pver)) = state%q(:ncol,:pver,ixcldice)
    input(:ncol,(4*pver+1)) = state%ps(:ncol)
    input(:ncol,(4*pver+2)) = nn_solin(:ncol) ! WARNING this is being lazily mined from part of SP solution... should be avoidable in future when bypassing SP totally but will take work.
    input(:ncol,(4*pver+3)) = cam_in%shf(:ncol)
    input(:ncol,(4*pver+4)) = cam_in%lhf(:ncol) 
 
    real(rk) :: input(inputlength)!,x1(width), x2(width)
    real(r8) :: output (outputlength)
    integer :: k, i

#ifdef BRAINDEBUG
      if (masterproc) then
        write (6,*) 'BRAINDEBUG input pre norm=',input(1,:)
      endif
#endif

    ! 2. Normalize input
    do k=1,inputlength
      input(k) = (input(k) - inp_sub(k))/inp_div(k)
    end do
#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG input post norm=',input
      endif
#endif

    do i=1,ncol
      output(i,:) = cloudbrain_net % output(input(i,:))
    end do

#ifdef BRAINDEBUG
      if (masterproc) then
        write (6,*) 'BRAINDEBUG output = ',output(1,:)
      endif
#endif

! INSERT output normalization

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG out post scale = ',output
      endif
#endif

! INSERT wiring to ptend, cam_out

  end subroutine neural_net

  subroutine init_keras_matrices()    
#ifdef NEURALLIB
#ifdef ENSEMBLE
    write (6,*) '------- NEURAL-FORTRAN: ensemble loading -------'
    cloudbrain_ensemble = ensemble_type('./Models/', noise)
    write (6,*) '------- NEURAL-FORTRAN: ensemble loaded -------'
#else
    call cloudbrain_net % load('./keras_matrices/model.txt')
    write (6,*) '------- NEURAL-FORTRAN: loaded network from txt file -------'
#endif
#else
#ifdef BRAINDEBUG
    do n = 1, size(cloudbrain_net % layers)
        write (6,*) 'BRAINDEBUG layer weight', n, cloudbrain_net % layers(n) % w(1$
    end do
#endif
    ! load weights and bias for each layer
    integer :: n, k, ios
    character(len=1) :: str
    character(len=2) :: str10  ! Yeah, I don't know how to deal with one vs two characters
    character(len=22) :: pref
    character(len=9) :: suf_bias
    character(len=11) :: suf_kernel

    pref = './keras_matrices/layer'
    suf_bias = '_bias.txt'
    suf_kernel = '_kernel.txt'

    ! 1. Input layer
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: reading layer1_bias'
    endif
    open (unit=555,file='./keras_matrices/layer1_bias.txt',status='old',action='read',iostat=ios)
    if (ios .ne. 0) then
      write (6,*) 'CLOUDBRAIN keras matrices unable to load, abort.'
      stop
    endif
    read(555,*) bias_inp(:)
    close (555)
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: reading layer1_kernel'
    endif
    open (unit=555,file='./keras_matrices/layer1_kernel.txt',status='old',action='read') 
    do k=1,width
      read(555,*) weights_inp(k,:)
    end do
    close (555)
#ifdef BRAINDEBUG
      if (masterproc) then
        write (6,*) 'BRAINDEBUG weights_inp = ',weights_inp
        write (6,*) 'BRAINDEBUG bias_inp = ',bias_inp
      endif
#endif

    ! 2. Intermediate layers
    do n=1,nn_nint
      write (str, '(I1)') n+1
      if (masterproc) then
        write (6,*) 'CLOUDBRAIN: reading layer*_bias', n+1, pref//trim(str)//suf_bias
      endif
      open (unit=555, file=pref//trim(str)//suf_bias, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
        write (6,*) 'CLOUDBRAIN keras matrices unable to load, abort.'
        stop
      endif
      read(555,*) bias_int(n, :)
      close (555)
      if (masterproc) then
        write (6,*) 'CLOUDBRAIN: reading layer*_kernel', n+1, pref//trim(str)//suf_kernel
      endif
      open (unit=555, file=pref//trim(str)//suf_kernel, status='old', action='read') 
      do k=1,width
        read(555,*) weights_int(n, k,:)
      end do
      close (555)
#ifdef BRAINDEBUG
        if (masterproc) then
          write (6,*) 'BRAINDEBUG weights_* = ',n+1, weights_int(n, k,:)
          write (6,*) 'BRAINDEBUG bias_* = ',n+1, bias_int(n, :)
        endif
#endif
    end do

    ! 3. Output layer
    if (nn_nint+2 .lt. 10) then
      write (str, '(I1)') nn_nint+2
      if (masterproc) then
        write (6,*) 'CLOUDBRAIN: reading layer*_bias = output', nn_nint+2, pref//trim(str)//suf_bias
      endif
      open (unit=555, file=pref//trim(str)//suf_bias, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
        write (6,*) 'CLOUDBRAIN keras matrices unable to load, abort.'
        stop
      endif
      read(555,*) bias_out(:)
      close (555)
      if (masterproc) then
        write (6,*) 'CLOUDBRAIN: reading layer*_kernel = output', nn_nint+2, pref//trim(str)//suf_kernel
      endif
      open (unit=555, file=pref//trim(str)//suf_kernel, status='old', action='read')
    else
      write (str10, '(I2)') nn_nint+2
      if (masterproc) then
        write (6,*) 'CLOUDBRAIN: reading layer*_bias = output', nn_nint+2, pref//trim(str10)//suf_bias
      endif
      open (unit=555, file=pref//trim(str10)//suf_bias, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
        write (6,*) 'CLOUDBRAIN keras matrices unable to load, abort.'
        stop
      endif
      read(555,*) bias_out(:)
      close (555)
      if (masterproc) then
        write (6,*) 'CLOUDBRAIN: reading layer*_kernel = output', nn_nint+2, pref//trim(str10)//suf_kernel
      endif
      open (unit=555, file=pref//trim(str10)//suf_kernel, status='old', action='read')
    endif
    do k=1,outputlength
      read(555,*) weights_out(k,:)
    end do
    close (555)
#ifdef BRAINDEBUG
      if (masterproc) then
        write (6,*) 'BRAINDEBUG weights_out = ',weights_out
        write (6,*) 'BRAINDEBUG bias_out = ',bias_out
      endif
#endif
#endif

  end subroutine init_keras_matrices
    

subroutine init_keras_norm()

  ! 1. Read sub
  if (masterproc) then
    write (6,*) 'CLOUDBRAIN: reading inp_sub'
  endif
  open (unit=555,file='./keras_matrices/inp_sub.txt',status='old',action='read')
  read(555,*) inp_sub(:)
  close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG inp_sub = ',inp_sub
    endif
#endif

  ! 2. Read div
  if (masterproc) then
    write (6,*) 'CLOUDBRAIN: reading inp_div'
  endif
  open (unit=555,file='./keras_matrices/inp_div.txt',status='old',action='read')
  read(555,*) inp_div(:)
  close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG inp_div = ',inp_div
    endif
#endif

  ! 3. Read out_scale
  if (masterproc) then
    write (6,*) 'CLOUDBRAIN: reading out_scale'
  endif
  open (unit=555,file='./keras_matrices/out_scale.txt',status='old',action='read')
  read(555,*) out_scale(:)
  close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG out_scale = ',out_scale
    endif
#endif

  end subroutine init_keras_norm

end module cloudbrain
