#include <misc.h>
#include <params.h>
!#define BRAINDEBUG

module cloudbrain
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use history,         only: outfld, addfld, add_default, phys_decomp
use physconst,       only: gravit,cpair,latvap,latice
use pmgrid, only: masterproc
!use runtime_opts, only: nn_nint, inputlength, outputlength, activation_type, width

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
  integer, parameter :: nn_nint = 8
  integer, parameter :: inputlength = 94
  integer, parameter :: outputlength = 65
  integer, parameter :: activation_type = 1
  integer, parameter :: width = 256

#ifdef NEURALLIB
#ifdef ENSEMBLE
  real(rk) :: noise = 0.0
  type(ensemble_type) :: cloudbrain_ensemble
#else
  type(network_type) :: cloudbrain_net
  integer(ik) :: fileunit, num_layers
  integer(ik) :: n
#endif
#else
  real :: weights_inp(width, inputlength)
  real :: bias_inp(width)
  real :: weights_int(nn_nint, width, width)
  real :: bias_int(nn_nint, width)
  real :: weights_out(outputlength, width)
  real :: bias_out(outputlength)
#endif

  real :: inp_sub(inputlength)
  real :: inp_div(inputlength)
  real :: out_scale(outputlength)

  public neural_net, init_keras_matrices, init_keras_norm
  contains

  subroutine neural_net (QBP, TBP, VBP, PS, SOLIN, SHFLX, LHFLX, &
                         PHQ, TPHYSTND, FSNT, FSNS, FLNT, FLNS, PRECT, &
                         icol)
    ! PNAS version: First row = inputs, second row = outputs
    ! icol is used for debugging to only output one colum
    ! Allocate inputs
    real(r8), intent(in) :: QBP(:)
    real(r8), intent(in) :: TBP(:)   
    real(r8), intent(in) :: VBP(:)
    real(r8), intent(in) :: PS
    real(r8), intent(in) :: SOLIN
    real(r8), intent(in) :: SHFLX
    real(r8), intent(in) :: LHFLX
    ! Allocate outputs
    real(r8), intent(out) :: PHQ(:)
    real(r8), intent(out) :: TPHYSTND(:)
    real(r8), intent(out) :: FSNT
    real(r8), intent(out) :: FSNS
    real(r8), intent(out) :: FLNT
    real(r8), intent(out) :: FLNS
    real(r8), intent(out) :: PRECT
    ! Allocate utilities
#ifdef NEURALLIB
    real(rk) :: input(inputlength),x1(width), x2(width)
#else
    real(r8) :: input(inputlength),x1(width), x2(width)
#endif
    real(r8) :: output (outputlength)
    integer :: k, nlev, n
    integer, intent(in) :: icol

    ! 1. Concatenate input vector to neural network
    nlev=30
#ifdef NEURALLIB
    input(1:nlev) = TBP(:) !QBP(:) 
    input((nlev+1):2*nlev) = QBP(:) !TBP(:)
#else
    input(1:nlev) = QBP(:)
    input((nlev+1):2*nlev) = TBP(:)
#endif
    input((2*nlev+1):3*nlev)=VBP(:)
    input(3*nlev+1) = PS
    input(3*nlev+2) = SOLIN
    input(3*nlev+3) = SHFLX
    input(3*nlev+4) = LHFLX
#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG input pre norm=',input
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

! 3. Neural network matrix multiplications and activations
#ifdef NEURALLIB
#ifdef ENSEMBLE
    output = cloudbrain_ensemble % average(input)
#else
    ! use neural fortran library
    output = cloudbrain_net % output(input)
#endif
#else
    ! hard code layers
    ! 3.1 1st layer: input length-->256.
    call matmul(input, x1, inputlength, width, weights_inp, bias_inp, activation_type)
#ifdef BRAINDEBUG
        if (masterproc .and. icol .eq. 1) then
          write (6,*) 'BRAINDEBUG layer = ',1
          write (6,*) 'BRAINDEBUG output = ',x1
        endif
#endif
    ! 3.2 Intermediate layers
    do n=1,nn_nint
      call matmul(x1, x2, width, width, weights_int(n, :, :), bias_int(n, :), activation_type)
      x1(:) = x2(:)
#ifdef BRAINDEBUG
        if (masterproc .and. icol .eq. 1) then
          write (6,*) 'BRAINDEBUG layer = ',n+1
          write (6,*) 'BRAINDEBUG output = ',x1
        endif
#endif
    end do
    ! 3.3 Last layer with linear activation
    call matmul(x1, output, width, outputlength, weights_out, bias_out, 0)
#endif

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG output = ',output
      endif
#endif

    ! 4. Unnormalize output
    do k=1,outputlength
      output(k) = output(k) / out_scale(k)
    end do

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG out post scale = ',output
      endif
#endif

    ! 5. Split output into components
#ifdef NEURALLIB
    TPHYSTND(:) =      output(1:nlev) * cpair! JORDAN SWAPPED PHQ(:)
    PHQ(:) = output((nlev+1):2*nlev)! This is still the wrong unit, needs to be converted to W/m^2
#else
    PHQ(:) = output(1:nlev) 
    TPHYSTND(:) = output((nlev+1):2*nlev) * cpair
#endif
    FSNT =        output(2*nlev+1)
    FSNS =        output(2*nlev+2)
    FLNT =        output(2*nlev+3)
    FLNS =        output(2*nlev+4)
    PRECT =       output(2*nlev+5)

  end subroutine neural_net

#ifdef NEURALLIB
#else
  subroutine matmul(inp, out, len_in, len_out, weights, bias, act_type)
    ! Also do LeakyReLU in here
    ! Neural network matrix multiplication
    ! Activation type: 0 = linear, 1 = LeakyReLU
    integer, intent(in) :: len_in, len_out, act_type
    real, intent(in)    :: inp(len_in)
    real, intent(out)   :: out(len_out)
    real, intent(in)    :: bias(len_out)
    real, intent(in)    :: weights(len_out,len_in)
    integer             :: k, j
  
    out(:) = 0.
    do k=1,len_out
      do j=1,len_in
        out(k) = out(k) + weights(k,j)* inp(j)
      end do
      out(k) = out(k) + bias(k)
    end do
  
    if (act_type .eq. 1) then
      call leaky_relu(out, len_out, 0.3)
    end if
  
  end subroutine matmul
  
  
  subroutine leaky_relu(x, len, alpha)
    real, intent(inout)    :: x(len)
    integer, intent(in)    :: len
    real, intent(in)       :: alpha
    integer                :: k
    do k=1,len
      x(k) = max(alpha * x(k), x(k))  ! Leaky ReLU
    end do
  end subroutine leaky_relu
#endif


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
