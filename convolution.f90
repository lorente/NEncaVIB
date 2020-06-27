module convolution

USE declfft
USE declarations

implicit none

CONTAINS
!
! By making use of the convolution theorem
! we compute the convolution of two functions
! g(x) and h(x) to give f(x).
!

  subroutine convolution_fft ( N_dim, N_padded, step, g, h, matrix )
  integer, intent (in) :: N_dim, N_padded
  real (q), intent (in) :: step
  complex (qc), intent (in) :: g (N_dim), h (N_dim)
  complex (qc), intent (out) :: matrix (N_dim)
  complex (qc) :: g_t (N_padded), h_t(N_padded), matrix_t (N_padded)
  integer :: i, j, i_ini, i_fin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                 !
! We zero-pad g and h to make g_t and h_t                         !
!                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i_ini = int(N_padded/4); i_fin = 3*int(N_padded/4)
    g_t = 0._qc; H_t  = 0._qc; matrix_t = 0._qc
    g_t (i_ini:i_fin) = g(1:N_dim)
    h_t (i_ini:i_fin) = h(1:N_dim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                 !
! We use FFT ordering of indices, remember:                  !
! j = 0,1,2,...,N/2,-N/2,-N/2+1,...,-2,-1                         !
! or                                                              !
! i = N/2+1, N/2+2,...,N,1,2,...,N/2                              !
! where g(N/2)=g(-N/2) due to periodicity (here padded to zero)   !
! or g(N)=g(0)                                                    !
! this ordering is achieved using                                 !
! j = mod (i+N/2-1,N) + 1                                         !
! with i = 1, 2,..., N                                            !
!                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! FFT index ordering
  do i = 1, N_padded
    j = mod (i+N_padded/2-1,N_padded) + 1 
    matrix_t (i) = g_t (j)
  enddo

  call fft1d (matrix_t, N_padded)

   g_t = matrix_t

! FFT index ordering
  do i = 1, N_padded
    j = mod (i+N_padded/2-1,N_padded) + 1 
    matrix_t (i) = h_t (j)
  enddo

  call fft1d (matrix_t, N_padded)

   h_t = matrix_t

! product of functions is transformed back

   matrix_t = g_t*h_t

   call fft1d_r (matrix_t, N_padded)

! and we need to redo the index ordering
! FFT index ordering
    g_t = (0._q,0._q)
  do i = 1, N_padded
    j = mod (i+N_padded/2-1,N_padded) + 1 
    g_t (j) = matrix_t (i)
  enddo

! and FFT normalization:
 
   matrix(1:N_dim) = g_t(i_ini:i_fin) * step / N_padded

   

   return
  end subroutine convolution_fft

  subroutine correlation_fft ( N_dim, N_padded, step, g, h, matrix )
  integer, intent (in) :: N_dim, N_padded
  real (q), intent (in) :: step
  complex (qc), intent (in) :: g (N_dim), h (N_dim)
  complex (qc), intent (out) :: matrix (N_dim)
  complex (qc) :: g_t (N_padded), h_t(N_padded), matrix_t (N_padded)
  integer :: i, j, i_ini, i_fin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                 !
! We zero-pad g and h to make g_t and h_t                         !
!                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i_ini = int(N_padded/4); i_fin = 3*int(N_padded/4)
    g_t = 0._qc; H_t  = 0._qc; matrix_t = 0._qc
    g_t (i_ini:i_fin) = g(1:N_dim)
    h_t (i_ini:i_fin) = h(1:N_dim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                 !
! We use FFT ordering of indices, remember:                  !
! j = 0,1,2,...,N/2,-N/2,-N/2+1,...,-2,-1                         !
! or                                                              !
! i = N/2+1, N/2+2,...,N,1,2,...,N/2                              !
! where g(N/2)=g(-N/2) due to periodicity (here padded to zero)   !
! or g(N)=g(0)                                                    !
! this ordering is achieved using                                 !
! j = mod (i+N/2-1,N) + 1                                         !
! with i = 1, 2,..., N                                            !
!                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! FFT index ordering
  do i = 1, N_padded
    j = mod (i+N_padded/2-1,N_padded) + 1 
    matrix_t (i) = g_t (j)
  enddo

  call fft1d (matrix_t, N_padded)

   g_t = matrix_t

! FFT index ordering
  do i = 1, N_padded
    j = mod (i+N_padded/2-1,N_padded) + 1 
    matrix_t (i) = h_t (j)
  enddo

! call reversed fft (changed 4 November 2010)
  call fft1d_r (matrix_t, N_padded)

   h_t = matrix_t

! product of functions is transformed back

! call reversed fft (changed 4 November 2010)
!  matrix_t = g_t*conjg(h_t) only true if original real functions
   matrix_t = g_t*h_t

   call fft1d_r (matrix_t, N_padded)

! and we need to redo the index ordering
! FFT index ordering
    g_t = (0._q,0._q)
  do i = 1, N_padded
    j = mod (i+N_padded/2-1,N_padded) + 1 
    g_t (j) = matrix_t (i)
  enddo

! and FFT normalization:
 
   matrix(1:N_dim) = g_t(i_ini:i_fin) * step / N_padded

   

   return
   end subroutine correlation_fft

!!!!!!!!!!!!!!!!!!!!!!!!!
!                       !
!      forward FFT      !
!                       !
!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft1d (temp, n)
   integer   n
   complex (q) temp (0:n-1)

   call dfftw_plan_dft_1d (plan,n,temp,temp,FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   return
   end subroutine fft1d

!!!!!!!!!!!!!!!!!!!!!!!!!
!                       !
!     backward FFT      !
!                       !
!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft1d_r (temp, n)
   integer  n
   complex (q) temp (0:n-1)

   call dfftw_plan_dft_1d (plan,n,temp,temp,FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   return
   end subroutine fft1d_r
! invert indices for the convolution in Self_0
      subroutine invert_indices (N_omega, h)
      integer,      intent(in)  :: N_omega
      integer             :: i
      complex (qc) :: buffer (N_omega), h (N_omega)

! We assume the intergal in omega is symmetric about zero

      do i = 1, N_omega
       buffer (i) = h (N_omega-i+1)
      enddo

       h = buffer

      return
      end subroutine invert_indices


end module convolution
