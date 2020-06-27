module  io ! input output module
Use declarations
Use used_files

! input output module
! with warning and error messages
!

CONTAINS

! general input file with file names

      subroutine  input (omega, step_omega, Temperature, eta, N_spin, &
      &        tol, Max_loop_index, B_au, gamma_file, output_file, &
      &        impurity_hamiltonian_file, Bias_ini, Bias_fin, N_bias, &
      &        bias_fraction, N_vib, Freq, W, vib_loop )
      implicit none

      character ( len = 100 ) :: impurity_hamiltonian_file
      character ( len = 100 ) :: gamma_file
      character ( len = 100 ) :: output_file

      integer :: ios, N_vib
      integer :: Max_loop_index, N_spin, N_omega, N_bias
      real (q) :: eta_input, omega_ini, step_omega, tol
      real (q) :: B_au, B_field, gyro, Temperature
      real (q) :: Bias_ini, Bias_fin, bias_fraction
      real (q), dimension (:), allocatable, intent (out):: omega, Freq, W
      complex (qc), intent (out) :: eta
      logical, intent (out) :: vib_loop

      ios = 0
      open (unit_input, file='WNnca.input', status='old')
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) omega_ini !eV from E_fermi
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) step_omega !eV
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) eta_input !eV
             if (eta_input < step_omega) then
                 print *, 'Attention :'
                 print *, '  broadening needs to be larger than step!!'
                 print *, ' broadening ::', eta_input, 'step::', step_omega
                 print *, ' This will cause wrong results or poor convergence.'
                 print *, ' we do not stop tough.'
                 print *, '________________________________________________________'
                 print *, ' '
             endif
                 
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) B_field  ! Teslas
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) gyro ! gyromagnetic factor
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) Temperature ! Kelvin
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) tol ! resolvent's convergence
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) Max_loop_index ! Maximum number of iterations for each consistence loop (there are two)
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) Bias_ini ! in eV
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) Bias_fin
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) N_bias
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) bias_fraction
           if (bias_fraction > 1.0 .or. bias_fraction <0) then
              print *, 'Attention :'
              print *, 'Bias fraction has to be a positive number below 1.0'
             print *, 'STOP.'
             stop
           endif
!
! input names of used files:
!
!  file containing one-electron impurity level (only one)
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) impurity_hamiltonian_file
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) gamma_file
!  output file with bias (V) and current in nA
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) output_file
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) N_vib ! Number of different modes
      allocate (Freq (N_vib))
      allocate (W (N_vib))
      do i = 1, N_vib
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) Freq (i) !  mode Frequency in meV
      Freq (i) = Freq (i)/(1000.*hartree) ! atomic units
      enddo
      do i = 1, N_vib
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) W (i) ! e-ph coupling for mode i (meV)
      W (i) = W (i) / (1000.*hartree) 
      enddo
      read (unit_input,*,IOSTAT=ios, ERR=100,END=200) vib_loop ! if true do self-consistence

      close (unit_input)

100  if (ios > 0) then
     print *,'ERROR reading WNnca.input, IOS=',ios
     stop
     endif
200  if (ios < 0) then
     print *,'END of WNnca.input, IOS=',ios
     stop
     endif

! symmetrized about zero due to finite arrays and omega --> - omega
! transformations:
!     Final omega = -omega_ini
! then the number of points in omega is:

       N_omega = 1 + 2 * int(-omega_ini / step_omega)
       allocate (omega (N_omega))
          
! Zeeman energy in atomic units :
       B_au =  gyro * B_field * 2.127159016E-6! Bohr magneton included

       if (B_field == 0) then
            N_spin=1
       else
            N_spin=2
       endif
! change all to atomic units :
       omega_ini = omega_ini / hartree
       step_omega =step_omega/hartree
       do i=1, N_omega
            omega (i) = omega_ini + (i-1)*step_omega
       enddo
       Temperature = Temperature * 25.852_q / (27211.6_q * 300._q)
       eta = ui * eta_input / hartree
       Bias_ini = Bias_ini / hartree
       Bias_fin = Bias_fin / hartree

       do i = 1, N_vib
       if (step_omega > Freq (i)) then
          print *, ' Problem: step_omega is larger than mode Freq!!'
          print *, ' We stop the calc, take a smaller step_omega.'
          stop
       endif
       enddo

       return
       end subroutine input

     subroutine memory (N_omega, N_spin, Gamma_t, Gamma_s, Fermi_t, Fermi_s,   &
      &   Dr, Gr, SelfD, selfG, Dl, Gl, Db, Gb, Hamiltonian, &
      &   SelfDl, SelfGl, &
      &   SelfDb, SelfGb, GGR, Gr_store, Dr_store, Dl_store, Gl_store)

     integer :: N_omega
     real (q), allocatable ::  Fermi_t (:), Fermi_s (:), omega (:)
     real (q), allocatable :: Gamma_t (:), Gamma_s (:)
     real (q), allocatable :: Hamiltonian (:)

     complex (qc), dimension (:), allocatable :: Dr, Dl, Db, Dr_store, Dl_store
     complex (qc), dimension (:,:), allocatable :: Gr, Gl, Gb, Gr_store, Gl_store
     complex (qc), dimension (:), allocatable :: SelfD, SelfDb, SelfDl
     complex (qc), dimension (:), allocatable :: SelfG
     complex (qc), dimension (:,:), allocatable :: SelfGb, SelfGl 
     complex (qc), dimension (:,:), allocatable :: GGR


     allocate (Dr (N_omega))
     allocate (Dr_store (N_omega))
     allocate (Gr (N_spin, N_omega))
     allocate (Gr_store (N_spin, N_omega))
     allocate (Dl (N_omega))
     allocate (Dl_store (N_omega))
     allocate (Gl (N_spin, N_omega))
     allocate (Gl_store (N_spin, N_omega))
     allocate (Db (N_omega))
     allocate (Gb (N_spin, N_omega))
     allocate (SelfD (N_omega))
     allocate (SelfG (N_omega))
     allocate (SelfDl (N_omega))
     allocate (SelfGl (N_spin, N_omega))
     allocate (SelfDb (N_omega))
     allocate (SelfGb (N_spin, N_omega))
     allocate (GGR (N_spin, N_omega))
     allocate (Hamiltonian (N_spin))
     allocate (Fermi_t (N_omega))
     allocate (Fermi_s (N_omega))
     allocate (Gamma_t (N_omega))
     allocate (Gamma_s (N_omega))

     return
     end subroutine memory


    end module io
