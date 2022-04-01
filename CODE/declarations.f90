module declarations

  implicit none

! PARAMETERS

  integer, parameter :: kk = SELECTED_INT_KIND (10)
  integer, parameter :: q = SELECTED_REAL_KIND(10)
  integer, parameter :: qs = SELECTED_REAL_KIND(5)
  integer, parameter :: qc = SELECTED_REAL_KIND(10)
  real (q), parameter :: pi_d = 3.14159265358979323846_q
  real (q), parameter :: sqrtpi_d = 1.77245385090551602729_q
  real (q), parameter :: hartree = 27.2116_q
  complex (qc), parameter :: zero = (0._q, 0._q)
  complex (qc), parameter :: ui = (0._q, 1._q)

! CHARACTERS

  character ( len = 100 ) :: impurity_hamiltonian_file
  character ( len = 100 ) :: gamma_file
  character ( len = 100 ) :: output_file

! LOGICALS

 logical :: converged, vib_loop

! integer loop indices

  integer :: i, j, k, i_spin, i_V
  integer :: dummy, ios
  integer :: self_index, Max_loop_index, i_omega
  integer :: N_vib

! array dimensions

  integer :: N_omega, N_spin, N_bias

! real variables

  real (q) :: omega_ini, omega_fin, step_omega,  tol, tol_a
  real (q) :: eta_input, gyro, B_field, B_au, E_fermi, Temperature
  real (q) :: Bias_ini, Bias_fin, bias, bias_fraction, curr

! complex variables

   complex (qc) :: eta, constant

! ALLOCATABLES
  integer, allocatable :: indexvib (:)

  real (q), allocatable ::  Fermi_t (:), Fermi_s (:), omega (:)
  real (q), allocatable :: Gamma_t (:), Gamma_s (:)
  real (q), allocatable :: Hamiltonian (:)
  real (q), allocatable, dimension (:) :: Freq, W, Bose

  complex (qc), dimension (:), allocatable :: Dr, Dl, Db, Dr_store, Dl_store
  complex (qc), dimension (:,:), allocatable :: Gr, Gl, Gb, Gr_store, Gl_store
  complex (qc), dimension (:), allocatable :: SelfD, SelfDb, SelfDl
  complex (qc), dimension (:), allocatable :: SelfG
 complex (qc), dimension (:,:), allocatable :: SelfGb, SelfGl
  complex (qc), dimension (:,:), allocatable :: GGR

end module declarations
