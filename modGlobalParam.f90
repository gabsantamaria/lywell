module modGlobalParam
    implicit none

	integer, parameter :: dp = kind(0.d0)
	integer, parameter :: dpo = kind(0.d0)
	integer, parameter :: dp_prec = 4!dp
	integer, parameter :: dp_near = 4!dp
	integer, parameter :: qp = 16
	integer, parameter :: il = 4
	integer, parameter :: ilo = 4

    real (kind = dp), parameter :: mu = 1.2566370614E-6_dp
    real (kind = dp), parameter ::  epsi = 8.8541878176E-12_dp
	real (kind = dp), parameter :: vp = 1./sqrt(mu*epsi)

    complex (kind = dp), parameter :: jj = cmplx (0.,1.)
    real (kind = dp), parameter :: pi = 3.141592653589793_dp
    real (kind = dp), parameter :: alphaCFIE = 1._dp!0.5

    real (kind = dp) :: lambda
    real (kind = dp) :: w
    real (kind = dp) :: kappa
	real (kind = dp) :: f

    type :: config
        logical :: msj
        integer (kind = il) :: Lforz
        integer (kind = il) :: numThetaRCS, numPhiRCS
        real (kind = dp) :: errorSolverCGS
        real (kind = dp) :: dirTheta, dirPhi
        real (kind = dp) :: rcsMonoFMin, rcsMonoFMax
        integer (kind = il) :: rcsMonoSamples
        character (len = 60) :: usar_precond, analisis_rcs
        logical :: prec_activo
        real (kind = dp) :: fill_in_prc, drop_tolerance
        integer (kind = il) :: n_pasos_rcs_mono
        character (len = 60) :: tipo_polarizac
        character (len = 10) :: variable_espec
        character (len = 30) file_name, file_script, test_name
        real (kind = dp) :: criterio_ladomax
    !	integer (kind = il) :: polVH !1 = Vertical, 2 = Horizontal
    end type
    type(config) :: conf


    contains

	function getlambda() result (r)
		!dummy
		real (kind = dp) :: r
		!
		r = lambda
	end function getlambda
	subroutine setlambda(v)
		!dummy
		real (kind = dp), intent(in) :: v
		!
		lambda = v
		w = 2.*pi*vp/lambda
		kappa = 2.*pi/lambda
		f = vp/lambda
	end subroutine


	function getw() result (r)
		!dummy
		real (kind = dp) :: r
		!
		r = w
	end function getw
	subroutine setw(v)
		!dummy
		real (kind = dp), intent(in) :: v
		!
		w = v
		lambda = 2.*pi*vp/w
		kappa = 2.*pi/lambda
		f = vp/lambda
	end subroutine



	function getkappa() result (r)
		!dummy
		real (kind = dp) :: r
		!
		r = kappa
	end function getkappa
	subroutine setkappa(v)
		!dummy
		real (kind = dp), intent(in) :: v
		!
		kappa = v
		lambda = 2.*pi/kappa
		w = 2.*pi*vp/lambda
		f = vp/lambda
	end subroutine


	function getf() result (r)
		!dummy
		real (kind = dp) :: r
		!
		r = f
	end function getf
	subroutine setf(v)
		!dummy
		real (kind = dp), intent(in) :: v
		!
		f = v
		lambda = vp/f
		w = 2.*pi*vp/lambda
		kappa = 2.*pi/lambda
	end subroutine

	subroutine msj(mensaje, id, dato)
		character (len = *), intent(in) :: mensaje
		integer (kind = il), intent(in), optional :: id
		real (kind = dp), intent(in), optional :: dato

		print*, mensaje
	end subroutine
end module modGlobalParam
