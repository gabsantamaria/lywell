module modRCS
    use modGlobalParam
    use modGlobalMetod
    implicit none

    private
    !Variables globales
    public :: inicializar_RCS, sigmaRCS
    !real ( kind = dp ), dimension (3) :: r_obs
    complex ( kind = dp ), pointer,  dimension (:) :: Vi
    real ( kind = dp ), pointer, dimension(:) :: e_long
    real ( kind = dp ), pointer, dimension(:,:) :: e_centro
    real ( kind = dp ), pointer, dimension (:,:) :: t_baric
    integer ( kind = il ), pointer, dimension(:,:) :: e_t
    integer ( kind = il ), pointer :: num_e
    integer ( kind = il ), pointer :: num_t

    complex (kind = dp), dimension(3) :: polariz
    complex (kind = dp) :: ctteOnda

contains

    subroutine inicializar_RCS (arg_centros_comun, arg_Vi, arg_long_l, arg_baric_t, arg_t_comun, arg_num_lcom, arg_num_t, arg_polariz, arg_ctteOnda)

        !dummy
        real ( kind = dp ),  dimension(:,:), target :: arg_centros_comun
        !real ( kind = dp ), dimension (3), intent(in) :: arg_r_obs
        complex ( kind = dp ),  dimension (:), intent(in), target :: arg_Vi
        real ( kind = dp ),  dimension(:), intent(in), target :: arg_long_l
        real ( kind = dp ),  dimension (:,:),intent(in), target :: arg_baric_t
        integer ( kind = il ),  dimension(:,:),intent(in), target :: arg_t_comun
        integer ( kind = il ), intent(in), target :: arg_num_lcom
        integer ( kind = il ), intent(in), target :: arg_num_t

        complex (kind = dp), dimension(3) :: arg_polariz
        complex (kind = dp) :: arg_ctteOnda

        !!!!!!!

        !r_obs = arg_r_obs
        num_e => arg_num_lcom
        num_t => arg_num_t
        !allocate (Vi(num_e))
        Vi => arg_Vi
        !allocate (e_centro(3, num_e))
        e_centro => arg_centros_comun
        !allocate (e_long(num_e))
        e_long => arg_long_l
        !allocate (t_baric(3,num_t))
        t_baric => arg_baric_t
        !allocate (e_t(2,num_e))
        e_t => arg_t_comun

        polariz = arg_polariz
        ctteOnda = arg_ctteOnda

    end subroutine


    function calc_E (r_obs) result (E)

        !dummy
        real ( kind = dp ), dimension (3), intent(in) :: r_obs
        complex (kind = dp), dimension(3) :: E
        !!!!!!!!
        complex (kind = dp), dimension(3) :: acumE
        complex (kind = dp), dimension(3) :: m
        complex (kind = dp), dimension(3) :: MM
        real (kind = dp), dimension(3) :: rm
        real (kind = dp), dimension(3) :: rrm
        integer (kind = il) :: i

        acumE = 0


        do i = 1,num_e

            rm = 0.5*(  t_baric(:, e_t(1,i)) + t_baric(:, e_t(2,i))   )
            m  = (  t_baric(:, e_t(2,i)) - t_baric(:, e_t(1,i))   )*Vi(i)*e_long(i)
            rrm = r_obs - rm
            !MM = producto_escalar(m,rrm)*rrm/(modulov(rrm)**2)
            MM = sum(m*rrm)*rrm/(dot_product(rrm,rrm))
            acumE=acumE + ( green(r_obs,rm)*(MM(:)-m(:)) )
        !print*, acumE
        end do

        E = acumE*jj*w*mu/(4*pi)
        !print*, E
    end function calc_E

    function sigmaRCS(theta, phi) result (sigma)
        !dummy
        real (kind = dp), intent(in) :: theta, phi
        real (kind = dp) :: sigma
        !
        !local
        real (kind = dp) :: r_zl
        complex (kind = dp), dimension(3) :: E
        E = calc_RCS(theta, phi, r_zl)

        sigma = dot_product(E,E)

        sigma = sigma*4*pi*(r_zl**2)/(dot_product(polariz*ctteOnda,polariz*ctteOnda))

    end function sigmaRCS



    function calc_RCS(theta, phi, r_zl) result(ERCS)
        !dummy
        real (kind = dp), intent(in) :: theta, phi
        complex (kind = dp), dimension(3) :: ERCS
        real (kind = dp), intent(out) :: r_zl
        !!!!!
        !local
        real (kind = dp), dimension (3) :: roo
        integer (kind = il), dimension(6) :: limits
        real (kind = dp) :: D, maxX, minX, maxY, minY, maxZ, minZ, rzl, lambda
        real (kind = dp), dimension(3) :: r_obs
        !

        limits = limites_dispersor(e_centro, num_e)

        maxX = e_centro(1,limits(1))
        minX = e_centro(1,limits(2))
        maxY = e_centro(2,limits(3))
        minY = e_centro(2,limits(4))
        maxZ = e_centro(3,limits(5))
        minZ = e_centro(3,limits(6))
        lambda = 2*pi/(w*sqrt(mu*epsi))
        D = sqrt(    ((maxX-minX)**2) + ((maxY-minY)**2) + ((maxZ-minZ)**2)   )

        rzl = 2*(D**2)/lambda
        if (50*D>rzl) then
            rzl = 50*D
        end if
        if (20*lambda>rzl) then
            rzl = 20*lambda
        end if
        rzl = rzl*10 !multiplicador de zona lejana
        r_zl = rzl
        roo(1) = 0.5*(maxX+minX)
        roo(2) = 0.5*(maxY+minY)
        roo(3) = 0.5*(maxZ+minZ)

        r_obs(1) = roo(1)+ ( rzl*sin(theta)*cos(phi) )
        r_obs(2) = roo(2)+ ( rzl*sin(theta)*sin(phi) )
        r_obs(3) = roo(3)+ ( rzl*cos(theta) )

        !return
        ERCS = calc_E(r_obs)

    end function calc_RCS







end module modRCS
