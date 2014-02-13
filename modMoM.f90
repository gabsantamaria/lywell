module modMoM
    use modGlobalParam
    use modGlobalMetod
    implicit none
    private
    public :: inicializar_MoM, llenarZ_MoM, Zmn_MoM, MatxVec_MoM, setZ_MatxVec_MoM
    complex (kind = dp), pointer, dimension (:,:) :: ZMoM_MxV !puntero copia de Z para la funcion de MatxVec_MoM
    real ( kind = dp ),  pointer, dimension (:,:) :: p_coord
    real ( kind = dp ), pointer, dimension (:,:) :: t_baric
    real ( kind = dp ), pointer, dimension (:,:,:) :: t_baric_sub
    real ( kind = dp ), pointer, dimension(:) :: t_area, e_long
    integer ( kind = il ), pointer, dimension (:,:) :: t_p, e_p, e_t, e_po
    integer ( kind = il ), pointer :: num_e, num_p, num_t

contains



    subroutine inicializar_MoM (arg_p_xyz, arg_baric_t, arg_baric_sub_t, arg_area_t, arg_long_l, arg_triang_p, arg_l_comun, arg_t_comun, arg_vo, arg_num_lcom, arg_num_p, arg_num_t)
        !dummy

        real ( kind = dp ),  dimension (:,:), intent(in), target :: arg_p_xyz
        real ( kind = dp ),  dimension (:,:), intent(in), target :: arg_baric_t
        real ( kind = dp ),  dimension (:,:,:), intent(in), target :: arg_baric_sub_t
        real ( kind = dp ),  dimension(:), intent(in), target :: arg_area_t, arg_long_l
        integer ( kind = il ),  dimension (:,:), intent(in), target :: arg_triang_p, arg_l_comun, arg_t_comun, arg_vo
        integer ( kind = il ), intent(in), target :: arg_num_lcom, arg_num_p, arg_num_t
        !real ( kind = dp ), intent(in) :: arg_w
print*, 'se ha inicializado MoM'
        !allocate(p_coord(3,arg_num_p))
        p_coord => arg_p_xyz
        !allocate(t_baric(3,arg_num_t))
        t_baric => arg_baric_t
        !allocate(t_baric_sub(3,9,arg_num_t))
        t_baric_sub => arg_baric_sub_t
        !allocate(t_area(arg_num_t))
        t_area => arg_area_t
        !allocate(e_long(arg_num_lcom))
        e_long => arg_long_l
        !allocate(t_p(3,arg_num_t))
        t_p => arg_triang_p
        !allocate(e_p(2,arg_num_lcom))
        e_p => arg_l_comun
        !allocate(e_t(2,arg_num_lcom))
        e_t => arg_t_comun
        !allocate(e_po(2,arg_num_lcom))
        e_po => arg_vo

        num_e => arg_num_lcom
        num_p => arg_num_p
        num_t => arg_num_t
        !w = arg_w

    end subroutine


    subroutine llenarZ_MoM(Z, mensaje)
        !dummy
        complex (kind = dp), dimension(num_e,num_e), intent(out) :: Z
        integer (kind = il), optional, intent (in) :: mensaje
        !local
        integer (kind = il) :: m, n
        !integer (kind = il) :: prevPrc, Prc
        !prevPrc=-1
        call setprc(num_e)
        do m = 1, num_e  !punto pbservación
            !    print*, m
            do n = 1, num_e  !punto fuente
                Z(m,n) = Zmn_MoM(m,n)
            end do

            if (present(mensaje)) then
                call updateprc(m)
                !Prc = (100*m/num_e)
                !if (prevPrc < Prc) then
                !    !Write(strprint, '(I5)') Prc
                !    print*, inttostr(Prc) // ' %'
                !    !call imprPRC(Prc)

                !    prevPrc = Prc
                !end if
            end if

        end do
		call setZ_MatxVec_MoM(Z)
    end subroutine llenarZ_MoM


    function Zmn_MoM(m,n) result (valorz)
        !dummy
        integer (kind = il), intent (in) :: m, n
        complex (kind = dp) :: valorz
        !local
        integer (kind = il) :: i
        complex (kind = dp) :: PhiP, PhiM, PintA, PintPhi
        complex (kind = dp), dimension(3) :: AmnP, AmnM
        complex (kind = dp) :: g11, g12, g21, g22 !mn: n=fuente y m=observacion del argumento de la funcion de green
        real (kind = dp), dimension(3) :: rho1, rho2
        AmnP = 0
        AmnM = 0
        PhiP = 0
        PhiM = 0

        do i=1,9
            g11 = green(t_baric(:,e_t(1,m)),t_baric_sub(:,i,e_t(1,n)))
            g12 = green(t_baric(:,e_t(1,m)),t_baric_sub(:,i,e_t(2,n)))
            g21 = green(t_baric(:,e_t(2,m)),t_baric_sub(:,i,e_t(1,n)))
            g22 = green(t_baric(:,e_t(2,m)),t_baric_sub(:,i,e_t(2,n)))
            !print*, g11
            !print*, g12
            !print*, g21
            !print*, g22
            rho1 = t_baric_sub(:,i,e_t(1,n))-p_coord(:,e_po(1,n))
            rho2 = p_coord(:,e_po(2,n))-t_baric_sub(:,i,e_t(2,n))
            !
            AmnP = AmnP + (g11*rho1) + (g12*rho2)
            AmnM = AmnM + (g21*rho1) + (g22*rho2)
            !
            PhiP = PhiP + (g11-g12)
            PhiM = PhiM + (g21-g22)
        end do
        !Constantes multiplicativas faltantes
        AmnP = AmnP*(mu/pi)*(1./72.)*e_long(n)
        AmnM = AmnM*(mu/pi)*(1./72.)*e_long(n)
        !
        PhiP = (-1.)*PhiP*e_long(n)/(jj*36*pi*w*epsi)
        PhiM = (-1.)*PhiM*e_long(n)/(jj*36*pi*w*epsi)

        !Productos internos
        !PintA = (jj*0.5*w*e_long(m))*(producto_escalar(AmnP, t_baric(:,e_t(1,m))- p_coord(:,e_po(1,m))) + producto_escalar(AmnM, -t_baric(:,e_t(2,m))+p_coord(:,e_po(2,m))))
        PintA = (jj*0.5*w*e_long(m))*(sum(AmnP*( t_baric(:,e_t(1,m))- p_coord(:,e_po(1,m)) )  ) + sum(AmnM*( -t_baric(:,e_t(2,m))+p_coord(:,e_po(2,m)) )  ))
        PintPhi = e_long(m)*(PhiM-PhiP)

        valorz = PintPhi+PintA
    end function Zmn_MoM

	subroutine setZ_MatxVec_MoM(matZ)
		complex ( kind = dp ), intent(in), target, dimension(:,:) :: matZ

        ZMoM_MxV => matZ
	end subroutine

    function MatxVec_MoM(vecX) result(res)
    	complex (kind = dp), intent(in), dimension(:) :: vecX
		complex (kind = dp), dimension(size(vecX,1)) :: res
		res = matmul(ZMoM_MxV,vecX)
    end function MatxVec_MoM

end module modMoM
