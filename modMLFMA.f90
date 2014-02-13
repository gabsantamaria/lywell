module modMLFMA
    use modGlobalParam
    use modGlobalMetod
    use modMoM
    use modIterativo
    implicit none



    private
    !Variables globales
    public :: inicializar_MLFMA, MatxVec_MLFMA, Zesparcida_MLFMA, ia_MLFMA, ja_MLFMA, n_esparcida_MLFMA, cerarZ, destruir_MLFMA

    !Variables puntero de geometria del dispersor
    real ( kind = dp ),  pointer, dimension (:,:) :: p_coord
    integer ( kind = il ), pointer, dimension (:,:) :: t_p, e_p
    integer ( kind = il ), pointer :: num_p, num_t
    real ( kind = dp ), pointer, dimension(:,:) :: e_centro
    integer (kind = il), pointer :: num_e
    real ( kind = dp ), pointer, dimension (:,:) :: t_baric
    real ( kind = dp ), pointer, dimension (:,:,:) :: t_baric_sub
    real ( kind = dp ), pointer, dimension(:) :: t_area, e_long
    real ( kind= dp ), pointer, dimension(:,:) :: t_normal
    integer ( kind = il ), pointer, dimension (:,:) :: e_t, e_po
    !
    !parametros de entrada del metodo
    integer (kind = il), parameter :: dig_multipolos = 6!10
    real (kind = dp) :: criterio_ladoW
    !
    integer (kind = il) :: n_esparcida_MLFMA
    integer (kind = il), allocatable, dimension(:) :: ia_MLFMA,ja_MLFMA
    complex (kind = dp_near), allocatable, dimension(:) :: Zesparcida_MLFMA
    complex (kind = dp), allocatable, dimension(:,:,:) :: funcRAD, funcREC
    integer (kind = il) :: nivel_max

    type :: cubo
        !propiedades intrinsecas
        real (kind = dp), dimension (3) :: centro
        integer (kind = il) :: padre
        !
        !funciones bases dentro del cubo
        integer (kind = il), allocatable, dimension(:) :: bases_propias
        integer (kind = il) :: n_propias
        !
        !funciones bases cercanas al cubo (incluyendo propias)
        integer (kind = il), allocatable, dimension(:) :: bases_cercanas
        integer (kind = il) :: n_cercanas
        !
        !Matriz Zcercana
        complex (kind = dp), allocatable, dimension(:,:) :: Zcerc
        !
        !cubos lejanos e indices de funcion de transferencia asociada
        integer (kind = il), allocatable, dimension(:,:) :: cubLej_indFuncTransf
        integer (kind = il) n_cubos_lejanos
        !
    end type



    type :: leveloctree

        type(cubo), allocatable, dimension(:) :: Cubos
        integer (kind = il) :: tamano_aux !tamaño actual del vector Cubos
        real (kind = dp) :: lado_cubo
        integer (kind = il) :: n_cubos
        integer (kind = il) :: n_padres
        integer (kind = il) :: n_muestras_theta, n_muestras_phi
        integer (kind = il) :: n_multipolos
        integer (kind = il) :: n_kappas
        complex (kind = dp), allocatable, dimension(:,:) :: functransf
        integer (kind = il) :: n_functransf
        real (kind = dp), allocatable, dimension(:,:) :: kappas
        real (kind = dp), allocatable, dimension(:,:) :: atheta, aphi
        real (kind = dp), allocatable, dimension(:,:) :: ang_kappas
        real (kind = dp), allocatable, dimension(:,:) :: mat_interp, mat_anterp
        real (kind = dp), allocatable, dimension(:) :: pesos_integ
        complex (kind = dp), allocatable, dimension(:,:,:) :: campoC, campoB, campoA

    end type

    type(leveloctree), allocatable, dimension(:) :: octree
integer :: num_iter_current = 0 !#prueba
contains

    subroutine destruir_MLFMA()
        if (allocated(octree)) then
            deallocate(octree)
        end if

        if (allocated(ia_MLFMA)) then
            deallocate(ia_MLFMA)
        end if

        if (allocated(ja_MLFMA)) then
            deallocate(ja_MLFMA)
        end if

        if (allocated(Zesparcida_MLFMA)) then
            deallocate(Zesparcida_MLFMA)
        end if

        if (allocated(funcRAD)) then
            deallocate(funcRAD)
        end if

        if (allocated(funcREC)) then
            deallocate(funcREC)
        end if

    end subroutine


    subroutine inicializar_MLFMA (arg_e_centro, arg_num_e, arg_e_t, arg_e_po, arg_e_long, arg_t_area, arg_t_baric, arg_t_baric_sub, arg_p_coord, arg_t_p, arg_e_p, arg_num_p, arg_num_t, arg_t_normal)
        !dummy
        real ( kind = dp ),  dimension(:,:), target :: arg_e_centro
        integer (kind = il), target :: arg_num_e
        real ( kind = dp ),  dimension (:,:), intent(in), target :: arg_t_baric
        real ( kind = dp ),  dimension (:,:,:), intent(in), target :: arg_t_baric_sub
        real ( kind = dp ),  dimension(:), intent(in), target :: arg_t_area, arg_e_long
        integer ( kind = il ),  dimension (:,:), intent(in), target :: arg_e_t, arg_e_po
        real ( kind = dp ),  target, dimension (:,:) :: arg_p_coord
        integer ( kind = il ), target, dimension (:,:) :: arg_t_p, arg_e_p
        integer ( kind = il ), target :: arg_num_p, arg_num_t
        real ( kind= dp ), target, dimension(:,:) :: arg_t_normal
        !
        !local
        integer (kind = il), allocatable, dimension(:) :: nMultipolos
        integer (kind = il), allocatable, dimension(:) :: BusquedaInicial
        integer (kind = il) ::  i
        integer (kind = il), dimension(6) :: limites
        real (kind = dp) :: ladomayor, deltaX, deltaY, deltaZ, maxX, minX, maxY, minY, maxZ, minZ, d
        real (kind = dp), dimension (3) :: roo
        real (kind = dp) :: tol
        integer (kind = il) :: criterio_incremento
                !

        criterio_incremento = 100
        criterio_ladoW = conf%criterio_ladomax!0.5_dp

        num_e => arg_num_e
        e_centro => arg_e_centro
        e_t => arg_e_t
        e_po => arg_e_po
        e_long => arg_e_long
        t_area => arg_t_area
        t_baric => arg_t_baric
        t_baric_sub => arg_t_baric_sub
		t_normal => arg_t_normal

        p_coord => arg_p_coord
        t_p => arg_t_p
        e_p => arg_e_p
        num_p => arg_num_p
        num_t => arg_num_t
        ! Inicializar MoM
        call inicializar_MoM(p_coord, t_baric, t_baric_sub, t_area, e_long, t_p, e_p, e_t, e_po, num_e, num_p, num_t)
        !

        limites = limites_dispersor(e_centro, num_e)
        maxX = e_centro(1,limites(1))
        minX = e_centro(1,limites(2))
        maxY = e_centro(2,limites(3))
        minY = e_centro(2,limites(4))
        maxZ = e_centro(3,limites(5))
        minZ = e_centro(3,limites(6))
        deltaX = maxX-minX
        deltaY = maxY-minY
        deltaZ = maxZ-minZ
        roo(1) = 0.5*(maxX+minX)
        roo(2) = 0.5*(maxY+minY)
        roo(3) = 0.5*(maxZ+minZ)


        if (  (deltaX >= deltaY) .and. (deltaX >= deltaZ) ) then
            ladomayor = deltaX
        else if (  (deltaY >= deltaX) .and. (deltaY >= deltaZ) ) then
            ladomayor = deltaY
        else
            ladomayor = deltaZ
        end if
        !Tolerancia para las funciones bases en los bordes
        tol = getlambda()/10000.
        ladomayor = ladomayor + tol
        tol = tol/2.
        print*, 'El cubo que contiene al dispersor tiene lado igual a ' // realtostr (ladomayor) // ' metros'
        print*, 'Asignando espacio en memoria y estableciendo muestras sobre la esfera unitaria...'
        !Busqueda del nivel maximo
        !nivel_max = ceiling(2. + ( log10(ladomayor/getlambda())/log10(2.) ))
        nivel_max = ceiling(1. + (    log10( ladomayor/(criterio_ladoW*getlambda()) ) / log10(2.)    ))



        print*, 'El dispersor se dividira en ' // inttostr(nivel_max) // ' niveles'

        if (nivel_max <= 3) then
            print*, 'El dispersor es demasiado pequeño eléctricamente para aplicar MLFMA'
            stop "El programa ha finalizado"
        end if

        allocate(octree(nivel_max))
        do i = 1, nivel_max
            octree(i)%lado_cubo = ladomayor/(2.**(i-1))
            d = sqrt(3.)*octree(i)%lado_cubo
            octree(i)%n_multipolos = floor(((2.*pi/getlambda())*d) + dig_multipolos*(((2.*pi/getlambda())*d)**(1./3.)))
            octree(i)%n_muestras_theta = octree(i)%n_multipolos
            octree(i)%n_muestras_phi = 2.*octree(i)%n_multipolos
            octree(i)%n_kappas = octree(i)%n_muestras_theta*octree(i)%n_muestras_phi
            if (i >= 3) then
                allocate(octree(i)%pesos_integ(octree(i)%n_kappas))
                allocate(octree(i)%kappas(3,octree(i)%n_kappas))
                allocate(octree(i)%atheta(3,octree(i)%n_kappas))
                allocate(octree(i)%aphi(3,octree(i)%n_kappas))
                allocate(octree(i)%ang_kappas(2,octree(i)%n_kappas))
            end if
            octree(i)%n_cubos = 0
            octree(i)%tamano_aux = 0
            octree(i)%n_padres = 0
        end do
        !
        do i = 1, nivel_max
            if (i >= 3) then
                allocate(octree(i)%Cubos(criterio_incremento))
                octree(i)%tamano_aux = criterio_incremento
            else
                allocate(octree(i)%Cubos(8**(i-1)))
                octree(i)%tamano_aux = 8**(i-1)
            end if
        end do

        octree(1)%n_cubos = 1
        octree(1)%Cubos(1)%centro = roo

        allocate(BusquedaInicial(num_e))
        do i = 1, num_e
            !La primera vez todas las funciones bases estaran dentro del cubo
            Busquedainicial(i) = i
        end do

        allocate(funcRAD(octree(nivel_max)%n_kappas,2,num_e))
        allocate(funcREC(octree(nivel_max)%n_kappas,2,num_e))

        call llenar_kappas()
        print*, 'Hecho'
        print*, 'Creando geometria OCTREE y calculando funciones de radiacion...'
        call geomMLFMA (roo, ladomayor, BusquedaInicial, num_e,  1, nivel_max,criterio_incremento)
		!funcREC = conjg(funcRAD)
        print*, 'Hecho.'

        do i = 3, nivel_max

            print*, 'Nivel ' // inttostr(i) // ' tiene ' // inttostr(octree(i)%n_cubos) // ' cubos, ' // inttostr(octree(i)%n_multipolos) // 'multipolos y ' // inttostr(octree(i)%n_kappas) // ' muestras sobre la esfera unitaria'

            allocate(octree(i)%campoC(octree(i)%n_kappas,2,octree(i)%n_cubos))
            allocate(octree(i)%campoB(octree(i)%n_kappas,2,octree(i)%n_cubos))
            allocate(octree(i)%campoA(octree(i)%n_kappas,2,octree(i)%n_cubos))

        end do

        print*, 'Llenando funciones de transferencia...'
        call llenar_functransf (nivel_max)
        print*,'Hecho.'

        print*, 'Llenando Z cercana...'
        call llenarZesparcida()
        !call llenarZesparcida()
        print*, 'Hecho'

        print*, 'Llenando matrices de interpolacion'
        call llenarMinterp(nivel_max)
        print*, 'Hecho'

	call monitor_memoria()

    
    end subroutine




    function cerarZ(Zin) result(Zout)
    !dummy

    complex (kind = dp), dimension(:,:), intent(in) :: Zin
    complex (kind = dp), allocatable, dimension(:,:) :: Zout
    !
    !local
    integer (kind = il) :: cub, cubcer, indCentro, m, n, indm, indn, cont, f, nf, indj, indice
    !
    allocate(Zout(num_e,num_e))
    Zout(:,:) = 0.!Zin(:,:)
cont = 0

        do cub = 1, octree(nivel_max)%n_cubos
            nf = octree(nivel_max)%Cubos(cub)%n_propias
            do f = 1, nf
                indice = octree(nivel_max)%Cubos(cub)%bases_propias(f)

                do m = 1, octree(nivel_max)%Cubos(cub)%n_cercanas
                    indj = octree(nivel_max)%Cubos(cub)%bases_cercanas(m)
                    cont = cont + 1
                    Zout(indice,indj) = 1.
                end do

            end do
        end do

    print*, 'NO CERADOS********************************** ', cont
    Zout(:,:) = Zout(:,:)*Zin(:,:)
    end function cerarZ



    function cerarZ_old(Zin) result(Zout)
    !dummy

    complex (kind = dp), dimension(:,:), intent(in) :: Zin
    complex (kind = dp), allocatable, dimension(:,:) :: Zout
    !
    !local
    integer (kind = il) :: cub, cublej, indCentro, m, n, indm, indn, cont
    !
    allocate(Zout(num_e,num_e))
    Zout(:,:) = Zin(:,:)
cont = 0
    do cub = 1,octree(nivel_max)%n_cubos
        do cublej = 1,octree(nivel_max)%Cubos(cub)%n_cubos_lejanos
            indCentro = octree(nivel_max)%Cubos(cub)%cubLej_indFuncTransf(1,cublej)
            do m = 1, octree(nivel_max)%Cubos(cub)%n_propias
                indm = octree(nivel_max)%Cubos(cub)%bases_propias(m)
                do n = 1, octree(nivel_max)%Cubos(indCentro)%n_propias                  
                    indn = octree(nivel_max)%Cubos(indCentro)%bases_propias(n)
                    Zout(indm, indn) = 0.
                    cont = cont + 1
                end do
            end do
        end do
    end do

    print*, 'CERADOS********************************** ', cont
    end function cerarZ_old


	subroutine  monitor_memoria()
	!local
	integer (kind = il) :: mem, memt, i, c, memp
	!
	memt = 0
	memp = 0
	print*, '****** Memoria ocupada por el simulador *******'
	memt = memt + sizeof(p_coord) + sizeof(t_p) + sizeof(e_p) + sizeof(e_centro) + sizeof(t_baric) + sizeof(t_baric_sub) + sizeof(t_area) + sizeof(e_long) + sizeof(t_normal) + sizeof(e_t) + sizeof(e_po)
	print*, 'Memoria de geometria: ', memt, 'bytes = ' , memt/1000000., ' MB'

	mem = sizeof(funcRAD)
	memt = memt + mem
	print*, 'funcRAD:  ', mem, 'bytes = ' , mem/1000000., ' MB'

	mem = sizeof(funcREC)
	memt = memt + mem
	print*, 'funcREC: ', mem, 'bytes = ' , mem/1000000., ' MB'

	mem = sizeof(octree)
	memt = memt + mem
	print*, '*****octree: ', mem, 'bytes = ' , mem/1000000., ' MB'
	do i = 1, nivel_max
		print*, 'Memoria contabilizada hasta ahora: ', memt, 'bytes = ' , memt/1000000., ' MB'
		print*, '	octree, nivel ' // inttostr(i)
		mem = sizeof(octree(i)%functransf)
		memt = memt + mem
		!print*, '		functransf: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%kappas)
		memt = memt + mem
		!print*, '		kappas: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%atheta)
		memt = memt + mem
		!print*, '		atheta: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%aphi)
		memt = memt + mem
		!print*, '		aphi: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%ang_kappas)
		memt = memt + mem
		!print*, '		ang_kappas: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%mat_interp)
		memt = memt + mem
		!print*, '		mat_interp: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%mat_anterp)
		memt = memt + mem
		!print*, '		mat_anterp: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%pesos_integ)
		memt = memt + mem
		!print*, '		pesos_integ: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%campoC)
		memt = memt + mem
		!print*, '		campoC: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%campoB)
		memt = memt + mem
		!print*, '		campoB: ', mem, 'bytes = ' , mem/1000000., ' MB'
		mem = sizeof(octree(i)%campoA)
		memt = memt + mem
		!print*, '		campoA: ', mem, 'bytes = ' , mem/1000000., ' MB'

		!print*, '		*****Memoria acumulada por cubos...'

		memp = 0
		do c = 1, octree(i)%n_cubos
			mem = sizeof(octree(i)%Cubos(c)%bases_propias)
			memp = memp + mem
		end do
		!print*, '			bases_propias: ', memp, 'bytes = ' , memp/1000000., ' MB'
		memt = memt + memp

		memp = 0
		do c = 1, octree(i)%n_cubos
			mem = sizeof(octree(i)%Cubos(c)%bases_cercanas)
			memp = memp + mem
		end do
		!print*, '			bases_cercanas: ', memp, 'bytes = ' , memp/1000000., ' MB'
		memt = memt + memp

		memp = 0
		do c = 1, octree(i)%n_cubos
			mem = sizeof(octree(i)%Cubos(c)%Zcerc)
			memp = memp + mem
		end do
		!print*, '			Zcerc: ', memp, 'bytes = ' , memp/1000000., ' MB'
		memt = memt + memp

		memp = 0
		do c = 1, octree(i)%n_cubos
			mem = sizeof(octree(i)%Cubos(c)%cubLej_indFuncTransf)
			memp = memp + mem
		end do
		!print*, '			cubLej_indFuncTransf: ', memp,  'bytes = ' , memp/1000000., ' MB'
		memt = memt + memp


	end do
        mem = sizeof(Zesparcida_MLFMA)
        memt = memt + mem
        print*, 'Zesparcida_MLFMA: ', mem, 'bytes = ' , mem/1000000., ' MB'

        mem = sizeof(ia_MLFMA)
        memt = memt + mem
        print*, 'ia_MLFMA (Zesparcida_MLFMA): ', mem, 'bytes = ' , mem/1000000., ' MB'

        mem = sizeof(ja_MLFMA)
        memt = memt + mem
        print*, 'ja_MLFMA (Zesparcida_MLFMA): ', mem, 'bytes = ' , mem/1000000., ' MB'
	print*, ' MEMORIA TOTAL: ', memt ,  'bytes = ' , memt/1000000., ' MB'

	end subroutine


    subroutine llenar_kappas()
        !local
        integer (kind = il) :: i,j,k, cont
        real (kind = dp) :: deltatheta, deltaphi, theta, phi
        !

        do i = 3, nivel_max

            cont = 0
            deltaphi = 2.*pi/octree(i)%n_muestras_phi
            deltatheta = pi/octree(i)%n_muestras_theta

            do j = 1, octree(i)%n_muestras_phi

                phi = (deltaphi/2.)*((2*j)-1)

                do k = 1,octree(i)%n_muestras_theta

                    theta = (deltatheta/2.)*((2*k)-1)

                    cont = cont + 1

                    octree(i)%kappas(1,cont) = sin(theta)*cos(phi)
                    octree(i)%kappas(2,cont) = sin(theta)*sin(phi)
                    octree(i)%kappas(3,cont) = cos(theta)

                    octree(i)%atheta(1,cont) = cos(theta)*cos(phi)
                    octree(i)%atheta(2,cont) = cos(theta)*sin(phi)
                    octree(i)%atheta(3,cont) = -sin(theta)

                    octree(i)%aphi(1,cont) = -sin(phi)
                    octree(i)%aphi(2,cont) = cos(phi)
                    octree(i)%aphi(3,cont) = 0.

                    octree(i)%ang_kappas(1,cont) = theta
                    octree(i)%ang_kappas(2,cont) = phi
                    octree(i)%pesos_integ(cont) = (cos((k-1)*deltatheta)-cos(k*deltatheta))*deltaphi
                end do
            end do

        end do

    end subroutine

    recursive subroutine geomMLFMA (centro, lcubo, vBusqueda, n_busq,  nivel, niveles, criterio_incremento)
        !dumy
        real (kind = dp), dimension(3), intent(in) :: centro
        real (kind = dp), intent(in) :: lcubo
        integer (kind = il), dimension(:), intent(in) :: vBusqueda
        integer (kind = il), intent(in) :: nivel, niveles
        integer (kind = il), intent(in) :: n_busq
        integer (kind = il), intent(in) :: criterio_incremento
        !
        !local
        integer (kind = il), allocatable, dimension(:) :: vEncontrados, vRestantes
        integer (kind = il) :: n_encon, i, n_rest

        logical :: encontro !fin
        real (kind = dp), dimension(3,8) :: centro_hijo
        real (kind = dp) :: lcubo_hijo
        integer (kind = il), allocatable, dimension(:) :: vBusquedaNew
        integer (kind = il) :: n_busqNew
        !integer (kind = il) :: prc
        type(cubo), allocatable, dimension(:) :: auxCubos

        !
        if (nivel == 1) then
            call setprc(8)
        end if

        if ( nivel == niveles ) then

            call calcRad (centro, vBusqueda, n_busq, niveles, .false.)

        else
            call dividir_en_8_hijos(centro, lcubo, centro_hijo, lcubo_hijo)

            allocate(vBusquedaNew(n_busq))
            vBusquedaNew = vBusqueda
            n_busqNew = n_busq

            octree(nivel+1)%n_padres = octree(nivel+1)%n_padres + 1

            do i = 1,8

                !Buscar funciones base
                call buscarFuncBases(centro_hijo(:,i), lcubo_hijo, vBusquedaNew,n_busqNew, vEncontrados,n_encon,vRestantes,n_rest, encontro)

                if (encontro) then
                    !guardar en octree

                    !si hay que redimensionar, hacerlo
                    if (octree(nivel+1)%n_cubos == octree(nivel+1)%tamano_aux) then

                        allocate(auxCubos(octree(nivel+1)%n_cubos))
                        auxCubos = octree(nivel+1)%Cubos
                        deallocate(octree(nivel+1)%Cubos)
                        allocate(octree(nivel+1)%Cubos(octree(nivel+1)%n_cubos + criterio_incremento))
                        octree(nivel+1)%Cubos(1:octree(nivel+1)%n_cubos) = auxCubos(:)
                        octree(nivel+1)%tamano_aux = octree(nivel+1)%n_cubos + criterio_incremento
                        deallocate(auxCubos)

                    end if
                    !
                    octree(nivel+1)%n_cubos = octree(nivel+1)%n_cubos + 1
                    octree(nivel+1)%Cubos(octree(nivel+1)%n_cubos)%centro= centro_hijo(:,i)
                    octree(nivel+1)%Cubos(octree(nivel+1)%n_cubos)%padre = octree(nivel+1)%n_padres

                    allocate(octree(nivel+1)%Cubos(octree(nivel+1)%n_cubos)%bases_propias(n_encon))
                    octree(nivel+1)%Cubos(octree(nivel+1)%n_cubos)%bases_propias = vEncontrados
                    octree(nivel+1)%Cubos(octree(nivel+1)%n_cubos)%n_propias = n_encon

                    call geomMLFMA (centro_hijo(:,i),lcubo_hijo,vEncontrados,n_encon,nivel+1,niveles, criterio_incremento)

                    deallocate(vBusquedaNew)
                    allocate(vBusquedaNew(n_rest))
                    vBusquedaNew = vRestantes
                    n_busqNew = n_rest
                end if
                call updateprc(i)
                !if (nivel == 1) then
                !    prc = ((i/8.)*100.)
                !    print*, inttostr(prc) // ' %'
                !end if
            end do
            deallocate(vBusquedaNew)
        end if

    end subroutine


    subroutine calcRad (centro, funcBases, nBases,niveles, maxPrecision)

        !dummy
        real (kind = dp), dimension(3), intent(in) :: centro
        integer (kind = il), dimension(nBases), intent(in) :: funcBases
        integer (kind = il), intent(in) :: nBases , niveles
        logical, optional :: maxPrecision
        !
        !local
        integer (kind = il) :: ii,j,i
        real (kind = dp), dimension(3) :: baricentroplus, vopuestoplus, rhoplus, baricentrominus, vopuestominus, rhominus, vkappa
        real (kind = dp) :: areaplus,areaminus,longitud, ka
        complex (kind = dp), dimension(3) :: vRadplus, vRadminus, vRadtotal, vRec
        real (kind = dp), dimension(3) :: vnormalPlus, vnormalMinus
        complex (kind = dp) :: UNO
        !
		UNO = (jj)/(jj)
        ka = (2.*pi)/getlambda()

        if (present(maxPrecision)) then
        else
            maxPrecision = .false.
        end if

        if (maxPrecision) then

            do i = 1,nBases
                !triangulo mas

                vopuestoplus = p_coord(:,e_po(1,funcBases(i)))
                areaplus = t_area(e_t(1, funcBases(i)))
                longitud = e_long(funcBases(i))

                !triangulo menos

                vopuestominus = p_coord(:,e_po(2,funcBases(i)))
                areaminus = t_area(e_t(2, funcBases(i)))

                !vector radiacion
                do j=1,octree(niveles)%n_kappas
                    vkappa = octree(niveles)%kappas(:,j)
                    vRadplus=0.
                    vRadminus=0.
                    do ii = 1, 9
                        baricentroplus = t_baric_sub(:, ii, e_t(1, funcBases(i)))
                        rhoplus = baricentroplus - vopuestoplus
                        baricentrominus = t_baric_sub(:, ii, e_t(2, funcBases(i)))
                        rhominus = vopuestominus - baricentrominus


                        vRadplus = vRadplus + exp(jj*ka*dot_product(vkappa,baricentroplus-centro))*rhoplus
                        vRadminus = vRadminus + exp(jj*ka*dot_product(vkappa,baricentrominus-centro))*rhominus
                    end do

					vnormalPlus = t_normal(:,e_t(1, funcBases(i)))
					vnormalPlus = vnormalPlus/norm2(vnormalPlus)
					vnormalMinus = t_normal(:,e_t(2, funcBases(i)))
					vnormalMinus = vnormalMinus/norm2(vnormalMinus)


                    vRadtotal = (vRadplus + vRadminus)*longitud*(1./18.)

                    funcRAD(j,1,funcBases(i)) = sum(vRadtotal*octree(niveles)%atheta(:,j))
                    funcRAD(j,2,funcBases(i)) = sum(vRadtotal*octree(niveles)%aphi(:,j))

					vRec = cross_productC ( -vkappa*UNO , cross_productC( conjg(vRadplus*longitud/18.) , vnormalPlus*UNO) ) + cross_productC ( -vkappa*UNO , cross_productC( conjg(vRadminus*longitud/18.) , vnormalMinus*UNO) )

                    funcREC(j,1,funcBases(i)) = alphaCFIE*conjg(funcRAD(j,1,funcBases(i))) + (1-alphaCFIE)*sum(vRec*octree(niveles)%atheta(:,j))
                    funcREC(j,2,funcBases(i)) = alphaCFIE*conjg(funcRAD(j,2,funcBases(i))) + (1-alphaCFIE)*sum(vRec*octree(niveles)%aphi(:,j))
                end do


            end do




        else
            do i = 1,nBases
                !triangulo mas
                baricentroplus = t_baric(:, e_t(1, funcBases(i)))
                vopuestoplus = p_coord(:,e_po(1,funcBases(i)))
                areaplus = t_area(e_t(1, funcBases(i)))
                longitud = e_long(funcBases(i))
                rhoplus = baricentroplus - vopuestoplus
                !triangulo menos
                baricentrominus = t_baric(:, e_t(2, funcBases(i)))
                vopuestominus = p_coord(:,e_po(2,funcBases(i)))
                areaminus = t_area(e_t(2, funcBases(i)))
                rhominus = vopuestominus - baricentrominus
                !vector radiacion
                do j=1,octree(niveles)%n_kappas
                    vkappa = octree(niveles)%kappas(:,j)
                    vRadplus = exp(jj*ka*dot_product(vkappa,baricentroplus-centro))*rhoplus
                    vRadminus = exp(jj*ka*dot_product(vkappa,baricentrominus-centro))*rhominus

					vnormalPlus = t_normal(:,e_t(1, funcBases(i)))
					vnormalPlus = vnormalPlus/norm2(vnormalPlus)
					vnormalMinus = t_normal(:,e_t(2, funcBases(i)))
					vnormalMinus = vnormalMinus/norm2(vnormalMinus)

                    vRadtotal = (vRadplus + vRadminus)*longitud*0.5

                    funcRAD(j,1,funcBases(i)) = sum(vRadtotal*octree(niveles)%atheta(:,j))
                    funcRAD(j,2,funcBases(i)) = sum(vRadtotal*octree(niveles)%aphi(:,j))

					vRec = cross_productC ( -vkappa*UNO , cross_productC( conjg(vRadplus*longitud*0.5) , vnormalPlus*UNO) ) + cross_productC ( -vkappa*UNO , cross_productC( conjg(vRadminus*longitud*0.5) , vnormalMinus*UNO) )

                    funcREC(j,1,funcBases(i)) = alphaCFIE*conjg(funcRAD(j,1,funcBases(i))) + (1-alphaCFIE)*sum(vRec*octree(niveles)%atheta(:,j))
                    funcREC(j,2,funcBases(i)) = alphaCFIE*conjg(funcRAD(j,2,funcBases(i))) + (1-alphaCFIE)*sum(vRec*octree(niveles)%aphi(:,j))


                end do
            end do

        end if



    end subroutine


    subroutine dividir_en_8_hijos (centro, lcubo, centro_hijo, lcubo_hijo)
         !dumy
        real (kind = dp), dimension(3), intent(in) :: centro
        real (kind = dp), intent(in) :: lcubo
        real (kind = dp), dimension(3,8), intent(out) :: centro_hijo
        real (kind = dp), intent(out) :: lcubo_hijo
        !logical, intent(out) :: fin
        !local
        real (kind = dp) :: lq
        lcubo_hijo = lcubo/2.
        !if (lcubo_hijo <= lambda/2) then
        !    fin = .true.
        !else
        !    fin = .false.
        lq = lcubo/4.
        !centro 1
        centro_hijo(1,1) = centro(1)-lq
        centro_hijo(2,1) = centro(2)-lq
        centro_hijo(3,1) = centro(3)-lq
        !centro 2
        centro_hijo(1,2) = centro(1)-lq
        centro_hijo(2,2) = centro(2)+lq
        centro_hijo(3,2) = centro(3)-lq
        !centro 3
        centro_hijo(1,3) = centro(1)-lq
        centro_hijo(2,3) = centro(2)-lq
        centro_hijo(3,3) = centro(3)+lq
        !centro 4
        centro_hijo(1,4) = centro(1)-lq
        centro_hijo(2,4) = centro(2)+lq
        centro_hijo(3,4) = centro(3)+lq
        !centro 5
        centro_hijo(1,5) = centro(1)+lq
        centro_hijo(2,5) = centro(2)-lq
        centro_hijo(3,5) = centro(3)-lq
        !centro 6
        centro_hijo(1,6) = centro(1)+lq
        centro_hijo(2,6) = centro(2)+lq
        centro_hijo(3,6) = centro(3)-lq
        !centro 7
        centro_hijo(1,7) = centro(1)+lq
        centro_hijo(2,7) = centro(2)-lq
        centro_hijo(3,7) = centro(3)+lq
        !centro 8
        centro_hijo(1,8) = centro(1)+lq
        centro_hijo(2,8) = centro(2)+lq
        centro_hijo(3,8) = centro(3)+lq
        !end if

    end subroutine



    subroutine buscarFuncBases(centro, lcubo, vBusqueda,n_busq,vEncontrados,n_encon,vRestantes,n_rest, encontro)
        !el vector vBusqueda contiene los indices de las funciones bases contenidas (e_centro)
        !dumy
        real (kind = dp), dimension(3), intent(in) :: centro
        real (kind = dp), intent(in) :: lcubo
        integer (kind = il), intent(in), dimension(:) :: vBusqueda
        integer (kind = il), intent(in) :: n_busq
        integer (kind = il), allocatable, intent(out), dimension(:) :: vEncontrados, vRestantes
        integer (kind = il), intent(out) :: n_encon, n_rest
        logical, intent(out) :: encontro

        !local
        real (kind = dp), dimension(6) :: limites
        real (kind = dp) :: lq
        integer (kind = il) :: i
        integer (kind = il), allocatable, dimension(:) :: aux1, aux2
        integer (kind = il) :: n_aux
        logical :: fuera
        !calcular limites (Xmax,Xmin,Ymax,Ymin,Zmax,Zmin) del cubo
        lq = lcubo/2.
        !Xmax
        limites(1) = centro(1)+lq
        !Xmin
        limites(2) = centro(1)-lq
        !Ymax
        limites(3) = centro(2)+lq
        !Ymin
        limites(4) = centro(2)-lq
        !Zmax
        limites(5) = centro(3)+lq
        !Zmin
        limites(6) = centro(3)-lq


        allocate(aux1(n_busq))
        allocate(aux2(n_busq))
        n_encon= 0
        n_rest = 0

        do i=1,n_busq

            fuera = .true.
            if (e_centro(1,vBusqueda(i)) <= limites(1) .and. &
                e_centro(1,vBusqueda(i)) >= limites(2) ) then
                if (e_centro(2,vBusqueda(i)) <= limites(3) .and. &
                    e_centro(2,vBusqueda(i)) >= limites(4) ) then
                    if (e_centro(3,vBusqueda(i)) <= limites(5) .and. &
                        e_centro(3,vBusqueda(i)) >= limites(6) ) then
                        !si esta dentro del cubo
                        !numero de funciones base enonctradas
                        n_encon = n_encon + 1
                        aux1(n_encon) = vBusqueda(i)
                        fuera = .false.
                    end if
                end if
            end if
            if (fuera) then
                !si esta fuera del cubo
                n_rest = n_rest + 1
                aux2(n_rest) = vBusqueda(i)
            end if
        end do
        if (n_encon > 0) then
            encontro = .true.
            allocate(vEncontrados(n_encon))
            allocate(vRestantes(n_rest))
            vEncontrados(:) = aux1(1:n_encon)
            vRestantes(:) = aux2(1:n_rest)
        else
            encontro = .false.
        end if
        deallocate(aux1)
        deallocate(aux2)

    end subroutine

    subroutine llenar_functransf (niveles)
        !dummy
        integer (kind = il), intent(in) :: niveles
        !
        !local
        integer (kind = il) :: n,i,j,k,cont,iPadre,jPadre, contlejT, contCercanas, numJ, iZ, jZ, n_propias
        real (kind = dp), dimension(3) :: rabPadre, rabHijo, resta
        complex (kind = dp), allocatable, dimension(:,:) :: auxTransferencias
        real (kind = dp), allocatable, dimension(:,:) :: auxrabTrans
        integer (kind = il), allocatable, dimension(:,:) :: auxlejT
        integer (kind = il), allocatable, dimension(:) :: auxVcercanos
        real (kind = dp) :: tol
        !integer (kind = il) :: prc, prcN
        !
        n_esparcida_MLFMA = 0
        tol = lambda/100.

        do n = 3, niveles
            print*, 'Llenando funciones de transferencia de nivel ' // inttostr(n)
            !prcN = -1
            call setprc(octree(n)%n_cubos)
            allocate(auxTransferencias(octree(n)%n_kappas,316))
            allocate(auxrabTrans(3,316))
            allocate(auxlejT(2,316))! en realidad el maximo seria menor a 26*8 = 208, pero dejemos 316

            auxrabTrans(:,:) = 0.
            cont = 1
            do i = 1, octree(n)%n_cubos

                contlejT=0
                if (n == niveles) then
                    allocate(auxVcercanos(num_e)) !lo maximo que creemos que habra de funciones bases cercanas
                end if
                contCercanas = 1
                do j = 1, octree(n)%n_cubos

                    !padres
                    iPadre = octree(n)%Cubos(i)%padre
                    jPadre = octree(n)%Cubos(j)%padre
                    rabPadre = octree(n-1)%Cubos(iPadre)%centro-octree(n-1)%Cubos(jPadre)%centro
                    if (norm2(rabPadre) < (2*octree(n-1)%lado_cubo-tol)) then
                        !hijos
                        rabHijo = octree(n)%Cubos(i)%centro-octree(n)%Cubos(j)%centro
                        if (norm2(rabHijo)>= (2*octree(n)%lado_cubo-tol)) then
                            !son cubos lejanos de padres cercanos
                            contlejT = contlejT + 1
                            do k = 1, cont
                                resta = auxrabTrans(:,k) - rabHijo
                                if ((abs(resta(1))<tol).and.(abs(resta(2))<tol).and.(abs(resta(3))<tol)) then
                                    !ya existia
                                    auxlejT(1,contlejT) = j
                                    auxlejT(2,contlejT) = k
                                    exit
                                else if (k == cont) then
                                    !no existe
                                    auxrabTrans(:,cont) = rabHijo
                                    auxTransferencias(:,cont) = calc_functransf(n,rabHijo)
                                    auxlejT(1,contlejT) = j
                                    auxlejT(2,contlejT) = cont
                                    cont = cont + 1
                                end if
                            end do
                        else
                            !son centros cercanos. Aprovechamos de llenar las func bases cercanas para el calculo de Znear
                            if (n == niveles) then
                                !si estoy en el ultimo nivel...

                                numJ = octree(n)%Cubos(j)%n_propias
                                auxVcercanos(contCercanas:(contCercanas + numJ - 1)) = octree(n)%Cubos(j)%bases_propias(:)
                                contCercanas = contCercanas + numJ
                            end if

                        end if

                    end if

                end do

                !guardar lo que tenia en auxlejT
                allocate(octree(n)%Cubos(i)%cubLej_indFuncTransf(2,contlejT))
                octree(n)%Cubos(i)%cubLej_indFuncTransf(:,:) = auxlejT(:,1:(contlejT))
                octree(n)%Cubos(i)%n_cubos_lejanos = contlejT
                !

                !guardar Vcercanos y llenar Zcercana
                if (n == niveles) then
                    !ordenar de forma ascendente bases cercanas
                    call Sort(auxVcercanos, contCercanas-1)

                    allocate(octree(niveles)%Cubos(i)%bases_cercanas(contCercanas-1))
                    octree(niveles)%Cubos(i)%n_cercanas = contCercanas-1
                    octree(niveles)%Cubos(i)%bases_cercanas(1:(contCercanas-1)) = auxVcercanos(1:(contCercanas-1))
                    !allocate(octree(niveles)%Cubos(i)%Zcerc(octree(niveles)%Cubos(i)%n_propias,octree(niveles)%Cubos(i)%n_cercanas))
                    deallocate(auxVcercanos)
                    !tamano de la matriz esparcida
                    n_esparcida_MLFMA = n_esparcida_MLFMA + octree(niveles)%Cubos(i)%n_propias*octree(niveles)%Cubos(i)%n_cercanas

                    !do iZ = 1, octree(niveles)%Cubos(i)%n_propias
                    !    do jZ = 1, octree(niveles)%Cubos(i)%n_cercanas
                    !        octree(niveles)%Cubos(i)%Zcerc(iZ,jZ) = Zmn_MoM(octree(niveles)%Cubos(i)%bases_propias(iZ),octree(niveles)%Cubos(i)%bases_cercanas(jZ))
                    !    end do
                    !end do

                end if

                call updateprc(i)
                !prc = (i*100./octree(n)%n_cubos)
                !if ((prcN/10) /= (prc/10)) then
                !    prcN = prc
                !    print*, inttostr(prc) // ' %'
                !end if
            end do

            allocate(octree(n)%functransf(octree(n)%n_kappas,cont-1))

            octree(n)%functransf(:,:) = auxTransferencias(:,1:(cont-1))

            deallocate(auxTransferencias)
            deallocate(auxrabTrans)
            deallocate(auxlejT)
            print*, 'Finalizado con ' // inttostr(cont-1) // ' funciones de transferencia unicas'

        end do


    end subroutine



    function calc_functransf(nivel, rab) result(v)
        !dummy
        integer (kind = il), intent(in) :: nivel
        real (kind = dp), dimension(3), intent(in) :: rab
        complex (kind =8), dimension(octree(nivel)%n_kappas) :: v
        !
        !local
        integer (kind = il) :: i, j, Ltrue, kDint
        real (kind = dp), dimension(3) :: vkappa
        real (kind = dp) :: kD, modul
        complex(kind = dp) :: valor
        !
        Ltrue = octree(nivel)%n_multipolos
        modul = norm2(rab)
        kD = (2.*pi/getlambda())*modul
        kDint = kD
        if (kDint < Ltrue) then
            Ltrue = kDint
        end if

        do i = 1, octree(nivel)%n_kappas
            vkappa = octree(nivel)%kappas(:,i)
            valor = 0.
            do j=0,Ltrue
                valor = valor + (((-jj)**(j+1))*((2.*j)+1)*hankel(j,kD)*legendre(j,dot_product(vkappa,rab/modul)))
            end do
            valor = valor*(2.*pi/getlambda())/(4.*pi)!/**/
            v(i) = valor
        end do

    end function calc_functransf



    subroutine llenarMinterp(niveles)
        !dummy
        integer (kind = il), intent(in) :: niveles
        !
        !local
        real (kind = dp), dimension(3) :: Wtheta, Wphi, muestTheta, muestPhi !muestTheta/Phi [0,<,>]
        real (kind = dp), dimension(9) :: Wtot
        real (kind = dp) :: XiTheta ,XiPhi
        integer (kind = il) :: i,n,m1,m2, nKfiner, nKcoarser, cont!, nWint
        integer (kind = il), dimension (9) :: imuestCerc
        real (kind  = 8) :: delta, tole
        real (kind  = 8) :: deltaTheta, deltaPhi, dt0,dt1,dt2,dp0,dp1,dp2
        !
        tole = lambda/1000000.

        do n = 3, niveles-1
            nKfiner = octree(n+1)%n_kappas
            nKcoarser = octree(n)%n_kappas
            allocate(octree(n)%mat_interp(2,9*nKcoarser))
            cont = 0
            do i = 1, nKcoarser

                !itera por cada kappa vector del nivel n
                XiTheta = octree(n)%ang_kappas(1,i)
                XiPhi = octree(n)%ang_kappas(2,i)

                imuestCerc = hallarKcercanos(n, i)

                muestTheta(1) = octree(n+1)%ang_kappas(1,imuestCerc(1))
                muestTheta(2) = octree(n+1)%ang_kappas(1,imuestCerc(2))
                muestTheta(3) = octree(n+1)%ang_kappas(1,imuestCerc(3))
                muestPhi(1) = octree(n+1)%ang_kappas(2,imuestCerc(1))
                muestPhi(2) = octree(n+1)%ang_kappas(2,imuestCerc(4))
                muestPhi(3) = octree(n+1)%ang_kappas(2,imuestCerc(7))

                deltaTheta = octree(n+1)%ang_kappas(1,2) - octree(n+1)%ang_kappas(1,1)
                deltaPhi = octree(n+1)%ang_kappas(2,1+ octree(n+1)%n_muestras_theta) - octree(n+1)%ang_kappas(2,1)
                dt0 = abs(xiTheta - muestTheta(1))
                dt1 = abs(xiTheta - muestTheta(2))
                dt2 = abs(xiTheta - muestTheta(3))
                dp0 = abs(XiPhi - muestPhi(1))
                dp1 = abs(XiPhi - muestPhi(2))
                dp2 = abs(XiPhi - muestPhi(3))


                !Prueba de validacion de hallarKcercanos (no valida si los deltas no son constantes (propio de Gauss-Legendre))
                if ( (dt0-tole) >deltaTheta*0.5  .or. (dp0-tole)>deltaPhi*0.5 .or. (dt1-tole)>1.5*deltaTheta .or. (dt2-tole)>1.5*deltaTheta .or. (dp1-tole)>1.5*deltaPhi .or. (dp2-tole)>1.5*deltaPhi) then


                    if ( abs(  muestTheta(1) - pi  )<= (tole + deltaTheta/2)  .or.    abs(  muestTheta(1) - 0.  )<= (tole + deltaTheta/2)   .or.            abs(  muestPhi(1) - 2*pi  )<= (tole + deltaPhi/2)  .or.    abs(  muestPhi(1) - 0.  )<= (tole + deltaPhi/2)             )  then

                    else
                        print*, 'ERROR en hallarKcercanos #######################, nivel: ', n
                        print*, 'kappa ', i

                        print*, 'tita: ', xiTheta
                        print*, 'fi: ', XiPhi
                        print*, 'delta tita: ', deltaTheta
                        print*, 'delta fi: ', deltaPhi

                        print*, 'tita0: ', muestTheta(1)
                        print*, 'tita<: ', muestTheta(2)
                        print*, 'tita>: ', muestTheta(3)
                        print*, 'fi0: ', muestPhi(1)
                        print*, 'fi<: ', muestPhi(2)
                        print*, 'fi>: ', muestPhi(3)

                        print*, 'dtita0: ', dt0
                        print*, 'dtita<: ', dt1
                        print*, 'dtita>: ', dt2
                        print*, 'dfi0: ', dp0
                        print*, 'dfi<: ', dp1
                        print*, 'dfi>: ', dp2
                    end if

                end if

                Wtheta = Wlagrange3(XiTheta,muestTheta)

                WPhi = Wlagrange3(XiPhi,muestPhi)

                do m1 = 1,3
                    do m2 = 1,3
                        cont = cont + 1
                        !Wtot(cont) = Wtheta(m1)*WPhi(m2)
                        octree(n)%mat_interp(1,cont) = Wtheta(m2)*WPhi(m1)
                    end do
                end do
                octree(n)%mat_interp(2,(cont-8):(cont)) = imuestCerc

            end do

        end do

    end subroutine



    function hallarKcercanos(nivel, indice, nivelmuestras) result(V)


        !dummy
        integer (kind = il), intent(in) :: indice, nivel !NIVEL INDICA EL NIVEL MAS GRUESO(bajo) A DONDE QUEREMOS INTERPOLAR
        integer (kind = il), intent(in), optional :: nivelmuestras
        integer (kind = il), dimension (9) :: V
        !
        !local
        real (kind = dp) :: Dphi, Dtheta, phi, theta, residuo, del, valorpos
        integer (kind = dp) :: Lmult, posPhi, posPhimay, posPhimen, n_muestras_phi, n_muestras_theta, posTheta, posThetamay, posThetamen, techo, piso, posmid, Lmuest
        !
        if (present(nivelmuestras)) then
            Lmuest = nivelmuestras
        else
            Lmuest = nivel + 1
        end if
        Lmult = octree(Lmuest)%n_multipolos
        n_muestras_phi = octree(Lmuest)%n_muestras_phi
        n_muestras_theta = octree(Lmuest)%n_muestras_theta
        Dphi = (2.*pi)/n_muestras_phi

        phi = octree(nivel)%ang_kappas(2,indice)
        theta = octree(nivel)%ang_kappas(1,indice)

        techo = 1
        piso = n_muestras_phi
        del = octree(Lmuest)%ang_kappas(2,n_muestras_theta + 1) -  octree(Lmuest)%ang_kappas(2,1)

        do
            posmid = (techo + piso)/2
            valorpos = octree(Lmuest)%ang_kappas(2,(posmid-1)*n_muestras_theta + 1)

            if ( (posmid == techo) .or. (posmid == piso) ) then
                if (  abs( phi - octree(Lmuest)%ang_kappas(2,(techo-1)*n_muestras_theta + 1) )   <   abs( phi - octree(Lmuest)%ang_kappas(2,(piso-1)*n_muestras_theta + 1) )      ) then
                    !techo es mejor que piso
                    posmid = techo
                else
                    posmid = piso
                end if
                exit
            end if


            if (valorpos == phi) then
                !cayo justo en el valor exacto
                exit
            end if
            if (del > 0) then
                !vector creciente

                if (valorpos > phi) then
                    piso = posmid
                elseif (valorpos < phi) then
                    techo = posmid
                end if

            else
                !vector decreciente

                if (valorpos > phi) then
                    techo = posmid
                elseif (valorpos < phi) then
                    piso = posmid
                end if

            end if


        end do

        posPhi = posmid
        !TRABAJAR AQUI CON EL posmid DE PHI

        if (del > 0) then
            !vector creciente
            if (posPhi == 1) then
                posPhimay = posPhi + 1
                posPhimen = n_muestras_phi
            else if (posPhi == n_muestras_phi) then
                posPhimay = 1
                posPhimen = posPhi - 1
            else
                posPhimen = posPhi - 1
                posPhimay = posPhi + 1
            end if


        else
            !vector decreciente
            if (posPhi == 1) then
                posPhimen = posPhi + 1
                posPhimay = n_muestras_phi
            else if (posPhi == n_muestras_phi) then
                posPhimen = 1
                posPhimay = posPhi - 1
            else
                posPhimay = posPhi - 1
                posPhimen = posPhi + 1
            end if


        end if





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






        !Busqueda del THETA más cercano
        techo = 1
        piso = n_muestras_theta
        del = octree(Lmuest)%ang_kappas(1,2) -  octree(Lmuest)%ang_kappas(1,1)
        !si del > 0 entonces el vector va creciendo (tipico de phi)
        !si del < 0 entonces el vector va decreciendo (tipico de theta)
        do
            posmid = (techo + piso)/2
            valorpos = octree(Lmuest)%ang_kappas(1,posmid)

            if ( (posmid == techo) .or. (posmid == piso) ) then
                if (  abs( theta - octree(Lmuest)%ang_kappas(1,techo) )   <   abs( theta - octree(Lmuest)%ang_kappas(1,piso) )      ) then
                    !techo es mejor que piso
                    posmid = techo
                else
                    posmid = piso
                end if
                exit
            end if


            if (valorpos == theta) then
                !cayo justo en el valor exacto
                exit
            end if
            if (del > 0) then
                !vector creciente

                if (valorpos > theta) then
                    piso = posmid
                elseif (valorpos < theta) then
                    techo = posmid
                end if

            else
                !vector decreciente

                if (valorpos > theta) then
                    techo = posmid
                elseif (valorpos < theta) then
                    piso = posmid
                end if

            end if


        end do
        posTheta = posmid
        !TRABAJAR AQUI CON EL posmid DE THETA

        if (del > 0) then
            !vector creciente
            if (posTheta == 1) then
                posThetamen = n_muestras_theta
                posThetamay = posTheta + 1
            else if (posTheta == n_muestras_theta) then
                posThetamen = posTheta - 1
                posThetamay = 1
            else
                posThetamay = posTheta + 1
                posThetamen = posTheta - 1
            end if


        else
            !vector decreciente
            if (posTheta == 1) then
                posThetamay = n_muestras_theta
                posThetamen = posTheta + 1
            else if (posTheta == n_muestras_theta) then
                posThetamay = posTheta - 1
                posThetamen = 1
            else
                posThetamen = posTheta + 1
                posThetamay = posTheta - 1
            end if


        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !YA TENEMOS LOS INDICES MEDIOS MAYORES Y MENORES DE THETA Y PHI INDEPENDIENTEMENTE

        !Guardar indices:
        !guardar phi0
        V(1) = n_muestras_theta*(posPhi-1) + posTheta !iphi0theta0
        V(2) = n_muestras_theta*(posPhi-1) + posThetamen !iphi0theta<
        V(3) = n_muestras_theta*(posPhi-1) + posThetamay !iphi0theta>
        !
        !guardar phi<
        V(4) = n_muestras_theta*(posPhimen-1) + posTheta !iphi<theta0
        V(5) = n_muestras_theta*(posPhimen-1) + posThetamen !iphi<theta<
        V(6) = n_muestras_theta*(posPhimen-1) + posThetamay !iphi<theta>
        !
        !guardar phi>
        V(7) = n_muestras_theta*(posPhimay-1) + posTheta !iphi>theta0
        V(8) = n_muestras_theta*(posPhimay-1) + posThetamen !iphi>theta<
        V(9) = n_muestras_theta*(posPhimay-1) + posThetamay !iphi>theta>
        !
    !print*, V


    end function hallarKcercanos


    function Wlagrange3(x, xmuestra) result(Wi)
        !dummy
        real (kind = dp), intent(in) :: x
        real (kind = dp), intent(in), dimension(3) :: xmuestra
        real (kind = dp), dimension(3) :: Wi
        !
        !local
        integer (kind = il) :: i, j
        real (kind = dp) :: ww
        !
        do i = 1, 3

            !w = 0.
            ww = 1.
            do j = 1, 3

                if (j /= i) then
                    ww = ww*( (x-xmuestra(j)) / (xmuestra(i) - xmuestra(j)) )
                end if

            end do
            Wi(i) = ww
        end do
    end function Wlagrange3

    function Interpolar(Mint, nivel, Datos) result(Vint)
        !dummy
        real (kind = dp), dimension(:,:),intent(in) :: Mint
        integer (kind = il),intent(in) :: nivel
        complex(kind = dp), intent(in), dimension(:,:) :: Datos
        complex(kind = dp), dimension(octree(nivel)%n_kappas,2) :: Vint
        !
        !local
        integer (kind = il) :: nFilas, cont, i, j
        !
        nFilas = octree(nivel)%n_kappas
        cont = 0
        Vint(:,:) = 0.

        do i = 1, nFilas

            do j = 1, 9
                cont = cont + 1
                Vint(i,:) = Vint(i,:) + ( Mint(1, cont)*Datos(Mint(2,cont),:) )
            end do

        end do
    end function Interpolar




    function Anterpolar(Mint, nivel, Datos) result(Vant)
        !dummy
        real (kind = dp), dimension(:,:),intent(in) :: Mint
        integer (kind = il),intent(in) :: nivel
        complex(kind = dp), intent(in), dimension(:,:) :: Datos
        complex(kind = dp), dimension(octree(nivel + 1)%n_kappas,2) :: Vant
        !
        !local
        integer (kind = il) :: nFil, nCol,cont, fila, col
        !
        cont = 0
        nFil = octree(nivel)%n_kappas
        Vant(:,:) = 0.
        do fila = 1, Nfil
            do col = 1,9
                cont = cont + 1
                Vant(Mint(2,cont),:) = Vant(Mint(2,cont),:) + ( Mint(1,cont)*Datos(fila,:) )
            end do
        end do

    end function Anterpolar


    subroutine llenarZesparcida_gabriel3()
        !local
        integer (kind = il) :: i,j,c,bp,bc,nc,m,n,fila,cont, prc, prca
        complex (kind = dp), allocatable, dimension(:) :: a_znear
        integer (kind = il), allocatable, dimension (:) :: i_znear, j_znear
        integer (kind = il) :: cub, nf, f, indj, indice, cublej, indcentro, indn, indm, mm
        !

        prca = -1
        fila = 1
        cont = 0

        allocate(a_znear(n_esparcida_MLFMA))
        allocate(i_znear(n_esparcida_MLFMA))
        allocate(j_znear(n_esparcida_MLFMA))


        a_znear(:) = 0.
        i_znear(:) = 0
        j_znear(:) = 0


        do cub = 1, octree(nivel_max)%n_cubos
            nf = octree(nivel_max)%Cubos(cub)%n_propias
            do f = 1, nf
                indice = octree(nivel_max)%Cubos(cub)%bases_propias(f)

                do m = 1, octree(nivel_max)%Cubos(cub)%n_cercanas
                    indj = octree(nivel_max)%Cubos(cub)%bases_cercanas(m)
                    cont = cont + 1
                    a_znear(cont) = Zmn_MoM(indice, indj)
                    i_znear(cont) = indice
                    j_znear(cont) = indj
                end do

            end do

        prc = 100*cub/octree(nivel_max)%n_cubos
        if (prc/10 /= prca/10) then
            print*, inttostr(prc) // ' %'
            prca = prc
        end if

        end do
        

        allocate(Zesparcida_MLFMA(n_esparcida_MLFMA))
        allocate(ja_MLFMA(n_esparcida_MLFMA))
        allocate(ia_MLFMA(num_e + 1))

        Zesparcida_MLFMA(:) = 0.
        ja_MLFMA(:) = 0
        ia_MLFMA(:) = 0
print*, 'cont: ', cont
print*, 'n_esparcida_MLFMA: ', n_esparcida_MLFMA
        call coocsrC ( num_e, n_esparcida_MLFMA, a_znear, i_znear, j_znear, Zesparcida_MLFMA, ja_MLFMA, ia_MLFMA )

        deallocate(a_znear)
        deallocate(i_znear)
        deallocate(j_znear)

    end subroutine llenarZesparcida_gabriel3



    subroutine llenarZesparcida_gabriel2()
        !local
        integer (kind = il) :: i,j,c,bp,bc,nc,m,n,fila,cont, prc, prca
        complex (kind = dp), allocatable, dimension(:) :: a_znear
        integer (kind = il), allocatable, dimension (:) :: i_znear, j_znear
        integer (kind = il) :: cub, nf, f, indj, indice, cublej, indcentro, indn, indm, mm
        !

        prca = -1
        fila = 1
        cont = 0

        n_esparcida_MLFMA = num_e*num_e

        allocate(a_znear(n_esparcida_MLFMA))
        allocate(i_znear(n_esparcida_MLFMA))
        allocate(j_znear(n_esparcida_MLFMA))


        a_znear(:) = 0.
        i_znear(:) = 0
        j_znear(:) = 0


        do cub = 1, octree(nivel_max)%n_cubos
            nf = octree(nivel_max)%Cubos(cub)%n_propias
            do f = 1, nf
                indice = octree(nivel_max)%Cubos(cub)%bases_propias(f)

                do m = 1, octree(nivel_max)%Cubos(cub)%n_cercanas
                    indj = octree(nivel_max)%Cubos(cub)%bases_cercanas(m)
                    cont = (indice - 1)*num_e + indj
                    a_znear(cont) = Zmn_MoM(indice, indj)
                end do

            end do

        prc = 100*cub/octree(nivel_max)%n_cubos
        if (prc/10 /= prca/10) then
            print*, inttostr(prc) // ' %'
            prca = prc
        end if

        end do

        do m = 1, num_e
            do f = 1, num_e
                cont = (m - 1)*num_e + f
                i_znear(cont) = m
                j_znear(cont) = f
            end do
        end do
            
        

        allocate(Zesparcida_MLFMA(n_esparcida_MLFMA))
        allocate(ja_MLFMA(n_esparcida_MLFMA))
        allocate(ia_MLFMA(num_e + 1))


        Zesparcida_MLFMA(:) = 0.
        ja_MLFMA(:) = 0
        ia_MLFMA(:) = 0
print*, 'cont: ', cont
print*, 'n_esparcida_MLFMA: ', n_esparcida_MLFMA
        call coocsrC ( num_e, n_esparcida_MLFMA, a_znear, i_znear, j_znear, Zesparcida_MLFMA, ja_MLFMA, ia_MLFMA )

        deallocate(a_znear)
        deallocate(i_znear)
        deallocate(j_znear)

    end subroutine llenarZesparcida_gabriel2




    subroutine llenarZesparcida_gabriel()
        !local
        integer (kind = il) :: i,j,c,bp,bc,nc,m,n,fila,cont, prc, prca
        complex (kind = dp), allocatable, dimension(:) :: a_znear
        integer (kind = il), allocatable, dimension (:) :: i_znear, j_znear
        !

        prca = -1
        fila = 1
        cont = 0

        allocate(a_znear(n_esparcida_MLFMA))
        allocate(i_znear(n_esparcida_MLFMA))
        allocate(j_znear(n_esparcida_MLFMA))
        a_znear(:) = 0.
        i_znear(:) = 0
        j_znear(:) = 0

        do m = 1, num_e
            Busq: do c = 1,octree(nivel_max)%n_cubos
                do bp = 1, octree(nivel_max)%Cubos(c)%n_propias
                    if (octree(nivel_max)%Cubos(c)%bases_propias(bp) == m) then
                        nc = octree(nivel_max)%Cubos(c)%n_cercanas

                        do j = 1, nc
                            cont = cont + 1
                            n = octree(nivel_max)%Cubos(c)%bases_cercanas(j)
                            a_znear(cont) = Zmn_MoM(m,n)
                            j_znear(cont) = n
                            i_znear(cont) = m
                        end do
                        exit Busq
                    end if
                end do
            end do Busq
            prc = 100*m/num_e
            if (prc/10 /= prca/10) then
                print*, inttostr(prc) // ' %'
                prca = prc
            end if
        end do

        allocate(Zesparcida_MLFMA(n_esparcida_MLFMA))
        allocate(ja_MLFMA(n_esparcida_MLFMA))
        allocate(ia_MLFMA(num_e + 1))
        Zesparcida_MLFMA(:) = 0.
        ja_MLFMA(:) = 0
        ia_MLFMA(:) = 0
print*, 'cont: ', cont
print*, 'n_esparcida_MLFMA: ', n_esparcida_MLFMA
        call coocsrC ( num_e, n_esparcida_MLFMA, a_znear, i_znear, j_znear, Zesparcida_MLFMA, ja_MLFMA, ia_MLFMA )

        deallocate(a_znear)
        deallocate(i_znear)
        deallocate(j_znear)

    end subroutine llenarZesparcida_gabriel



    subroutine llenarZesparcida()
        !local
        integer (kind = il) :: i,j,c,bp,bc,nc,m,n,fila,cont
        !

        fila = 1
        cont = 0

        allocate(Zesparcida_MLFMA(n_esparcida_MLFMA))
        allocate(ja_MLFMA(n_esparcida_MLFMA))
        allocate(ia_MLFMA(num_e + 1))

        call setprc(num_e)
        do m = 1, num_e
            Busq: do c = 1,octree(nivel_max)%n_cubos
                do bp = 1, octree(nivel_max)%Cubos(c)%n_propias
                    if (octree(nivel_max)%Cubos(c)%bases_propias(bp) == m) then
                        nc = octree(nivel_max)%Cubos(c)%n_cercanas
                        ia_MLFMA(m) = fila
                        do j =1, nc
                            cont = cont + 1
                            n = octree(nivel_max)%Cubos(c)%bases_cercanas(j)
                            Zesparcida_MLFMA(cont) = Zmn_MoM(m,n)
                            ja_MLFMA(cont) = n
                            fila = fila + 1
                        end do
                        exit Busq
                    end if
                end do
            end do Busq
            call updateprc(m)
            !prc = 100*m/num_e
            !if (prc/10 /= prca/10) then
            !    print*, inttostr(prc) // ' %'
            !    prca = prc
            !end if
        end do
        ia_MLFMA(num_e + 1) = fila


    end subroutine llenarZesparcida


    subroutine test_interp()!#prueba
        complex (kind = dp), allocatable, dimension(:,:) :: cOrig, cInter, cAnter, orig_t, inter_t, anter_t, orig_p, inter_p, anter_p, cMet1, cMet2, met1_t, met1_p, met2_t, met2_p
        real (kind = dp), allocatable, dimension(:,:) :: W
        real (kind = dp), allocatable, dimension(:) :: thetas_n, phis_n, thetas_n1, phis_n1
        real (kind = dp), allocatable, dimension(:,:) :: func1, func2, func2n
        
        integer (kind = il) :: niv, cub, i, j, cont

        num_iter_current = num_iter_current + 1

        if (num_iter_current == 10) then
                print*, 'iniciando prueba interpolacion'

                cont = 0
                niv = 4
                cub = 3

                allocate(W(octree(niv-1)%n_kappas,octree(niv)%n_kappas))
                W(:,:) = 0.
print*, 'a'
                do i = 1, octree(niv-1)%n_kappas
                    do j = 1,9
                        cont = cont + 1
                        W(i,octree(niv-1)%mat_interp(2,cont)) = octree(niv-1)%mat_interp(1,cont)
                    end do
                end do
print*, 'b'
                allocate(cOrig(octree(niv)%n_kappas,2))
                allocate(cInter(octree(niv-1)%n_kappas,2))
                allocate(cAnter(octree(niv)%n_kappas,2))

                cOrig = octree(niv)%campoC(:,:,cub)
                cInter = Interpolar(octree(niv-1)%mat_interp,niv-1,cOrig)
                cAnter = Anterpolar(octree(niv-1)%mat_interp, niv-1, cInter)
print*, 'c'
                allocate(orig_t(octree(niv)%n_muestras_theta,octree(niv)%n_muestras_phi))
                allocate(anter_t(octree(niv)%n_muestras_theta,octree(niv)%n_muestras_phi))
                allocate(inter_t(octree(niv-1)%n_muestras_theta,octree(niv-1)%n_muestras_phi))

                allocate(orig_p(octree(niv)%n_muestras_theta,octree(niv)%n_muestras_phi))
                allocate(anter_p(octree(niv)%n_muestras_theta,octree(niv)%n_muestras_phi))
                allocate(inter_p(octree(niv-1)%n_muestras_theta,octree(niv-1)%n_muestras_phi))

                allocate(thetas_n(octree(niv)%n_muestras_theta))
                allocate(thetas_n1(octree(niv - 1)%n_muestras_theta))
                allocate(phis_n(octree(niv)%n_muestras_phi))
                allocate(phis_n1(octree(niv - 1)%n_muestras_phi))



                allocate(func1(octree(niv - 1)%n_kappas,2))
                allocate(func2(octree(niv - 1)%n_kappas,2))
                allocate(func2n(octree(niv)%n_kappas,2))

                allocate(cMet1(octree(niv)%n_kappas,2))
                allocate(cMet2(octree(niv)%n_kappas,2))
                allocate(met1_t(octree(niv)%n_muestras_theta,octree(niv)%n_muestras_phi))
                allocate(met2_t(octree(niv)%n_muestras_theta,octree(niv)%n_muestras_phi))
                allocate(met1_p(octree(niv)%n_muestras_theta,octree(niv)%n_muestras_phi))
                allocate(met2_p(octree(niv)%n_muestras_theta,octree(niv)%n_muestras_phi))


                func1(:,1) = sin(3*octree(niv - 1)%ang_kappas(1,:)) + 2*cos(1.3*octree(niv - 1)%ang_kappas(2,:))
                func2(:,1) = 5*cos(7*octree(niv - 1)%ang_kappas(1,:)) + 2*sin(5.4*octree(niv - 1)%ang_kappas(2,:))
                func2n(:,1) = 5*cos(7*octree(niv)%ang_kappas(1,:)) + 2*sin(5.4*octree(niv)%ang_kappas(2,:))
                func1(:,2) = sin(3*octree(niv - 1)%ang_kappas(1,:)) + 2*cos(1.3*octree(niv - 1)%ang_kappas(2,:))
                func2(:,2) = 5*cos(7*octree(niv - 1)%ang_kappas(1,:)) + 2*sin(5.4*octree(niv - 1)%ang_kappas(2,:))
                func2n(:,2) = 5*cos(7*octree(niv)%ang_kappas(1,:)) + 2*sin(5.4*octree(niv)%ang_kappas(2,:))                

                cMet1 = Anterpolar(octree(niv-1)%mat_interp, niv-1, cInter*func1*func2)
                cMet2 = func2n*Anterpolar(octree(niv-1)%mat_interp, niv-1, cInter*func1)

                cont=0
                do i = 1, octree(niv)%n_muestras_phi
                    do j = 1, octree(niv)%n_muestras_theta
                        cont  = cont + 1
                        orig_t(j,i) = cOrig(cont,1)
                        orig_p(j,i) = cOrig(cont,2)

                        anter_t(j,i) = cAnter(cont,1)
                        anter_p(j,i) = cAnter(cont,2)

                        met1_t(j,i) = cMet1(cont,1)
                        met1_p(j,i) = cMet1(cont,2)

                        met2_t(j,i) = cMet2(cont,1)
                        met2_p(j,i) = cMet2(cont,2)                        

                        thetas_n(j) = octree(niv)%ang_kappas(1,cont)
                        phis_n(i) = octree(niv)%ang_kappas(2,cont)
                    end do
                end do

print*, 'd'
                cont=0
                do i = 1, octree(niv-1)%n_muestras_phi
                    do j = 1, octree(niv-1)%n_muestras_theta
                        cont  = cont + 1
                        inter_t(j,i) = cInter(cont,1)
                        inter_p(j,i) = cInter(cont,2)

                        thetas_n1(j) = octree(niv-1)%ang_kappas(1,cont)
                        phis_n1(i) = octree(niv-1)%ang_kappas(2,cont)
                    end do
                end do                
print*, 'e'
                call guardarVectorR(thetas_n,'_thetas_n.dat')
                call guardarVectorR(thetas_n1,'_thetas_n1.dat')
                call guardarVectorR(phis_n,'_phis_n.dat')
                call guardarVectorR(phis_n1,'_phis_n1.dat')

                call guardarMatrizR(real(orig_t),'_orig_t_R.dat');
                call guardarMatrizR(real(orig_p),'_orig_p_R.dat');
                call guardarMatrizR(real(anter_t),'_anter_t_R.dat');
                call guardarMatrizR(real(anter_p),'_anter_p_R.dat');
                call guardarMatrizR(real(inter_t),'_inter_t_R.dat');
                call guardarMatrizR(real(inter_p),'_inter_p_R.dat');
                call guardarMatrizR(real(met1_t),'_met1_t_R.dat');
                call guardarMatrizR(real(met1_p),'_met1_p_R.dat');
                call guardarMatrizR(real(met2_t),'_met2_t_R.dat');
                call guardarMatrizR(real(met2_p),'_met2_p_R.dat');


                call guardarMatrizR(aimag(orig_t),'_orig_t_I.dat');
                call guardarMatrizR(aimag(orig_p),'_orig_p_I.dat');
                call guardarMatrizR(aimag(anter_t),'_anter_t_I.dat');
                call guardarMatrizR(aimag(anter_p),'_anter_p_I.dat');
                call guardarMatrizR(aimag(inter_t),'_inter_t_I.dat');
                call guardarMatrizR(aimag(inter_p),'_inter_p_I.dat');
                call guardarMatrizR(aimag(met1_t),'_met1_t_I.dat');
                call guardarMatrizR(aimag(met1_p),'_met1_p_I.dat');
                call guardarMatrizR(aimag(met2_t),'_met2_t_I.dat');
                call guardarMatrizR(aimag(met2_p),'_met2_p_I.dat');                

                call guardarMatrizR(W, '_Winterp.dat')



                print*, 'fin guardado interp'

        end if

    end subroutine



    subroutine llenar_campoCyB(xI,niveles)
        !dummy
        complex (kind = dp), dimension(:) :: xI
        integer (kind = il) :: niveles
        !
        !local
        integer (kind = il) :: n,cub,nf,f,m,parent,indCentro,indTrans, numdekappas
        real (kind = dp), dimension(3) :: r
        real (kind = dp), allocatable, dimension(:) :: Kr
        complex (kind = dp), allocatable, dimension(:) :: expKr
        complex (kind = dp), allocatable, dimension(:,:) :: Agtras

        !

        do n = 3, niveles
            octree(n)%campoC(:,:,:) = 0.
            octree(n)%campoB(:,:,:) = 0.
        end do
        !

        do n = niveles,4,-1

            if (n == niveles) then
                do cub = 1, octree(niveles)%n_cubos
                    nf = octree(niveles)%Cubos(cub)%n_propias
                    do f = 1,nf
                        octree(niveles)%campoC(:,:,cub) = octree(niveles)%campoC(:,:,cub) + funcRAD(:,:,octree(niveles)%Cubos(cub)%bases_propias(f))*xI(octree(niveles)%Cubos(cub)%bases_propias(f))
                    end do
                end do
            end if
            numdekappas = octree(n-1)%n_kappas
            allocate(Kr(numdekappas))
            allocate(Agtras(octree(n-1)%n_kappas,2))
            allocate(expKr(octree(n-1)%n_kappas))

            do cub = 1,octree(n)%n_cubos

                parent = octree(n)%Cubos(cub)%padre
                r = octree(n)%Cubos(cub)%centro - octree(n-1)%Cubos(parent)%centro
                do m = 1,octree(n-1)%n_kappas
                    Kr(m) = dot_product(octree(n-1)%kappas(:,m),r)
                end do
                expKr = exp(jj*((2.*pi)/lambda)*Kr)
                Agtras = Interpolar(octree(n-1)%mat_interp,n-1,octree(n)%campoC(:,:,cub))

                octree(n-1)%campoC(:,1,parent) = octree(n-1)%campoC(:,1,parent) + Agtras(:,1)*expKr
                octree(n-1)%campoC(:,2,parent) = octree(n-1)%campoC(:,2,parent) + Agtras(:,2)*expKr

            end do
            deallocate(Kr)
            deallocate(Agtras)
            deallocate(expKr)
        end do
!print*, 'salida'
        !llenar campo B
        do n = 3, niveles
            do cub = 1,octree(n)%n_cubos
                do m = 1,octree(n)%Cubos(cub)%n_cubos_lejanos
                    indCentro = octree(n)%Cubos(cub)%cubLej_indFuncTransf(1,m)
                    indTrans = octree(n)%Cubos(cub)%cubLej_indFuncTransf(2,m)

                    octree(n)%campoB(:,1,cub) = octree(n)%campoB(:,1,cub) + octree(n)%campoC(:,1,indCentro)*octree(n)%functransf(:,indTrans)
                    octree(n)%campoB(:,2,cub) = octree(n)%campoB(:,2,cub) + octree(n)%campoC(:,2,indCentro)*octree(n)%functransf(:,indTrans)

                end do

                octree(n)%campoB(:,1,cub) = octree(n)%campoB(:,1,cub)*octree(n)%pesos_integ(:)
                octree(n)%campoB(:,2,cub) = octree(n)%campoB(:,2,cub)*octree(n)%pesos_integ(:)

            end do
        end do

    end subroutine

    subroutine llenar_campoA(niveles)
        !dummy
        integer (kind = il) :: niveles
        !
        !local
        integer (kind = il) :: n,cub,parent, m
        complex (kind = dp), allocatable, dimension(:,:) :: datos, Anter!!!temporal
        real (kind = dp), dimension(3) :: r
        real (kind = dp), allocatable, dimension(:) :: Kr
        complex (kind = dp), allocatable, dimension(:) :: expKr
        !
        do n = 4, niveles
            allocate(datos(octree(n-1)%n_kappas,2))
            allocate(Kr(octree(n-1)%n_kappas))
            allocate(expKr(octree(n-1)%n_kappas))
            do cub = 1, octree(n)%n_cubos
                parent = octree(n)%Cubos(cub)%padre
                r = octree(n)%Cubos(cub)%centro - octree(n-1)%Cubos(parent)%centro
                do m = 1, octree(n-1)%n_kappas
                    Kr(m) = dot_product(octree(n-1)%kappas(:,m),r)
                end do
                expKr = exp(-jj*((2.*pi)/lambda)*Kr)
                if (n == 4) then
                    datos(:,1) = octree(n-1)%campoB(:,1,parent)*expKr
                    datos(:,2) = octree(n-1)%campoB(:,2,parent)*expKr
                else
                    datos(:,1) = octree(n-1)%campoA(:,1,parent)*expKr
                    datos(:,2) = octree(n-1)%campoA(:,2,parent)*expKr
                end if

                octree(n)%campoA(:,:,cub) = octree(n)%campoB(:,:,cub) + Anterpolar(octree(n-1)%mat_interp, n-1, datos)

            end do
            deallocate(datos)
            deallocate(Kr)
            deallocate(expkr)
        end do



    end subroutine



    function MatxVec_MLFMA(vecX) result(res)
        !dummy
        complex (kind = dp), intent(in), dimension(:) :: vecX
        complex (kind = dp), dimension(size(vecX,1)) :: res
        !
        !local
        integer (kind = il) :: cub, nf, indice,f, q, m, indj, cont
        !

        call llenar_campoCyB(vecX, nivel_max)
        !print*, 'Listo agregación'
        call llenar_campoA(nivel_max)
        !print*, 'Listo disgregación'

        do cub = 1, octree(nivel_max)%n_cubos
            nf = octree(nivel_max)%Cubos(cub)%n_propias
            do f = 1, nf
                indice = octree(nivel_max)%Cubos(cub)%bases_propias(f)

                
                !res(indice) = sum(  sum( funcREC(:,:,indice)*octree(nivel_max)%campoA(:,:,cub)  ,DIM = 2)*octree(nivel_max)%pesos_integ(:), DIM=1  ) !#mosca   
                res(indice) = sum(  sum( funcREC(:,:,indice)*octree(nivel_max)%campoA(:,:,cub)  ,DIM = 2),DIM=1  )
                res(indice) = res(indice)*(jj*w*mu/(4.*pi))


!                do m = 1, octree(nivel_max)%Cubos(cub)%n_cercanas
!                    indj = octree(nivel_max)%Cubos(cub)%bases_cercanas(m)
!                    res(indice) = res(indice) + (  vecX(indj)*octree(nivel_max)%Cubos(cub)%Zcerc(f,m)  )
!                end do

            end do

        end do

        !MxV con Zcercana
        cont = 0
        do m = 1, num_e
            do nf = 1, ia_MLFMA(m+1)-ia_MLFMA(m)
                cont = cont + 1
                res(m) = res(m) +  vecX(ja_MLFMA(cont))*Zesparcida_MLFMA(cont)
            end do
        end do

        call test_interp() !#prueba

    end function MatxVec_MLFMA





end module modMLFMA

