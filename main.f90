program main
    use modGlobalParam
    use modGlobalMetod
    use modReadStl
    use modRWG
    use modMoM
    use modRCS
    use modIterativo
    use modFMM
    use modMLFMA
    implicit none
!NO FUNCIONAAAAAAAAAAAaa
!Mariely te va a espiar
!!!!!!CON MI DEVELOP
!SAPO
    complex (kind = dp_prec), allocatable, dimension(:) :: precond_alu
    integer (kind = il), allocatable, dimension(:) :: precond_jlu
    integer (kind = il), allocatable, dimension(:) :: precond_ju

    complex (kind = dp), allocatable, dimension(:,:) :: Z!, Zinv
    complex ( kind = dp ),  pointer, dimension (:,:) :: Z_MatxVec

    complex (kind = dp), allocatable, dimension(:) :: I
    complex (kind = dp), allocatable, dimension(:) :: b
    complex (kind = dp), allocatable, dimension(:) :: guess
    real ( kind = dp ), allocatable, dimension(:,:) :: p_coord
    integer ( kind = il ), allocatable, dimension(:,:) :: t_p
    integer ( kind = il ), allocatable, dimension(:,:) :: e_p
    integer ( kind = il ), allocatable, dimension(:,:) :: e_t
    integer ( kind = il ), allocatable, dimension(:,:) :: e_po
    real (kind = dp ), allocatable, dimension(:) :: e_long
    real ( kind = dp ), allocatable, dimension(:) :: t_area
	real ( kind= dp ), allocatable, dimension(:,:) :: t_normal
    real ( kind = dp ), allocatable, dimension(:,:) :: e_centro
    real ( kind = dp ), allocatable, dimension(:,:) :: t_baric
    real ( kind = dp ), allocatable, dimension(:,:,:) :: t_baric_sub

    real (kind = dp) :: frequency

    integer ( kind = il ) :: num_e
    integer ( kind = il ) :: num_p
    integer ( kind = il ) :: num_t


    complex (kind = dp), dimension(3) :: pol_onda
    complex (kind = dp) :: ctte_onda

    real (kind = dp) :: error_solver
    real (kind = dp), dimension(3) :: dir_onda
    
    character (len = 20) :: comm, val
    character (len = 3), parameter :: strdefault= '/%\'


    integer (kind = il) :: ent, numtest


!real (kind = 8), dimension (10) :: v1,v2
!v1(:) = 2.8735
!v2(:) = 8.3212999
!call test(90017288,4,v1,v2)
    !call testSAAD()

    !Inicializacion
!sapo mi prueba de committ
!claro que no me va a espiar sopo

    ctte_onda = 1.
    lambda = 1.
    call setlambda(lambda)
    !skipkD = .false.


    conf%test_name = strdefault
    conf%msj = .false.
    conf%numThetaRCS = 1
    conf%numPhiRCS = 360
    conf%errorSolverCGS = 0.001
    conf%Lforz = -1
    conf%dirTheta = 90.
    conf%dirPhi = 180.
    conf%rcsMonoFMin = 10
    conf%rcsMonoFMax = 275
    conf%rcsMonoSamples = 75
    conf%usar_precond = 'ILUT'
    conf%prec_activo = .true.
    conf%fill_in_prc = 5.
    conf%drop_tolerance = 0.001
    conf%analisis_rcs = 'MONO'
    conf%n_pasos_rcs_mono = 180
    conf%tipo_polarizac = 'VV'
    conf%variable_espec = 'FREC'
    conf%file_name = strdefault
    conf%file_script = strdefault
    conf%criterio_ladomax = 0.5_dp

    !conf%polVH = 1
    error_solver = conf%errorSolverCGS
    dir_onda(1) = sin(conf%dirTheta*pi/180.)*cos(conf%dirPhi*pi/180.)
    dir_onda(2) = sin(conf%dirTheta*pi/180.)*sin(conf%dirPhi*pi/180.)
    dir_onda(3) = cos(conf%dirTheta*pi/180.)
    call ajustar_polariz(dir_onda,pol_onda)






call waitcommands()

contains



subroutine destruir_todo()

    call destruir_MLFMA()
    call destruir_FMM()

    if (allocated(p_coord)) then
        deallocate(p_coord)
    end if

    if (allocated(guess)) then
        deallocate(guess)
    end if

    if (allocated(b)) then
        deallocate(b)
    end if

    if (allocated(I)) then
        deallocate(I)
    end if

    if (allocated(Z)) then
        deallocate(Z)
    end if

    if (allocated(precond_ju)) then
        deallocate(precond_ju)
    end if

    if (allocated(precond_jlu)) then
        deallocate(precond_jlu)
    end if

    if (allocated(precond_alu)) then
        deallocate(precond_alu)
    end if

    if (allocated(t_p)) then
        deallocate(t_p)
    end if

    if (allocated(e_p)) then
        deallocate(e_p)
    end if

    if (allocated(e_t)) then
        deallocate(e_t)
    end if

    if (allocated(e_po)) then
        deallocate(e_po)
    end if

    if (allocated(e_long)) then
        deallocate(e_long)
    end if

    if (allocated(t_area)) then
        deallocate(t_area)
    end if

    if (allocated(t_normal)) then
        deallocate(t_normal)
    end if

    if (allocated(e_centro)) then
        deallocate(e_centro)
    end if

    if (allocated(t_baric)) then
        deallocate(t_baric)
    end if

    if (allocated(t_baric_sub)) then
        deallocate(t_baric_sub)
    end if


end subroutine

subroutine loadscript()

    !local
    integer (kind = il) :: iunit, i
    character (len = 43) :: console
    character (len = 20) :: argcomm2, argval2
    !
    call getunit ( iunit )
    open(UNIT=iunit, FILE = trim(conf%file_script))
    print*, 'Corriendo script...'
    do i = 1, 9999999
        read(iunit,*, end=10) console
        !print*, i, '  >> ' // trim(console)
        call savetolog(inttostr(i) // '  >> ' // trim(console) // ' - - -' // gettiempo(), .true.)
        call getcommand(console, argcomm2, argval2)
        call waitcommands(.true.,argcomm2,argval2)
    end do
    
    


10 print*, 'Script finalizado'
close(iunit)
call savelog(trim(conf%test_name)// '_LOG.dat') 
end subroutine



subroutine waitcommands(justone, argcomm, argval)
    !dummy
    logical, optional, intent(in) :: justone
    character (len = 20), optional, intent(in) :: argcomm, argval
    !
    !local
    logical :: flag_script

    if (present(justone)) then
        flag_script = justone
    else
        flag_script = .false.
    end if

    if (present(argcomm)) then
    else
        flag_script = .false.
    end if
    if (present(argval)) then
    else
        flag_script = .false.
    end if

    if (flag_script == .false.) then
        print*, '* * * * * * *MLFMAxwell* * * * * * *'
        print*, ' '
        print*, 'ejecute help para ayuda'
    end if


    do
        if (flag_script) then
            comm = argcomm
            val = argval
        else
        call getcommand_console(comm, val)
        end if

        if (streq(comm,'run')) then
            if (streq(val,'mlfma')) then
                call testMLFMA()
            elseif (streq(val,'fmm')) then
                call testfmm()
            elseif (streq(val,'mom')) then
                call testMoM()
            elseif (streq(val,'wizard')) then
                call wizard
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif  (streq(comm,'help')) then
            print*, 'Debe introducir comandos con la siguiente sintaxis:'
            print*, 'comando=valor'
            print*, 'Ejecute el comando [config] para mostrar las configuraciones y las posibles simulaciones a ejecutar'
        elseif  (streq(comm,'loadscript')) then
            read(val,*) conf%file_script
            call loadscript()
        elseif  (streq(comm,'config')) then
            call showconfig()
        elseif  (streq(comm,'filename')) then
            read(val,*) conf%file_name
        elseif  (streq(comm,'testname')) then
            read(val,*) conf%test_name
        elseif  (streq(comm,'ladomax')) then
            read(val,*) conf%criterio_ladomax   
        elseif  (streq(comm,'errsolver')) then
            read(val,*) conf%errorSolverCGS
            error_solver = conf%errorSolverCGS
        elseif  (streq(comm,'rcsbi.ntheta')) then
            read(val,*) conf%numThetaRCS
        elseif  (streq(comm,'rcsbi.nphi')) then
            read(val,*) conf%numPhiRCS
        elseif  (streq(comm,'force_multipoles')) then
            read(val,*) conf%Lforz
        elseif  (streq(comm,'rcsbi.dirtheta')) then
            read(val,*) conf%dirTheta
            dir_onda(1) = sin(conf%dirTheta*pi/180.)*cos(conf%dirPhi*pi/180.)
            dir_onda(2) = sin(conf%dirTheta*pi/180.)*sin(conf%dirPhi*pi/180.)
            dir_onda(3) = cos(conf%dirTheta*pi/180.)
            call ajustar_polariz(dir_onda,pol_onda)
        elseif  (streq(comm,'rcsbi.dirphi')) then
            read(val,*) conf%dirPhi
            dir_onda(1) = sin(conf%dirTheta*pi/180.)*cos(conf%dirPhi*pi/180.)
            dir_onda(2) = sin(conf%dirTheta*pi/180.)*sin(conf%dirPhi*pi/180.)
            dir_onda(3) = cos(conf%dirTheta*pi/180.)
            call ajustar_polariz(dir_onda,pol_onda)
        elseif  (streq(comm,'rcsmonof.fmin')) then
            read(val,*) conf%rcsMonoFMin
        elseif  (streq(comm,'rcsmonof.fmax')) then
            read(val,*) conf%rcsMonoFMax
        elseif  (streq(comm,'rcsmonof.samples')) then
            read(val,*) conf%rcsMonoSamples
        elseif  (streq(comm,'precond')) then
            if (streq(val,'no')) then
                conf%usar_precond = 'NO'
            elseif (streq(val,'ilu0')) then
                conf%usar_precond = 'ILU0'
            elseif (streq(val,'ilut')) then
                conf%usar_precond = 'ILUT'
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif (streq(comm,'ilut.fillin')) then
            read(val,*) conf%fill_in_prc
        elseif (streq(comm,'ilut.droptol')) then
            read(val,*) conf%drop_tolerance
        elseif (streq(comm,'ilut.droptol')) then
            read(val,*) conf%drop_tolerance
        elseif (streq(comm,'rcsmono.active')) then
            if (streq(val,'true') .or. streq(val,'t')) then
                conf%analisis_rcs = 'MONO'
            elseif (streq(val,'false') .or. streq(val,'f')) then
                conf%analisis_rcs = 'BI'
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif (streq(comm,'rcsmono.samples')) then
            read(val,*) conf%n_pasos_rcs_mono
        elseif (streq(comm,'polariz')) then
            if (streq(val,'VV')) then
                conf%tipo_polarizac = 'VV'
            elseif (streq(val,'HH')) then
                conf%tipo_polarizac = 'HH'
            elseif (streq(val,'VH') .or. streq(val,'HV')) then
                conf%tipo_polarizac = 'VH'
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif (streq(comm,'lambda')) then
            conf%variable_espec = 'LAMBDA'
            read(val,*) lambda
            call setlambda(lambda)
        elseif (streq(comm,'frec')) then
            conf%variable_espec = 'FREC'
            read(val,*) frequency
            call setf(frequency*1000000000_dp)
        else
            print*, 'Comando no valido (recuerde no utilizar espacios)'
        end if

        if (flag_script) then
            exit
        end if

    end do


end subroutine


subroutine showconfig()

character (len = 30) :: charerrsolver, chardroptol, charlambda, charfrec

write (charerrsolver,'(ES10.2)') conf%errorSolverCGS
write (chardroptol,'(ES10.2)') conf%drop_tolerance
write (charlambda,'(ES10.2)') getlambda()
write (charfrec,'(ES10.2)') getf()

print*, 'Configuracion actual del simulador...'
print*, ' '
        print*, '> Nombre del test (guarda results automaticamente)  [testname=' , trim(conf%test_name) , ']'
        print*, '> Archivo a cargar                                  [filename=' , trim(conf%file_name) , ']'
        print*, '> # de muestras en THETA de RCS Biestatica          [rcsbi.ntheta=' // inttostr(conf%numThetaRCS) //']'
        print*, '> # de muestras en PHI de RCS Biestatica            [rcsbi.nphi=' // inttostr(conf%numPhiRCS) //']'
        !print*, '> Error del solver CGS               [errsolver=' , conf%errorSolverCGS , ']'
        print*, '> Error del solver CGS                              [errsolver=' , trim(charerrsolver) , ']'
        print*, '> Lado maximo de los cubos mas pequeños             [ladomax=' // realtostrE(conf%criterio_ladomax) // ']'
        print*, '> Numero forzado de multipolos (-1 para automatico) [force_multipoles=' // inttostr(conf%Lforz) // ']'
        print*, '> Ang(°) THETA de incidenc de onda (rcs biestatica) [rcsbi.dirtheta=' // realtostr(conf%dirTheta) //']'
        print*, '> Ang(°) PHI de incidenc de onda (rcs biestatica)   [rcsbi.dirphi=' // realtostr(conf%dirPhi) //']'
        print*, '> Frecuencia inicial barrido RCS mono frec (MHz)    [rcsmonof.fmin=' // realtostr(conf%rcsMonoFMin) //']'
        print*, '> Frecuencia final barrido RCS mono frec (MHz)      [rcsmonof.fmax=' // realtostr(conf%rcsMonoFMax) //']'
        print*, '> # de muestras para RCS mono frec                  [rcsmonof.samples=' // inttostr(conf%rcsMonoSamples) //']'
        print*, '> Tipo precondicionador (NO/ILU0/ILUT)              [precond=' // trim(conf%usar_precond) //']'
        print*, '> Fill in de ILUT (%)                               [ilut.fillin=' // realtostr(conf%fill_in_prc) //']'
        print*, '> Drop tolerance de ILUT                            [ilut.droptol=' , trim(chardroptol) , ']'
    if (streq(conf%analisis_rcs,'MONO')) then
        print*, '> Analisis rcs monostatica activo                   [rcsmono.active=true]'
    else
        print*, '> Analisis rcs monostatica activo                   [rcsmono.active=false]'
    end if
        print*, '> # de muestras RCS monostatica (0°<=phi<180°)      [rcsmono.samples=' // inttostr(conf%n_pasos_rcs_mono) //']'
        print*, '> Tipo de polarizacion (VV/HH/VH(solo rcsmono))     [polariz=' , trim(conf%tipo_polarizac) , ']'
        print*, '> Longitud de onda de trabajo (m)                   [lambda=' , trim(charlambda) , ']'
        print*, '> Frecuencia de trabajo (GHz)                       [frec=' , trim(charfrec) , ']'
print*, '- - - - - - - - - - - - - - - - - - - -'
print*, 'Comandos de ejecucion de simulaciones:'
print*, 'run=mlfma (Resolver con MLFMA)'
print*, 'run=fmm (Resolver con FMM)'
print*, 'run=mom (Resolver con MoM)'
print*, 'run=wizard (Ejecuta el asistente de Configuracion)'
print*, '[inhabilitado] run=rcsmonomom (Barrido de frec para rcs monostatica con MoM)'
print*, 'loadscript=nombre_archivo (ejecuta un script de comandos del simulador'
print*, ' '
print*, ' '

end subroutine

subroutine wizard()


    numtest = 0

    print*, '#Test a correr:'
    print*, '[0] Default'
    print*, '[1] testComparaImatxvecMLFMA (Comparacion I MoM vs MLFMA puro)'
    print*, '[2] testComparaImatxvec (comparacion I MoM vs FMM puro)'
    print*, '[3] testMLFMA (Realizar solo MLFMA puro)'
    print*, '[4] testFMM (Realizar solo FMM puro)'
    print*, '[5] testMoM (Realizar solo MoM puro)'
    print*, '[6] testcomparaZMLFMA (Comparacion Zfar MLFMA vs Zfar MoM)'
    print*, '[7] testComparaZfar (Comparacion Zfar FMM vs Zfar MoM)'
    print*, '[8] testRCSMoM_monof (RCS monostatica barrido frecuencia)'
    read*, numtest

    do
        print*, ' '
        print*, ' '
        print*, '###Seleccione correr [1], o configure la simulacion'
        print*, '[1] OK, correr'
        print*, '   [2] # de muestras en THETA de RCS Biestatica [default = ' // inttostr(conf%numThetaRCS) //']'
        print*, '   [3] # de muestras en PHI de RCS Biestatica [default = ' // inttostr(conf%numPhiRCS) //']'
        print*, '   [4] Error del solver CGS [default = ' , conf%errorSolverCGS , ']'
        print*, '   [5] Numero forzado de multipolos (-1 para automatico) [default = ' // inttostr(conf%Lforz) // ']'
        print*, '   [6] Ang(°) THETA de incidenc de onda (rcs biestatica) [default = ' // realtostr(conf%dirTheta) //']'
        print*, '   [7] Ang(°) PHI de incidenc de onda (rcs biestatica) [default = ' // realtostr(conf%dirPhi) //']'
        print*, '   [8] Frecuencia inicial barrido RCS mono frec (MHz) [default = ' // realtostr(conf%rcsMonoFMin) //']'
        print*, '   [9] Frecuencia final barrido RCS mono frec (MHz) [default = ' // realtostr(conf%rcsMonoFMax) //']'
        print*, '   [10] # de muestras para RCS mono frec [default = ' // inttostr(conf%rcsMonoSamples) //']'
        print*, '   [11] Tipo precondicionador (NO/ILU0/ILUT) [default = ' // trim(conf%usar_precond) //']'
        print*, '   [12] Fill in de ILUT (%) [default = ' // realtostr(conf%fill_in_prc) //']'
        print*, '   [13] Drop tolerance de ILUT [default = ' , conf%drop_tolerance , ']'
        print*, '   [14] Tipo de analisis RCS (mono/bi) [default = ' , trim(conf%analisis_rcs) , ']'
        print*, '   [15] # de muestras RCS monostatica (0°<=phi<180°) [default = ' // inttostr(conf%n_pasos_rcs_mono) //']'
        print*, '   [16] Tipo de polarizacion (VV/HH) [default = ' , trim(conf%tipo_polarizac) , ']'
        print*, '   [17] Variable a especificar (LAMBDA/FREC) [default = ' , trim(conf%variable_espec) , ']'

        read*, ent
        if (ent == 1) then
            exit
        else if (ent == 2) then

            print*, 'Num de muestras en THETA al calcular RCS:'
            read*, conf%numThetaRCS
        else if (ent == 3) then
            print*, 'Num de muestras en PHI al calcular RCS:'
            read*, conf%numPhiRCS
        elseif (ent == 4) then
            print*, 'Error del solver CGS:'
            read*, conf%errorSolverCGS
        elseif (ent == 5) then
            print*, 'Numero forzado de multipolos (-1 para automatico:'
            read*, conf%Lforz
        elseif (ent == 6) then
            print*, 'Angulo THETA de incidencia de la onda (Grados °)'
            read*, conf%dirTheta
        elseif (ent == 7) then
            print*, 'Angulo PHI de incidencia de la onda (Grados °)'
            read*, conf%dirPhi
        elseif (ent == 8) then
            print*, 'Frecuencia inicial de barrido RCS'
            read*, conf%rcsMonoFMin
        elseif (ent == 9) then
            print*, 'Frecuencia final de barrido RCS'
            read*, conf%rcsMonoFMax
        elseif (ent == 10) then
            print*, 'Numero de muestras RCS barrido frecuencia'
            read*, conf%rcsMonoSamples
        elseif (ent == 11) then
            print*, 'Usar precondicionador (NO/ILU0/ILUT)'
            read*, conf%usar_precond
            if (to_upper(trim(conf%usar_precond)) == 'NO') then
                conf%usar_precond = 'NO'
            else if (to_upper(trim(conf%usar_precond)) == 'ILUT') then
                conf%usar_precond = 'ILUT'
            else
                conf%usar_precond = 'ILU0'
            end if
        elseif (ent == 12) then
            print*, 'Fill in de ILUT (%)'
            read*, conf%fill_in_prc
        elseif (ent == 13) then
            print*, 'Drop tolerance de ILUT'
            read*, conf%drop_tolerance            
        elseif (ent == 14) then
            print*, 'Tipo de analisis RCS (mono/bi)'
            read*, conf%analisis_rcs
            if (to_upper(trim(conf%analisis_rcs)) == 'MONO') then
                conf%analisis_rcs = 'MONO'
            else
                conf%analisis_rcs = 'BI'
            end if
        elseif (ent == 15) then
            print*, 'Numero de muestras para RCS monostatica'
            read*, conf%n_pasos_rcs_mono
        elseif (ent == 16) then
            print*, 'Tipo de polarizacion (VV/HH)'
            read*, conf%tipo_polarizac
            if (to_upper(trim(conf%tipo_polarizac)) == 'HH') then
                conf%tipo_polarizac = 'HH'
            else
                conf%tipo_polarizac = 'VV'
            end if
        elseif (ent == 17) then
            print*, 'Variable a especificar (LAMBDA/FREC)'
            read*, conf%variable_espec
            if (to_upper(trim(conf%variable_espec)) == 'LAMBDA') then
                conf%variable_espec = 'LAMBDA'
            else
                conf%variable_espec = 'FREC'
            end if            
        end if

    end do




    !Lforzado = conf%Lforz
    error_solver = conf%errorSolverCGS
    dir_onda(1) = sin(conf%dirTheta*pi/180.)*cos(conf%dirPhi*pi/180.)
    dir_onda(2) = sin(conf%dirTheta*pi/180.)*sin(conf%dirPhi*pi/180.)
    dir_onda(3) = cos(conf%dirTheta*pi/180.)

    call ajustar_polariz(dir_onda,pol_onda)



    print*, 'Ingrese ' // trim(conf%variable_espec) // ' (GHz o m):'

    if (to_upper(trim(conf%variable_espec)) == 'LAMBDA') then
        read*, lambda
        call setlambda(lambda)
    else
        read*, frequency
        call setf(frequency*1000000000_dp)
    end if


    if (numtest == 0) then
        call testtest()
    elseif (numtest == 1) then
    !        call testComparaImatxvecMLFMA()
    elseif (numtest == 2) then
    !        call testComparaImatxvec()
    elseif (numtest == 3) then
        call testMLFMA()
    elseif (numtest == 4) then
        call testFMM()
    elseif (numtest == 5) then
        call testMoM()
    elseif (numtest == 6) then
    !        call testcomparaZMLFMA()
    elseif (numtest == 7) then
    !        call testComparaZfar()
    elseif (numtest == 8) then
    !        call testRCSMoM_monof(conf%rcsMonoFMin*1000000, conf%rcsMonoFMax*1000000, conf%rcsMonoSamples)
    end if


end subroutine



subroutine testtest()
    complex (kind = dp), allocatable, dimension(:,:) :: Zout
    complex (kind = dp), allocatable, dimension(:) :: Zmom_sparse, a1, a2
    complex (kind = dp_near), allocatable, dimension(:) :: Zmom_csr
    integer (kind = il), allocatable, dimension(:) :: i_Zmom_sparse, j_Zmom_sparse, i_Zmom_csr, j_Zmom_csr, i1, j1, i2, j2
    integer(kind = il), allocatable, dimension(:) :: iw
    integer(kind = il) :: ierr, cont, fil, col, tama
    complex (kind = dp), allocatable, dimension(:) :: vvv, yyy, yyy2, wk
    print*, '********testtest*********'

    call cap_geom()

    call timestamp

    call MoM

    call inicializar_MLFMA(e_centro, num_e, e_t, e_po, e_long, t_area, t_baric, t_baric_sub, p_coord, t_p, e_p, num_p, num_t, t_normal)

allocate(Zout(num_e,num_e))
    Zout = cerarZ(Z)

    !tama = 52992!num_e*num_e
    tama = num_e*num_e

    allocate(Zmom_sparse(tama))
    allocate(i_Zmom_sparse(tama))
    allocate(j_Zmom_sparse(tama))
    cont = 0

    allocate(vvv(num_e))
    allocate(yyy(num_e))
    allocate(yyy2(num_e))

    do fil = 1, num_e
        vvv(fil) = fil/1.3823 + (jj)*57.382/fil
    end do

    do fil = 1, num_e
        do col = 1, num_e
            if ((abs(Zout(fil,col)) < 0.0000000001) .and. .false.) then
            else
                cont = cont  + 1        
                Zmom_sparse(cont) = Zout(fil, col)
                i_Zmom_sparse(cont) = fil
                j_Zmom_sparse(cont) = col
            end if
        end do
    end do

    yyy = matmul(Zout, vvv)

print*, 'En testtest este es el numero de cercanas: ', cont
    deallocate(Zout)
    allocate(Zmom_csr(tama))
    allocate(i_Zmom_csr(num_e + 1))
    allocate(j_Zmom_csr(tama))

    call coocsrC ( num_e, tama, Zmom_sparse, i_Zmom_sparse, j_Zmom_sparse, Zmom_csr, j_Zmom_csr, i_Zmom_csr )

    deallocate(Zmom_sparse)
    deallocate(j_Zmom_sparse)
    deallocate(i_Zmom_sparse)

    !call amuxC ( num_e, vvv, yyy2, Zmom_csr, j_Zmom_csr, i_Zmom_csr )
    !call amuxC ( num_e, vvv, yyy2, Zesparcida_MLFMA, ja_MLFMA, ia_MLFMA )

print*, 'ERROR ENTRE yyy - yyy2: **********************'
print*, 100*real(sqrt(sum((yyy-yyy2)*conjg(yyy-yyy2))))/real(sqrt(sum(sqrt((yyy)*conjg(yyy)))))




print*, 'size:'
print*, size(Zesparcida_MLFMA,1)
print*, size(Zmom_csr)
print*, sum(Zmom_csr - Zesparcida_MLFMA)


    call crear_precond(tama, Zmom_csr, j_Zmom_csr, i_Zmom_csr, ierr)
    allocate(a1(num_e*num_e))
    allocate(j1(num_e*num_e))
    allocate(i1(num_e*num_e))
    allocate(wk(num_e))
    call msrcsrC ( num_e, precond_alu, precond_jlu, a1, j1, i1, wk )


    deallocate(wk)

    call crear_precond(tama, Zesparcida_MLFMA, ja_MLFMA, ia_MLFMA, ierr)
    allocate(a2(n_esparcida_MLFMA))
    allocate(j2(n_esparcida_MLFMA))
    allocate(i2(n_esparcida_MLFMA))
    allocate(wk(num_e))
    call msrcsrC ( num_e, precond_alu, precond_jlu, a2, j2, i2, wk )

call amuxC ( num_e, vvv, yyy, a1, j1, i1 )
call amuxC ( num_e, vvv, yyy2, a2, j2, i2 )

print*, 'DESPUES DE PRECONDICIONAR ERROR ENTRE yyy(Zllena) - yyy2(ZsparesMLFMA): **********************'
print*, 100*real(sqrt(sum((yyy-yyy2)*conjg(yyy-yyy2))))/real(sqrt(sum(sqrt((yyy)*conjg(yyy)))))



call setZ_MatxVec_MoM(Z)


    call CGSold(MatxVec_MLFMA,b, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)

    deallocate(precond_alu)
    deallocate(precond_jlu)
    deallocate(precond_ju)   
end subroutine



!    subroutine test(p1, p2, p3, p4)

!        integer (kind = 4) :: p1
!        integer (kind = 4), optional :: p2
!        real (kind = 8), optional, dimension (:) :: p3, p4

!        print*, p1
!        print*, p2

!        if (present(p2) .and. present(p3) .and. present(p4)) then

!            print*, p2*p3*p4


!        end if

!    end subroutine

    !subroutine testSAAD()

    !   integer (kind = il) :: ns, nj, ni
    !   integer (kind = il), allocatable :: ja(:), ia(:)
    !   real (kind = dp), allocatable :: are(:), aim(:)



    !   !Prueba de precondicionador
    !   complex (kind = 8), dimension(200) :: alu, ao, a
    !   complex (kind = 8), dimension(200) :: wa, wu, wl
    !   integer (kind = 4), dimension(200) :: jlu, ju, jr, jwl, jwu, ir, jc, jao, iao
    !   integer (kind = 4), dimension(200) :: wia
    !   integer (kind = 4) :: ierr, nnz, nrowin
    !   complex (kind = 8), dimension(10) :: bprueba, xsol
    !   a(:) = -888888.
    !   ia(:) = -1
    !   ja(:) = -1

    !   nrowin = 10
    !   open(1,file='j.dat')
    !   read(1,*) nj   ! your first line with 22
    !   allocate( ja(nj) )  ! further on you only have 21 elements
    !   read(1,*)ja          ! so, read them in
    !   !print*, 'vector ja'
    !   !print*, ja
    !   !deallocate(j)
    !   close(1)

    !   open(1,file='i.dat')
    !   read(1,*) ni   ! your first line with 22
    !   allocate( ia(ni) )  ! further on you only have 21 elements
    !   read(1,*)ia          ! so, read them in
    !   !print*, 'vector ia'
    !   !print*, ia
    !   !deallocate(i)
    !   close(1)

    !   open(1,file='s.dat')
    !   read(1,*) ns   ! your first line with 22
    !   allocate( are(ns) )  ! further on you only have 21 elements
    !   read(1,*)are          ! so, read them in
    !   !print*, 'matriz esparcida sreal'
    !   !print*, are
    !   !deallocate(s)
    !   close(1)

    !   open(1,file='si.dat')
    !   read(1,*) ns   ! your first line with 22
    !   allocate( aim(ns) )  ! further on you only have 21 elements
    !   read(1,*)aim          ! so, read them in
    !   !print*, 'matriz esparcida simag'
    !   !print*, aim
    !   !deallocate(s)
    !   close(1)

    !   a = are + (jj)*aim

    !   !a(1:7) = (/ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7/)
    !   !ja(1:7) = (/ 1, 3, 2, 4, 1, 2, 4 /)
    !   !ia(1:5) = (/1, 3, 5, 6, 8  /)
    !   call coocsrC ( nrowin, ns, a, ia, ja, ao, jao, iao )
    !   jlu = -1
    !   alu = -88888
    !   !call ilut ( nrowin, ao, jao, iao, 0*int(0.05*nrowin), 0.0000_8, alu, jlu, ju, 140, wu, wl, jr, jwl, jwu, ierr )

    !   call ilu0C ( nrowin, ao, jao, iao, alu, jlu, ju, jwl, ierr )



    !   bprueba =                 (/ 0.8147 ,   0.9058  ,  0.1270   , 0.9134  ,  0.6324  ,  0.0975  ,  0.2785  ,  0.5469  ,  0.9575  ,  0.9649/)
    !   bprueba = bprueba + (jj)* (/ 0.1576  ,  0.9706  ,  0.9572  ,  0.4854  ,  0.8003  ,  0.1419  ,  0.4218  ,  0.9157  ,  0.7922  ,  0.9595/)
    !   print*, 'bprueba: '
    !   print*, bprueba
    !   print*, '************************************'

    !   print*, 'alu: '
    !   print*, alu(1:60)
    !   print*, '************************************'

    !   print*, 'jlu: '
    !   print*, jlu(1:60)
    !   print*, '************************************'

    !   call lusol0C ( 10, bprueba, xsol, alu(1:47), jlu(1:47), ju(1:10) )


    !   print*, '******************************'
    !   print*, 'RESPUESTA FINAL XSOL:'
    !   print*, xsol
    !   print*, '******************************'
    !   a = -88888
    !   ja = -1
    !   ia = -1

    !   call msrcsrC ( nrowin, alu, jlu, a, ja, ia, wa )

    !   call csrcooC ( nrowin, 3, 200, a, ja, ia, ns, ao, ir, jc, ierr )

    !   !print*, 'alu:'
    !   !print*, alu
    !   !print*, 'jlu:'
    !   !print*, jlu
    !   call guardarVectorR(real(ao),'ilu0_ao_real')
    !   call guardarVectorR(aimag(ao),'ilu0_ao_imag')
    !   call guardarVectorI(jc,'ilu0_jc')
    !   call guardarVectorI(ir,'ilu0_ir')
    !end subroutine

    subroutine ajustar_polariz(vdirecc, vpolariz)
        !dummy
        real (kind = dp), dimension(3), intent(inout) :: vdirecc
        complex (kind = dp), dimension(3), intent(inout) :: vpolariz
        !
        !local
        real (kind = dp), dimension(3) :: vector_vertical
        !

        if (conf%tipo_polarizac == 'VV') then
            vpolariz(1) = 0.
            vpolariz(2) = 0.
            vpolariz(3) = 1.
        elseif (conf%tipo_polarizac == 'HH') then
            vector_vertical(1) = 0.
            vector_vertical(2) = 0.
            vector_vertical(3) = 1.
            vpolariz = cross_productR(vdirecc, vector_vertical)
            vpolariz = vpolariz/sqrt(sum(vpolariz*conjg(vpolariz)))
        end if

    end subroutine

    subroutine crear_precond(n_nozero, a_near, ja_near, ia_near, status)
        complex (kind = dp_near), dimension(:) :: a_near
        integer (kind = il), dimension(:) :: ja_near, ia_near
        integer (kind = il), allocatable, dimension(:) :: iw
        integer ( kind = il ), allocatable, dimension(:) :: jr
        integer ( kind = il ), allocatable, dimension(:) :: jwl
        integer ( kind = il ), allocatable, dimension(:) :: jwu
        complex ( kind = dp ), allocatable, dimension(:) :: wl
        complex ( kind = dp ), allocatable, dimension(:) :: wu

        integer (kind = il) :: n_nozero, tamano_esperado
        integer (kind = ilo) :: status

        integer (kind = il) :: fill_in
        real (kind = dp) :: drop_tol

        if (allocated(precond_jlu)) then
            deallocate(precond_jlu)
        end if 
        if (allocated(precond_alu)) then
            deallocate(precond_alu)
        end if
        if (allocated(precond_ju)) then
            deallocate(precond_ju)
        end if
        if (allocated(iw)) then
            deallocate(iw)
        end if


        !print*, 'Creando precondicionador a partir de la Z cercana la cual tiene ' // inttostr(n_nozero) // ' elementos'
        call savetolog('Creando precondicionador a partir de la Z cercana la cual tiene ' // inttostr(n_nozero) // ' elementos',.true.)
        
        if (trim(conf%usar_precond) == 'ILU0') then    
            allocate(precond_jlu(n_nozero + 1))!allocate(precond_jlu(int(n_nozero*1.5)))
            allocate(precond_alu(n_nozero + 1))!allocate(precond_alu(int(n_nozero*1.5)))
            allocate(precond_ju(num_e))!allocate(precond_ju(int(num_e*1.5)))
            allocate(iw(num_e))!allocate(iw(int(num_e*1.5)))
            precond_jlu(:) = 0
            precond_alu(:) = 0.
            precond_ju(:) = 0
            iw(:) = 0

            call ilu0C(num_e, a_near, ja_near, ia_near, precond_alu, precond_jlu, precond_ju, iw, status)
            deallocate(iw)
            if (status == 0) then
               !print*, 'Hecho'
               call savetolog('Hecho', .true.)
            else
               !print*, 'ERROR en precondicionador, se encontro pivote nulo en indice', status
               call savetolog('ERROR en precondicionador, se encontro pivote nulo en indice'// inttostr(status))
            end if

        else if (trim(conf%usar_precond) == 'ILUT') then

            fill_in = ceiling(0.01*conf%fill_in_prc*num_e)

!print*, 'Numero de elementos de fill_in: ', fill_in

call savetolog('Numero de elementos de fill_in: ' // inttostr(fill_in),.true.)
            drop_tol = conf%drop_tolerance
            tamano_esperado = 2*num_e*fill_in + n_nozero + 1
!print*, 'Tamaño de la matriz asignado para el precondicionador :', tamano_esperado

call savetolog('Tamaño de la matriz asignado para el precondicionador :' // inttostr(tamano_esperado), .true.)
            allocate(precond_jlu(tamano_esperado))!allocate(precond_jlu(int(n_nozero*1.5)))
            allocate(precond_alu(tamano_esperado))!allocate(precond_alu(int(n_nozero*1.5)))
            allocate(precond_ju(num_e))!allocate(precond_ju(int(num_e*1.5)))
            allocate(jwl(num_e))!allocate(iw(int(num_e*1.5)))
            allocate(jwu(num_e))!allocate(iw(int(num_e*1.5)))
            allocate(jr(num_e))
            allocate(wl(num_e + 1))
            allocate(wu(num_e + 1))
            precond_jlu(:) = 0
            precond_alu(:) = 0.
            precond_ju(:) = 0
            jwl(:) = 0
            jwu(:) = 0
            jr(:) = 0
            wl(:) = 0.
            wu(:) = 0.


            call ilutC (num_e, a_near, ja_near, ia_near, fill_in, drop_tol, precond_alu, precond_jlu, precond_ju, tamano_esperado, wu, wl, jr, jwl, jwu, status )

!print *, 'TODO EL PRECONDICIONADOR: '
!print*, precond_alu(:)
            deallocate(wu)
            deallocate(wl)
            deallocate(jr)
            deallocate(jwu)
            deallocate(jwl)
!print*, precond_jlu
!print*, precond_alu
            if (status == 0) then
                print*, 'Hecho'
            else if (status > 0) then
                print*, 'ERROR en precondicionador, se encontro pivote nulo en indice', status
            else if (status == -1) then
                print*, 'ERROR en precondicionador: input matrix may be wrong. (The elimination process has generated a row in L or U whose length is >  n.)'
            else if (status == -2) then
                print*, 'ERROR en precondicionador: The matrix L overflows the array alu.'
            else if (status == -3) then
                print*, 'ERROR en precondicionador: The matrix U overflows the array alu.'
            else if (status == -4) then
                print*, 'ERROR en precondicionador: Illegal value for lfil.'
            else if (status == -5) then
                print*, 'ERROR en precondicionador: zero pivot encountered.'                
            end if

        end if
        print*, 'Memoria ocupada por el precondicionador: '
        print*, (sizeof(precond_alu) + sizeof(precond_ju) + sizeof(precond_jlu))/1000000., ' MB'
    end subroutine


    subroutine gen_rcs_monostatica(test_no, prec_activo)
        !dummy
        integer (kind = ilo), intent(in) :: test_no !3 = mlfma, 4 = fmm, 5 = mom
        integer (kind = ilo), intent(in) :: prec_activo !0 = con precondicionador, otro = sin precondicionador
        !
        !local
        integer (kind = il) :: ang
        real (kind = dp) :: angphi, dephi
        real (kind = dp), dimension(3) :: direccion, prev_direccion, deltaR
        complex (kind = dp), dimension(3) :: polariz
        complex (kind = dp), allocatable, dimension(:) :: b_mono
        real (kind = dp), allocatable, dimension(:,:) :: mat_rcs_mono
        character (len = 60) :: fname
        logical :: flag_both
        character (len = 60) :: local_pol
        complex (kind = dp), allocatable, dimension(:) :: phaseCorrection
        integer (kind = dp) :: ibase
        !

        allocate(phaseCorrection(num_e))

        local_pol = conf%tipo_polarizac
        if (streq(conf%tipo_polarizac,'VH')) then
            flag_both = .true.
            conf%tipo_polarizac = 'VV'
        else
            flag_both=.false.
        end if

        do
            allocate(b_mono(num_e))
            allocate(mat_rcs_mono(conf%n_pasos_rcs_mono,2))
            mat_rcs_mono(:,:) = 0.
            guess(:) = 0.

            dephi = 1.*pi/conf%n_pasos_rcs_mono

            prev_direccion(3) = 0.
            prev_direccion(1) = -1.
            prev_direccion(2) = 0.
            do ang = 1, conf%n_pasos_rcs_mono
                print*, ' '
                print*, '*****Calculando RCS monostatica... Polarizacion: ' // conf%tipo_polarizac
                print*, '*****Angulo de incidencia ' // inttostr(ang) // ' de ' // inttostr(conf%n_pasos_rcs_mono)
                call savetolog('*****Angulo de incidencia ' // inttostr(ang) // ' de ' // inttostr(conf%n_pasos_rcs_mono))
                angphi = (ang - 1)*dephi
                direccion(3) = 0.
                direccion(1) = -cos(angphi)
                direccion(2) = -sin(angphi)

                !Ajuste de fase
                deltaR = direccion - prev_direccion
                do ibase = 1, num_e
                    phaseCorrection(ibase) = exp( (-jj)*getkappa()*dot_product(deltaR, e_centro(:,ibase) )  )
                end do
                guess = guess*phaseCorrection
                !



                call ajustar_polariz(direccion, polariz)

                b_mono = generar_onda(direccion, ctte_onda, polariz)

                if (test_no == 3) then
                    if (prec_activo == 0) then
                        call CGSold(MatxVec_MLFMA,b_mono, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)
                    else
                        call CGSold(MatxVec_MLFMA,b_mono, guess, error_solver, num_e, I, 1)
                    end if
                elseif (test_no == 4) then
                    if (prec_activo == 0) then
                        call CGSold(MatxVec_FMM,b_mono, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)
                    else
                        call CGSold(MatxVec_FMM,b_mono, guess, error_solver, num_e, I, 1)
                    end if
                elseif (test_no == 5) then
                    if (prec_activo == 0) then
                        call CGSold(MatxVec_MoM,b_mono, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)
                    else
                        call CGSold(MatxVec_MoM,b_mono, guess, error_solver, num_e, I, 1)
                    end if
                end if

                call inicializar_RCS (e_centro, I, e_long, t_baric, e_t, num_e, num_t, polariz, ctte_onda)
                mat_rcs_mono(ang,1) = angphi*180./pi
                mat_rcs_mono(ang,2) = 10*log10(sigmaRCS(pi/2.,angphi))
                guess = I
                prev_direccion = direccion
                call guardarMatrizR(mat_rcs_mono, 'last_rcs_mono_incremental.dat')
                
                print*, ' '
            end do

            if (streq(conf%test_name,strdefault)) then
                print*, 'Introduzca el nombre del archivo a guardar con la RCS monostatica:'
                read*, fname
                call guardarMatrizR(mat_rcs_mono, trim(fname))
            else
                call guardarMatrizR(mat_rcs_mono, trim(conf%test_name) // '_RCSMONO' // trim(conf%tipo_polarizac) // '.dat')
            end if

            deallocate(mat_rcs_mono)
            deallocate(b_mono)

            if (flag_both == .false.) then
                exit
            else
                flag_both = .false.
                conf%tipo_polarizac = 'HH'
            end if

        end do

        if (prec_activo == 0) then
            deallocate(precond_alu)
            deallocate(precond_jlu)
            deallocate(precond_ju)                        
        end if        

        conf%tipo_polarizac = local_pol
        deallocate(phaseCorrection)
    end subroutine


    subroutine testMLFMA()
    	!local
    	character (len = 50) :: file_i_out, ent
        integer(kind = il) :: ierr
    	!
        ierr = -1

        call cap_geom(.true.)
        print*, '********testMLFMA*********'
        call timestamp

        guess(:) = 0.
        call inicializar_MLFMA(e_centro, num_e, e_t, e_po, e_long, t_area, t_baric, t_baric_sub, p_coord, t_p, e_p, num_p, num_t, t_normal)
        b = generar_onda(dir_onda, ctte_onda, pol_onda)

        if (trim(conf%usar_precond) /= 'NO') then
            !crear precondicionador a partir de la Z cercana (Zespacida) de MLFMA
            call crear_precond(n_esparcida_MLFMA, Zesparcida_MLFMA, ja_MLFMA, ia_MLFMA, ierr)
            !
        end if

        if (conf%analisis_rcs == 'BI') then
            if (ierr == 0) then
                call CGSold(MatxVec_MLFMA,b, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)
                deallocate(precond_alu)
                deallocate(precond_jlu)
                deallocate(precond_ju)            
            else
                call CGSold(MatxVec_MLFMA,b, guess, error_solver, num_e, I, 1)
            end if
            call ComparaRCS(I,I,conf%numThetaRCS, conf%numPhiRCS)
        elseif (conf%analisis_rcs == 'MONO') then
            if (ierr == 0) then
                call gen_rcs_monostatica(3,0)
            else
                call gen_rcs_monostatica(3,-1)
            end if
        end if




        call timestamp



        if (streq(conf%test_name,strdefault)) then
            print*, 'Desea guardar el archivo de corrientes I? (S/N)'
            read*, ent
            if (trim(ent) == 'n' .or. trim(ent) == 'N') then
            else
                print*, 'Nombre de archivo:'
                read*, file_i_out
                call guardarVectorC(I, trim(file_i_out))
            end if
        else
            call guardarVectorC(I, trim(conf%test_name)//'_IMLFMA.dat')
        end if


		print*, '***FINALIZADO***'
        call destruir_todo()
    end subroutine


    subroutine testFMM()
    	!local
    	character (len = 50) :: file_i_out, ent
        integer(kind = il) :: ierr
    	!
        ierr = -1
        call cap_geom(.true.)
        print*, '********testFMM*********'
        call timestamp

        guess(:) = 0.
        call inicializar_FMM(e_centro, num_e, e_t, e_po, e_long, t_area, t_baric, t_baric_sub, p_coord, t_p, e_p, num_p, num_t)
        b = generar_onda(dir_onda, ctte_onda, pol_onda)


        if (trim(conf%usar_precond) /= 'NO') then
            !crear precondicionador a partir de la Z cercana (Zespacida) de MLFMA
            call crear_precond(n_esparcida_FMM, Zesparcida_FMM, ja_FMM, ia_FMM, ierr)
            !
        end if

        if (conf%analisis_rcs == 'BI') then
            if (ierr == 0) then
                call CGSold(MatxVec_FMM,b, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)
                deallocate(precond_alu)
                deallocate(precond_jlu)
                deallocate(precond_ju)            
            else
                call CGSold(MatxVec_FMM,b, guess, error_solver, num_e, I, 1)
            end if
            call ComparaRCS(I,I,conf%numThetaRCS, conf%numPhiRCS)
        elseif (conf%analisis_rcs == 'MONO') then
            if (ierr == 0) then
                call gen_rcs_monostatica(4,0)
            else
                call gen_rcs_monostatica(4,-1)
            end if
        end if

        call timestamp

        

        if (streq(conf%test_name,strdefault)) then
            print*, 'Desea guardar el archivo de corrientes I? (S/N)'
            read*, ent
            if (trim(ent) == 'n' .or. trim(ent) == 'N') then
            else
                print*, 'Nombre de archivo:'
                read*, file_i_out
                call guardarVectorC(I, trim(file_i_out))
            end if
        else
            call guardarVectorC(I, trim(conf%test_name)//'_IFMM.dat')
        end if


		print*, '***FINALIZADO***'
        call destruir_todo()
    end subroutine

    subroutine testMoM()
    	!local
    	character (len = 50) :: file_i_out, ent
        complex (kind = dp), allocatable, dimension(:) :: Zmom_sparse
        complex (kind = dp_near), allocatable, dimension(:) :: Zmom_csr
        integer (kind = il), allocatable, dimension(:) :: i_Zmom_sparse, j_Zmom_sparse, i_Zmom_csr, j_Zmom_csr
        integer(kind = il), allocatable, dimension(:) :: iw
        integer(kind = il) :: ierr
        integer (kind = il) :: fil, col, cont
        
    	!
        print*, '********testMoM*********'
        ierr = -1

        call cap_geom()

        call timestamp

        call MoM

        guess(:) = 0.


        call setZ_MatxVec_MoM(Z)

        if (trim(conf%usar_precond) /= 'NO') then
            !crear precondicionador a partir de la Z cercana (Zespacida) de MLFMA

            allocate(Zmom_sparse(num_e*num_e))
            allocate(i_Zmom_sparse(num_e*num_e))
            allocate(j_Zmom_sparse(num_e*num_e))
            cont = 0
            do fil = 1, num_e
                do col = 1, num_e
                    cont = cont  + 1
                    Zmom_sparse(cont) = Z(fil, col)
                    i_Zmom_sparse(cont) = fil
                    j_Zmom_sparse(cont) = col
                end do
            end do

            allocate(Zmom_csr(num_e*num_e))
            allocate(i_Zmom_csr(num_e + 1))
            allocate(j_Zmom_csr(num_e*num_e))

            call coocsrC ( num_e, num_e*num_e, Zmom_sparse, i_Zmom_sparse, j_Zmom_sparse, Zmom_csr, j_Zmom_csr, i_Zmom_csr )
            deallocate(Zmom_sparse)
            deallocate(j_Zmom_sparse)
            deallocate(i_Zmom_sparse)

            call crear_precond(num_e*num_e, Zmom_csr, j_Zmom_csr, i_Zmom_csr, ierr)

            !
            call destruir_todo()
        end if


        

        if (conf%analisis_rcs == 'BI') then
            if (ierr == 0) then
                call CGSold(MatxVec_MoM,b, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)
                deallocate(precond_alu)
                deallocate(precond_jlu)
                deallocate(precond_ju)   
            else
                call CGSold(MatxVec_MoM,b, guess, error_solver, num_e, I, 1)
            end if
            call ComparaRCS(I,I,conf%numThetaRCS, conf%numPhiRCS)
        elseif (conf%analisis_rcs == 'MONO') then
            if (ierr == 0) then
                call gen_rcs_monostatica(5,0)
            else
                call gen_rcs_monostatica(5,-1)
            end if
        end if

        !call CGS(MatxVec_MoM, b, guess, error_solver, num_e, I,1)
        !call imprimirVector (I,num_e,'I_MoM')
        call timestamp
        

        if (streq(conf%test_name,strdefault)) then
            print*, 'Desea guardar el archivo de corrientes I? (S/N)'
            read*, ent
            if (trim(ent) == 'n' .or. trim(ent) == 'N') then
            else
                print*, 'Nombre de archivo:'
                read*, file_i_out
                call guardarVectorC(I, trim(file_i_out))
            end if
        else
            call guardarVectorC(I, trim(conf%test_name)//'_IMoM.dat')
        end if


		print*, '***FINALIZADO***'
    end subroutine




    subroutine ComparaRCS(I1, I2, numTheta, numPhi)
        !dummy
        complex (kind = dp), intent(in), dimension(:) :: I1,I2
        integer (kind = il), intent(in) :: numTheta, numPhi
        !
        !local
        real (kind = dp), allocatable, dimension(:,:) :: rcsI1, rcsI2, rcsIDif
        real (kind = dp) :: normI1, normI2, normIRest, normIDif
        character (len = 30) :: nArch
        !

        allocate(rcsI1(numTheta, numPhi))
        allocate(rcsI2(numTheta, numPhi))
        allocate(rcsIDif(numTheta, numPhi))

        print*, 'Calculando RCS de las corrientes dadas por MoM en todas las direcciones...'
        rcsI1 = RCS_bi(I1,numTheta,numPhi)
        normI1 = norm2(rcsI1)
        print*, 'Calculando RCS de las corrientes dadas por FMM/MLFMA en todas las direcciones...'
        rcsI2 = RCS_bi(I2,numtheta,numphi)
        normI2 = norm2(rcsI2)
        normIRest = norm2(rcsI1-rcsI2)
        print*, 'Calculando RCS de las corrientes de diferencia en todas las direcciones...'
        rcsIDif = RCS_bi(I1-I2,numtheta,numphi)
        normIDif = norm2(rcsIDif)


        print*, 'Error relativo porcentual de matriz RCS de FMM/MLFMA respecto a MoM:'
        print*, 100*normIRest/normI1
        print*, 'Error relativo porcentual de matriz RCS de corrientes diferencia respecto a MoM: '
        print*, 100*normIDif/normI1


        if (streq(conf%test_name,strdefault)) then
            print*, 'Desea guardar los resultados de la RCS en archivos? (S/N)'
            read*, nArch
        else
            nArch = 'S'
        end if

        if (nArch == 'N' .or. nArch == 'n') then
        else

            if (streq(conf%test_name,strdefault)) then
                print*, 'Introduzca el nombre <base> de archivo: '
                read*, nArch
            else
                nArch = trim(conf%test_name) // '_RCSBI'
            end if



            print*, 'Guardando archivos'
            call guardarMatrizR(rcsI1, trim(nArch) // '_MoM.dat')
            call guardarMatrizR(rcsI2, trim(nArch) // '_FMM-MLFMA.dat')
            call guardarMatrizR(rcsIDif, trim(nArch) // '_I-DIF.dat')

        !            call guardarMatrizRGNUPLOT(rcsI1, trim(nArch) // '_GNUPLOT_MoM', numTheta, numPhi)
        !            call guardarMatrizRGNUPLOT(rcsI2, trim(nArch) // '_GNUPLOT_FMM-MLFMA', numTheta, numPhi)
        !            call guardarMatrizRGNUPLOT(rcsIDif, trim(nArch) // '_GNUPLOT_I-DIF', numTheta, numPhi)

        end if

    end subroutine






    function RCS_bi(Corriente, numtheta, numphi) result(MRCS)
        !dummy
        complex (kind = dp), intent(in), dimension(:) :: Corriente
        integer (kind = il), intent(in) :: numtheta, numphi
        real (kind = dp), allocatable, dimension(:,:) :: MRCS
        !
        !local
        integer (kind = il) :: p,t
        real (kind = dp) :: detheta, dephi, fi, tita
        integer (kind = il) :: porc, prevPorc
        !
        prevPorc = -1.
        detheta = pi/numtheta
        dephi = 2.*pi/numphi
        allocate(MRCS(numtheta,numphi))


        call inicializar_RCS (e_centro, Corriente, e_long, t_baric, e_t, num_e, num_t, pol_onda, ctte_onda)
        call setprc(numphi)
        do p = 1, numphi
            fi = (2*p-1)*dephi/2.
            do t = 1, numtheta
                tita = (2*t-1)*detheta/2.
                MRCS(t,p) = 10*log10(sigmaRCS(tita,fi))
            end do
            call updateprc(p)
            !porc = 100*p/numphi
            !if (porc/10 /= prevPorc/10) then
            !    prevPorc = porc
            !    print*, inttostr(porc) // ' %'
            !end if
        end do


    end function





    subroutine cap_geom(skipZ)
        !dummy
        logical, optional :: skipZ
        !


        ! Capturar Geometría en base a RWG
        if (streq(conf%file_name,strdefault)) then
            print*, 'Archivo a cargar:'
            read *, conf%file_name
        else
            print*, 'Archivo a cargar: ' // conf%file_name
        end if

        print*, 'Capturando geometria...'

        call obtener_p_t(conf%file_name, p_coord, t_p, num_p, num_t,t_normal)
        call parametros_rwg (p_coord, t_p, num_t, e_p, e_long, e_t, e_po, num_e, e_centro)
        call hallar_geometria_rwg (p_coord,t_p,num_p,num_t,t_area,t_baric,t_baric_sub)

        if (present(skipZ)) then
            if (skipZ) then
            else
                allocate (Z(num_e,num_e))
            end if
        else
            allocate (Z(num_e,num_e))
        end if

        allocate (I(num_e))
        allocate(guess(num_e))
        allocate (b(num_e))


        print*, 'Geometria Capturada. ' // inttostr(num_e) // ' incognitas'
        !call frecuenciaMaxMallado(t_area, e_long, num_t, num_e)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine


    subroutine MoM()

        !Usar modulo MoM
        call inicializar_MoM(p_coord, t_baric, t_baric_sub, t_area, e_long, t_p, e_p, e_t, e_po, num_e, num_p, num_t)
        print*, 'Llenando matriz Z...'
        call llenarZ_MoM(Z,1)
        print*, 'Hecho'
        print*, 'Generando onda...'
        !call guardarMatrizR(real(Z),'Zreal')
        !call guardarMatrizR(aimag(Z),'Zimag')
        b = generar_onda(dir_onda, ctte_onda, pol_onda)
        print*, 'Hecho'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine


    function generar_onda(direccion, ctte, polarizacion) result (res)
        !subroutine generar_onda(direccion, ctte, polarizacion, b)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: direccion
        complex (kind = dp), dimension(3), intent(in) :: polarizacion
        complex (kind = dp), intent(in) :: ctte
        complex (kind = dp), dimension(num_e) :: res
        !complex (kind = dp), dimension(num_e), intent(out) :: b
        !
        !local
        integer (kind = il) :: i
        real (kind = dp), dimension(3) :: baricentro
        complex (kind = dp), dimension(3) :: ondaE, ondaH
        real (kind = dp), dimension(3) :: Udir
        complex (kind = dp) :: UNO
        !
		UNO = jj/jj
        !Udir = direccion/modulov(direccion)
        Udir = direccion/sqrt(dot_product(direccion,direccion))
        do i = 1, num_e
            baricentro = t_baric(:,e_t(1,i))
            !onda = polarizacion*ctte*exp((-jj)*(w*sqrt(mu*epsi))*producto_escalarR(Udir,baricentro))
            ondaE = polarizacion*ctte*exp((-jj)*(w*sqrt(mu*epsi))*dot_product(Udir,baricentro))
            ondaH = cross_productC(Udir*UNO, ondaE)/sqrt(mu/epsi)
            !b(i) = producto_escalar(onda,baricentro-p_coord(:,e_po(1,i)))
            res(i) = alphaCFIE*sum(ondaE*( baricentro-p_coord(:,e_po(1,i)) ) ) + (1 - alphaCFIE)*sum(cross_productC(UNO*t_normal(:,e_t(1,i)),ondaH)*( baricentro-p_coord(:,e_po(1,i)) ) )

            baricentro = t_baric(:,e_t(2,i))
            !onda = polarizacion*ctte*exp((-jj)*(w*sqrt(mu*epsi))*producto_escalarR(Udir,baricentro))
            ondaE = polarizacion*ctte*exp((-jj)*(w*sqrt(mu*epsi))*dot_product(Udir,baricentro))
            ondaH = cross_productC(Udir*UNO, ondaE)/sqrt(mu/epsi)
            !b(i) = b(i) + producto_escalar(onda,p_coord(:,e_po(2,i))-baricentro)
            res(i) = res(i) + alphaCFIE*sum(ondaE*(p_coord(:,e_po(2,i))-baricentro) ) + (1 - alphaCFIE)*sum(cross_productC(UNO*t_normal(:,e_t(2,i)),ondaH)*( baricentro-p_coord(:,e_po(2,i)) ) )
            res(i) = res(i)*0.5*e_long(i)
        end do

    !end subroutine
    end function generar_onda



end program main
