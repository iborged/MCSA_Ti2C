PROGRAM lia
  implicit none
  integer :: iterations, m1, m2
  real (kind=8) :: bltz, metropolis=0.0, r_numb
  real (kind=8) :: energy1, energy2, delta
  !Subroutine
  real(kind=8) :: cal_eng, dev, op
  integer :: suma, imax, temp, ang,temp1,temp2
  !!
  integer :: i, j, c_num(432), n1, n2, flag, index_switch
  real(kind=8) :: crd_Mo(864,3), crd_C(864,3), x, a(3,3)
  real(kind=8) , dimension(:,:) :: coord_Mo(864,3),coord_C(864,4)
  !integer :: fid
  !character (len=5) :: fileno, charI
  !character (len=20):: foutname
  !!
  open(unit=86,file="tic-666.vasp")    !archivo que contiene las coordenadas de Mo
  !open(unit=87, file='initial.vasp')
  !open(unit=88, file='final.vasp')

  !write(*,*)'   energy    temp       op        dev       ang   suma'
  !25 FORMAT(1X,F10.7,1X,i6,1X,F10.4,1X,F10.6,1X,i6,1X,i6)
  !25 FORMAT(1X,i6,1X,F10.4)


  flag=0
  bltz=0.00008617333262      !(eV/K  --  electronvolt por Kelvin)
  !bltz=0.00000633339      !(Ry/K  --  electronvolt por Kelvin)

  imax=150000           !ciclos, cycles (iterations)
  temp1=1200
  temp2=1200

  read(86,*)
  read(86,*)
  do i=1,3
     read(86,*) a(i,1), a(i,2), a(i,3)
  enddo
  read(86,*)
  read(86,*)
  read(86,*)
  do i=1,864
     read(86,*)crd_Mo(i,:)
  enddo

  do i=1,864                         ! recordar la forma extrana de fortran de leer los files, ya que se pueden usar dos veces el mismo
     read(86,*)crd_C(i,:)            ! cont --> rango de lectura (de 1 a 192) y lee las primeras 192 lineas y luego "de 1 a 192" y lee las
  enddo                              ! cont --> segundas 192 lineas, es decir, de 129 a 384. | IB


  do i=1,864
     coord_C(i,1:3)=crd_C(i,1:3)     ! aca estoy guardando las coordenadas input en otro arreglo mas amplio (con 4 columnas) para ser utilizado | IB
     coord_Mo(i,1:3)=crd_Mo(i,1:3)   ! en la generacion de los arreglos. OJO que voy de 1 a 3 para tener arrays de igual dimension | IB
     coord_C(i,4)=1.0   !occupied    ! aca estoy diciendo que si en la columna cuatro hay un 1.0, eso significa que esta acupado por un carbono | IB
  enddo


!!!!!!!!!!!!!!!!!!! generating first random structure 16 C atoms vacancy and 16 occupency !!!!!!!!!!!!!!!!!!!
  do i=1,432
     c_num(i)=0
  enddo

  n2=1
  do i=1,3000
     call random_number (x)          ! esta parte de aca genera un numero random entre 1-384 | IB
     n1=int((x*864.0)) + 1           ! que se asigna a la variable n1 | IB
     !write(*,*) n1

     do j=1,n2
        if (c_num(j)==n1) then
           flag=1
        endif
     enddo

     if (flag.eq.1) then
        flag=0
        cycle
     else
        if (c_num(n2).eq.0) then
           c_num(n2)=n1
           n2=n2+1
        else
           flag=2
        endif

     endif
     if (flag.eq.2 .and. n2.eq.432) then
        exit
     endif
  enddo

  ! do i=1,192                   ! just checking that I am able to get a random mix of number between 1-384 without repetition
  !    write(*,*)c_num(i)        ! no need to uncomment this section of the code for the annealing | IB
  ! enddo


  do i=1,432                     ! aca las coordenadas que esten en la posiciones de c_num(i) van a ser cero
     coord_C(c_num(i),4)=0.0     ! cero significa que van a estar vacias, es decir, no las voy a usar en la primer estructura al azar
  enddo                          ! esto tiene sentido porque estoy "banqueando" 192 posiciones al azar al hacerlas cero en la 4ta posicion del array |IB

  !do i=1,384                           ! aca estoy imprimiendo todas las coordenadas que no sean cero en la 4ta posicion, es decir,
  !   if( coord_C(i,4).ne.0.0 ) then    ! las que se van a utilizar en la primer estructura 'input structure'
  !      write(*,*) i, coord_C(i,:)     ! IB
  !   endif
  !enddo

  !do i=1,192                    ! it prints out a list of 192 elements with 1's and 0's on the fourth column, lo cual esta muy bien
  !   write(*,*) coord_C(i,:)    ! now, which
  !enddo                         !

!!!!!!!!!! END !!!!!!!!!!


  !primer llamada de la subrutina
  call val(coord_Mo, coord_C, suma, dev, ang, op, cal_eng)

  energy1=cal_eng

  !print*,suma,dev,ang
  !write(*,*)"        energy1                   ", "energy2"

  energy2=0.0
  op=0.0
  do temp=1200,1200
     do iterations=1,imax
        energy2=0.0

        call random_number (x) !            ! 100
        n1=int(x*864.0) + 1  ! 1st index
        call random_number (x)
        n2=int(x*864.0) + 1  ! 2nd index

        if (n1==n2) cycle                   ! si n1 y n2 son iguales devuelvase a 100 y genere un nuevo par de numeros hasta q sean distintos


        m1=int(coord_C(n1,4))               ! index position of atom1
        m2=int(coord_C(n2,4))               ! index position of atom2
        if (m1==m2) cycle                   ! if occupency of atom1 and atoms2 both 0 or 1 then go back to do loop


        index_switch=int(coord_C(n1,4))     ! esta parte hace el intercambio de carbonos.
        coord_C(n1,4)=coord_C(n2,4)
        coord_C(n2,4)=index_switch


        !la segunda llamada es necesaria para calcular la energia2
        call val(coord_Mo, coord_C, suma, dev, ang, op, cal_eng)
        energy2=cal_eng

        delta=energy2-energy1
        !    =  (old - (new)
        !write(*,*) energy1, energy2
        !write(*,*)" "
        !14 format(f10.5)
        !write(*,14)delta

        if (delta.ge.0.0) then !if loop1
           !write(*,*) iterations, energy2
           energy1=energy2
           write(*,'(i6, f12.4)') temp, op
           !write(*,'(i6, f12.4)') 2200-temp, op
           cycle

        elseif (delta.lt.0.0) then
           call random_number (r_numb)
           !metropolis=exp(((delta)/((2200-temp)*bltz)))
           metropolis=exp(((delta)/(temp*bltz)))          !aca el delta y bltz estan en eV/K

           if ((r_numb-metropolis).lt.0.0) then ! if loop2
              energy2=energy1
              !cycle
              !write(*,'(i6, f12.4)') 2200-temp, op
              write(*,'(i6, f12.4)') temp, op
              !write(*,*) r_numb, metropolis
              !write(*,'(4f20.8)') delta, temp*bltz, metropolis
              !write(*,*) iterations, energy1

           else
              index_switch=int(coord_C(n1,4))    !esta es la parte nueva que inicia de la estructura previa.
              coord_C(n1,4)=coord_C(n2,4)        !aca parece estar la clave
              coord_C(n2,4)=index_switch
              !cycle
           endif ! end loop2

        endif ! end loop1


        !12 format(2f20.7)
        !write(*,12) r_numb, metropolis

        !el problema es que los deltasE son muy grandes y al dividir delta/KT obtengo numeros que al elevarlos a e^-(delta/KT) esto es cero
        !y metropolis siempre va a rechazar las estructuras.
     enddo
  enddo

  close(86)

end program lia
