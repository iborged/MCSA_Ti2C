subroutine val(coord_Ti, coord_C, suma, dev, ang, op, cal_eng)
  implicit none
  real(kind=8), intent(in) :: coord_C(864,4), coord_Ti(864,3)
  !
  real(kind=8), intent(out) :: dev, cal_eng
  integer, intent(out) :: suma, ang
  integer :: nc, nti
  !
  real(kind=8) :: crd(1296,3)
  !
  integer :: i,j,m,num(5000),nas,n,suma2(5000), ang_raw
  real(kind=8) :: dx=0.0, dy=0.0, dz=0.0, dist=0.0, sx=0.0, sy=0.0, sz=0.0
  real(kind=8) :: sum_bv_tic(5000), sum_bv_titi(5000), sum_bv_cc(5000), sum_bv_cti(5000)
  real(kind=8) :: bv=0.0, sum_bv_dev=0.0, op
  real(kind=8) :: pi, c(3), vectores(5000,4)
  !
  real(kind=8) :: alpha=0.0, norma1=0.0, norma2=0.0, normat=0.0, xx=0.0, yy=0.0, zz=0.0, pp=0.0
  !
  pi=4*atan(1.0)
  !
  nc=432
  nti=864
  !
  !suma: es la sumatoria de los carbonos al cuadrado
  !
  !
  !  open(30,file="struct")          !este file en teoria tiene una estructura completa con un arreglo en particular
  !33 FORMAT(1X,F10.5,1X,F10.5,1X,F10.5)

  do i=1,5000              !inicializando estos arreglos, es decir, metiendo ceros en cada espacio de la metiria.
     sum_bv_tic(i)=0.0
     sum_bv_titi(i)=0.0
     sum_bv_cc(i)=0.0
     sum_bv_cti(i)=0.0
     num(i)=0
     suma2(i)=0            !esto es para el parametro de orden
     vectores(i,:)=0.0
  enddo
  !
  do i=1,nti                     !esto estaba comentado
     crd(i,:)=coord_Ti(i,:)    !aca se vuelve turbio el asunto
  enddo                         !miedo, terror, espanto
  !

  m=nti
  do i=1,nti
     if (int(coord_C(i,4)).ne.0) then
        m=m+1
        do j=1,3
           crd(m,j)=coord_C(i,j)
        enddo
        !print*, (crd(m,j), j=1,3)
     endif
  enddo

  !print*, m

  c(1)=(-0.001979*13.605668)       !C cuadrados (eV/C^2) |  Ti2C - IB usando R, con un R2=0.80 que se obtuvieron de un sistema 2x2x2
  c(2)=(0.6144*13.605668)          !BV          (eV/dev) |  Ti2C - IB usando R, con un R2=0.80 que se obtuvieron de un sistema 2x2x2
  c(3)=(0.009577*13.605668)        !angles      (eV/ang) |  Ti2C - IB usando R, con un R2=0.80 que se obtuvieron de un sistema 2x2x2

  !c(1)=-0.004822*13.605668     !C cuadrados (eV/C^2) |  Ti2C - IG/DA
  !c(2)=0.5894*13.605668        !BV          (eV/dev) |  Ti2C - IG/DA
  !c(3)=0.0003095*13.605668     !angles      (eV/ang) |  Ti2C - IG/DA

  sx=25.9899272919         !lattice x
  sy=25.9899272919         !lattice y
  sz=25.9899272919         !lattice z


  !                               ###########################################################################
  !                               ###########################  carbonos al  #################################
  !                               ###########################   cuadrado    #################################
  !                               ###########################################################################

  do i=1,nti
     n=0
     do j=nti+1, nti+nc
        dx=crd(i,1)-crd(j,1)
        if (dx.gt.sx/2.0) dx=dx-sx
        if (dx.lt.-sx/2.0) dx=dx+sx

        dy=crd(i,2)-crd(j,2)
        if ( dy.gt.(sy/2.0) ) dy=dy-sy
        if ( dy.lt.(-sy/2.0)) dy=dy+sy

        dz=crd(i,3)-crd(j,3)
        if ( dz.gt.(sz/2.0) ) dz=dz-sz
        if ( dz.lt.(-sz/2.0)) dz=dz+sz

        dist=sqrt(dx**2+dy**2+dz**2)

        if (dist.lt.2.4) then
           n=n+1
        endif
     enddo
     num(i)=n                 ! estoy guardando/almacenando cada valor de n en el arreglo num(i), eso es lo que se hace aca.
  enddo

  suma=0
  do i=1,nti
     suma=suma+num(i)**2     !aca esta haciendo la sumatoria de los carbonos al cuadrado
     suma2(i)=num(i)**2
     !  write(*,*)i,num(i)**2
  enddo

  op=0.0
  op=sqrt(sum(suma2,dim=1)*(1.0/nti))  !parametro de orden, Root Mean Square (RMS)

  ! write(*,*)suma   !valor final de la sumatoria de los carbonos al cuadrado


  !                               ###########################################################################
  !                               ###########################              ##################################
  !                               ########################### BOND VALENCE ##################################
  !                               ###########################              ##################################
  !                               ###########################################################################


  !Bond valence for Ti-C -- I
  !R is the length of a bond between the two given atoms. In this code R is represented by "dist"
  !The bond valence has the property that its sum around each atom in a compound is equal to the valence (oxidation state) of that atom.
  !RTiC=1.877

  !OPEN(UNIT=5000, FILE="BV_TiC.dat", ACTION="WRITE")
  bv=0.0
  do i=1,nti                        !este bucle have BV para Ti-C, usando de 1-32 estoy usando los Ti
     !n=0
     bv=0.0
     do j=nti+1,nti+nc                    !aca estoy usando los C
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx

           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0 ) dy=dy+sy

           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz

           dist=sqrt(dx**2+dy**2+dz**2)

           if (dist.lt.2.3) then
              bv=EXP(((1.877-dist)*(1.0/0.37)))    !Ti-C bond-valence parameter
              sum_bv_tic(i)=sum_bv_tic(i)+bv
              !n=n+1
              !write(111,*)i,(j-32)
           endif
        endif
     enddo
     !write(*,*)sum_bv_tic(i),i
     !write(111,*)i,j
  enddo
  !write(*,*)


  !Bond valence for Ti-Ti -- II

  !OPEN(UNIT=1001, FILE="BV_TiTi.dat", ACTION="WRITE")
  !write(1001,*)"    BV TiTi         Ti atom"

  bv=0.0
  do i=1,nti    !este buble hace BV para Ti-Ti
     !n=0
     bv=0.0
     do j=1,nti
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx

           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0) dy=dy+sy

           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz

           dist=sqrt(dx**2+dy**2+dz**2)

           if (dist.lt.3.3) then
              bv=EXP( ((2.383-dist)*(1.0/0.37)) )    !Ti-Ti bond-valence parameter
              sum_bv_titi(i)=sum_bv_titi(i)+bv
              !n=n+1
           endif
        endif
     enddo
     !write(*,*)sum_bv_titi(i),i
  enddo
  !write(*,*)


  !BOND VALENCE FOR C-C -- III
  !RCC=1.540

  !open(unit=1001, file="bv_cc.dat", action="write")
  !write(1001,*)"    BV CC         C atom"

  bv=0.0
  do i=nti+1,nti+nc
     !n=0
     bv=0.0
     do j=nti+1, nti+nc                      !C atoms
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx

           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0) dy=dy+sy

           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz

           dist=sqrt(dx**2+dy**2+dz**2)

           if (dist.lt.3.4) then
              bv=EXP( ((1.540-dist)*(1.0/0.37)) )    !C-C bond-valence parameter
              sum_bv_cc(i)=sum_bv_cc(i)+bv
           endif
        endif
     enddo
     !write(*,*) sum_bv_cc(i),i
  enddo
  !write(*,*)

  !BV for C-Ti -- IV

  !open(unit=1001, file="bv_cti.dat", action="write")
  !write(1001,*)"    BV CTi         Ti atom"

  bv=0.0
  do i=nti+1,nti+nc             !atotis de carbono.
     !n=0
     bv=0.0
     do j=1,nti             !atotis de ti
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx

           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0) dy=dy+sy

           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz

           dist=sqrt(dx**2+dy**2+dz**2)

           if (dist.lt.2.3) then
              bv=EXP( ((1.877-dist)*(1.0/0.37)) )    !C-Ti bond-valence parameter
              sum_bv_cti(i)=sum_bv_cti(i)+bv
           endif
        endif
     enddo
     !write(*,*)sum_bv_cti(i),i
  enddo
  !write(*,*)

  sum_bv_dev=0.0
  do i=1,nti
     sum_bv_dev=sum_bv_dev + ((sum_bv_tic(i) + sum_bv_titi(i)-4)**2)             !sum of ALL BV values for Ti.
  enddo                                                                          !cual es el estado de oxidacion en el Ti2C ??

  !write(*,*)sum_bv_dev

  !do i=1,32
  !   write(*,*) sum_bv_tic(i) + sum_bv_titi(i)                                  !aca imprime la valencia de cada Ti individual
  !enddo

  !write(*,*)
  dev=0.0
  dev=sqrt((sum_bv_dev/nti))  !desviacion promedio de los atotis de Ti.
  !write(*,*)dev
  !write(*,*)sum_bv_dev,dev

  !write (*,*)suma,dev


  !                               ###########################################################################
  !                               ###########################              ##################################
  !                               ###########################  Angles 180  ##################################
  !                               ###########################              ##################################
  !                               ###########################################################################

  nas=0.0                       !el truco para el contador raro, tratar de entender.
  do i=1,nti                    !almacenando las coordenadas de Tilybdeno
     n=0
     xx=0.0
     yy=0.0
     zz=0.0
     do j=nti+1,nti+nc               !atotis de Carbono
        dx=crd(i,1)-crd(j,1)
        if (dx.gt.(sx/2.0)) dx=dx-sx
        if (dx.lt.(-sx/2.0)) dx=dx+sx

        dy=crd(i,2)-crd(j,2)
        if (dy.gt.(sy/2.0)) dy=dy-sy
        if (dy.lt.(-sy/2.0)) dy=dy+sy

        dz=crd(i,3)-crd(j,3)
        if (dz.gt.(sz/2.0)) dz=dz-sz
        if (dz.lt.(-sz/2.0)) dz=dz+sz

        dist=sqrt(dx**2+dy**2+dz**2)

        if (dist.lt.2.4) then
           n=n+1
           xx=dx !crd(i,1)-crd(j,1)         !al hacer esta resta de un punto menos otro punto, esto se convierte en un vector (imp. no olvidar).
           yy=dy !crd(i,2)-crd(j,2)         !por eso es que hago la resta de coordenada por coordenada, para obtener cada vector.
           zz=dz !crd(i,3)-crd(j,3)
           !42         FORMAT(F10.1,1X,F10.5,1X,F10.5,1X,F10.5)
           !write(*,42)i,xx,yy,zz
           nas=nas+1
           vectores(nas,:)=(/dble(i), xx, yy, zz/)   !aprender que este es el formato que se usa el / / para encerrar las cantidades deseadas.
           !write(*,*)i, n, vectores(nas,:)
           !write(*,42)vectores(nas,:)
        endif
     enddo
  enddo

  ang=0     !si quito el comentario para 'ang=0' los angulos en el annealing son iguales a cero, asi que no quitar nunca.
  ang_raw=0
  do i=1,5000
     do n=1,5000
        pp=0.0
        alpha=0.0
        norma1=0.0
        norma2=0.0
        normat=0.0
        if (int(vectores(i,1)).eq.int(vectores(n,1))) then
           pp=dot_product((vectores(i,2:4)),(vectores(n,2:4)))
           norma1=sqrt(((vectores(i,2)**2) + (vectores(i,3)**2) + (vectores(i,4)**2)))
           norma2=sqrt(((vectores(n,2)**2) + (vectores(n,3)**2) + (vectores(n,4)**2)))
           normat=(norma1*norma2)
           alpha=(180/pi)*acos((pp/normat))
           !print*,alpha
           if (alpha.gt.171.and.alpha.lt.185) then
              !print*,alpha
              ang_raw=ang_raw+1
           endif
           ang=ang_raw/2
        endif
     enddo
  enddo


  !OPEN(UNIT=5000, FILE="results.dat", ACTION="WRITE")

  cal_eng=0.0
  cal_eng=cal_eng + (c(1)*suma) + (c(2)*dev) + (c(3)*(ang) )
  !cal_eng=cal_eng + (c(1)*suma) + (c(2)*dev) + (c(3)*(ang/0.5) - 4643.7120 )

  !34 FORMAT(i4,1X,F10.5,1X,i6,1X,F10.6)
  !write(*,34)suma,dev,ang,energy
  !write(*,*)
  !write(*,*)suma, dev, ang/2, op
  !write(*,*)

  return
end subroutine val
