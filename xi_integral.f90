program list 
implicit none
!gfortran -fbacktrace -fbounds-check lll.f90 para chekear linea de error
!------------------------------------------------------------
integer, allocatable  :: head(:,:,:), tot(:,:,:), galaxybin(:)
integer, allocatable :: link(:)
integer np
real Lbox, Lg
integer Ng
integer bx,by,bz 
real t,x,y,z,vx,vy,vz,logM,M
!--------------------------cosas para periodicidad-----------
real, allocatable :: ll(:,:), masax(:), masay(:), masaz(:)
real tg, xg, yg, zg, vxg, vyg, vzg, logmg, mg
integer inx,iny,inz, p,i,j, ibin
integer imx,imy,imz, bines, nt, q
!--------------------distancias------------------------------
real dx, dy, dz, lmin, lminlog,lglog, bin, distance, distlog
real v, dpares, Npares, R1, R0, dilog, direal, pi, mass
logical, allocatable ::check(:)
!------------------------------------------------------------
np=4400009             !particulas en la tabla
!-----------------defino los box-----------------------------

xg= 500.
yg= 500.
zg= 500.

mass= 0

pi=acos(-1.)

Lbox=1000.

Lg=360.  
Lg=lbox/int(lbox/lg)
   
Lglog=log10(Lg)                 !editar
Ng=int(Lbox/Lg)
Lmin=.5
Lminlog=log10(Lmin)  

bines=60
bin=(lglog-lminlog)/float(bines)

nt=1
    
allocate(ll(4,np))
allocate(link(np))
allocate(head(Ng,Ng,Ng))
allocate(tot(Ng,Ng,Ng))
allocate(galaxybin(bines+1))					!LE AGREGUE UN INDICE PORQUE SI PONIA LG=200 HABIA UNOS PARES QUE ME CAIAN EN UN BIN 31
allocate(check(np))

allocate(masax(bines+1))
allocate(masay(bines+1))
allocate(masaz(bines+1))

check(:)=.false.


open(unit=14,file='check.dat',status='unknown')
!tot(ng,ng,ng)=0
!--------------busco head de cada box------------------------
open(unit=10,file='galaxias.dat',status='old',action='read')
do i=1,np
 read(10,*,end=2) t,x,y,z,vx,vy,vz,logm,m
!-------------- 
 ll(1,i)=x
 ll(2,i)=y
 ll(3,i)=z

 ll(4,i)=logm
	
!--------¿para?
 bx=int(x/Lg)+1
 by=int(y/Lg)+1
 bz=int(z/Lg)+1
 
 if (bx==ng+1) bx=ng
 if (by==ng+1) by=ng
 if (bz==ng+1) bz=ng
 
 head(bx,by,bz)=i        !-----------me queda definido un head para cada box
 	
end do
2 rewind(unit=10)
!------------------------------------------------------------
do j=1,np
 read(10,*,end=3) t,x,y,z,vx,vy,vz,logm,m
 bx=int(x/Lg)+1
 by=int(y/Lg)+1
 bz=int(z/Lg)+1
 
 if (bx==ng+1) bx=ng
 if (by==ng+1) by=ng
 if (bz==ng+1) bz=ng
 
 link(head(bx,by,bz))=j
 head(bx,by,bz)=j
 tot(bx,by,bz)=tot(bx,by,bz)+1
3 end do



!stop            !v 
!--------------PERIODICIDAD--------------------------------
!-----------------------------------------------------------
!open(11,file='box.dat',status='unknown')
open(12,file='integral.dat',status='unknown')
open(33,file='sumatoria.dat',status='unknown')
!distlog=0.
!ibin=0

!loop1: do j=1,nt

! do
! q=rand()*np
!  if (check(q).eqv. .true.) cycle 
!  exit 
! enddo 
 
! check(q)= .true.
	
! xg=ll(1,q)
! yg=ll(2,q)
! zg=ll(3,q)
 !xg=0.
 !yg=0.
 !zg=0.
 bx=int(xg/lg)+1    !---ingreso punto, veo a que box pertenece
 by=int(yg/lg)+1
 bz=int(zg/lg)+1


 
 do inx=bx-1,bx+1
  if (inx>=1 .and. inx<=ng) then
   imx=inx
  elseif (inx==0) then
   imx=ng
  elseif (inx==ng+1) then
   imx=1
  endif
 
  do iny=by-1,by+1
   if (iny>=1 .and. iny<=ng) then
    imy=iny
   elseif (iny==0) then
    imy=ng
   elseif (iny==ng+1) then
    imy=1
   endif
  
   do inz=bz-1,bz+1
    if (inz>=1 .and. inz<=ng) then
     imz=inz
    elseif (inz==0) then
     imz=ng
    elseif (inz==ng+1) then
     imz=1
    endif
        p=head(imx,imy,imz) 	
 !----------------defino distancias con sus arreglos para los extremos----
        do i=1,tot(imx,imy,imz)
         dx=(ll(1,p)-xg)
         dy=(ll(2,p)-yg)
         dz=(ll(3,p)-zg)
       
	 mass=10**(ll(4,p))

         if (bx==ng+1) bx=ng
	 if (by==ng+1) by=ng
	 if (bz==ng+1) bz=ng
       
         if (dx>lbox/2.) dx=lbox-abs(dx)	
   	 if (dy>lbox/2.) dy=lbox-abs(dy)
 	 if (dz>lbox/2.) dz=lbox-abs(dz)

 	  distance=(((dx**2.)+(dy**2.)+(dz**2.))**0.5)
	  if (distance<lg .and. distance>lmin) then

		distlog=(log10(distance))
		
		ibin=int((distlog-lminlog)/bin)+1
	        	
		galaxybin(ibin)=galaxybin(ibin)+1
		
		write(33,*) mass, dx, dy, dz, distance

			
          endif
         p=link(p)
        enddo
   end do 
  end do  
 end do
!enddo loop1

!do i=1,nt                                                                                              PREGUNTAR PORQUE SE GENERA ESTO !!!!!
! write(14,*) check(i)
!enddo

!-----------PARES TEORICOS-----------------------------------
V=1000.**3
dpares=float(nt)*float(np)/V
!print*, nt, np

do i=1,bines                    !int((lgnorm-lminnorm)/bin)
	R1=float(i)*bin+lminlog		
	R1=10.**R1
	R0=float(i-1)*bin+lminlog	
	R0=10.**R0
	Npares=dpares*(4./3.)*pi*(R1**3-R0**3)	
	!-----------DISTANCIAS---------------
	dilog=lminlog+(float(i)-.5)*bin 
	direal=10.**dilog
	
        !print*,dpares,R0,R1
							!###CAMBIE i POR 10**i por la escala							
	write(12,*) galaxybin(i), masax(i), masay(i), masaz(i)
	
enddo
 close(33)
 close(11)
 close(12)
end program list
