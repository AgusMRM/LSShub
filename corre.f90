program list 
!------------------------------------------------------------
integer, allocatable  :: head(:,:,:), tot(:,:,:)
integer, allocatable :: link(:)
integer np
real Lbox, Lg
integer Ng
integer bx,by,bz 
real t,x,y,z,vx,vy,vz,logM,M
!--------------------------cosas para periodicidad-----------
real, allocatable :: ll(:,:)
real tg, xg, yg, zg, vxg, vyg, vzg, logmg, mg
integer inx,iny,inz, p
integer imx,imy,imz
!--------------------distancias------------------------------
real dx, dy, dz, lmin
!------------------para funcion de correlacion---------------
!integer, allocatable :: vtrg(:,:), vtrp(:,:)
integer,allocatable :: Rb(:,:),tabla(:,:)
real distc
!integer galaxias, factorialg, factorialg2, kk
integer ibin, asd , q, qq

!------------------------------------------------------------
np=4400009                !particulas en la tabla
!-----------------defino los box-----------------------------
Lbox=1000.
Lg=200                      !editar
Ng=int(Lbox/Lg)
Lmin=.5       
bin=1			!bineado			!agregue un vtrg (vector ggalaxias)  y un bineado de .1 Mpc. 
allocate(ll(3,np))
allocate(link(np))
allocate(head(Ng,Ng,Ng))
allocate(tot(Ng,Ng,Ng))
!allocate(vtrg(1,int(Lg/bin)))
!allocate(vtrp(1,int(Lg/bin)))
asd=int(lg/bin)-int(lmin/bin)
allocate(Rb(1,asd))

!vtrp=0
!vtrg=0
tot(ng,ng,ng)=0
!--------------busco head de cada box------------------------
open(unit=10,file='galaxias.dat',status='old',action='read')
do i=1,np
 read(10,*,end=2) t,x,y,z,vx,vy,vz,logM,M
!-------------- 
 ll(1,i)=x
 ll(2,i)=y
 ll(3,i)=z
!--------¿para?
 bx=int(x/Lg)+1
 by=int(y/Lg)+1
 bz=int(z/Lg)+1
 head(bx,by,bz)=i        !-----------me queda definido un head para cada box
end do
2 rewind(unit=10)
!------------------------------------------------------------
do j=1,np
 read(10,*,end=3) t,x,y,z,vx,vy,vz,logM,M
 bx=int(x/Lg)+1
 by=int(y/Lg)+1
 bz=int(z/Lg)+1
 link(head(bx,by,bz))=j
 head(bx,by,bz)=j
 tot(bx,by,bz)=tot(bx,by,bz)+1
3 end do 
!--------------PERIODICIDAD--------------------------------
!-----------------------------------------------------------
!write(*,*) 'ingrese x,y,z'
!read(*,*) xg,yg,zg
allocate(tabla(3,int((lg-lmin)/bin)))
open(12,file='correlacion.dat',status='new')
open(11,file='9box.dat',status='new')
do qq=1,20
 q=int(rand()*4000000)
	
	xg=ll(1,q)
	yg=ll(2,q)
	zg=ll(3,q)

	bx=int(xg/lg)+1    !---ingreso punto, veo a que box pertenece
	by=int(yg/lg)+1
	bz=int(zg/lg)+1

!	
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
	        dx=abs(ll(1,p)-xg)
	        dy=abs(ll(2,p)-yg)
	        dz=abs(ll(3,p)-zg)
	        if (dx>lbox/2.) dx=lbox-abs(dx)	
		if (dy>lbox/2.) dy=lbox-abs(dy)
		if (dz>lbox/2.) dz=lbox-abs(dz)
		 distc=((dx**2.)+(dy**2.)+(dz**2.))**0.5	
		 if ( distc<lg .and. distc>Lmin) then 
		    write(11,*) ll(1,p), ll(2,p), ll(3,p), distc
	!            vtrg(1,int(distc/bin))=vtrg(1,int(distc/bin))+1   			!ahora tengo en el vector la cantidad de galaxias en cada bin
	!---------------------¿a que bin pertenece?------------------------------
		    ibin=int((distc-lmin)/bin)+1
		    Rb(1,ibin)=Rb(1,ibin)+1 	    
	        endif
	        p=link(p)
	       enddo
	  end do 
	 end do  
	end do 
	 !close(11)
	
!	do i=1,int((lg-lmin)/bin)
!	 tabla(qq,i)=Rb(1,i)	
!	 write(12,*) tabla(qq,i)
!	enddo 
	 
enddo
do i=1,int((lg-lmin)/bin)
	write(12,*) Rb(1,i)
enddo
	close(11)	
	close(12)
end program list
