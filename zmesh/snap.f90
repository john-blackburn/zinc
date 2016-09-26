! USAGE:
! subroutine getpart(part,iunit)   get next part from iunit. part of type(Tpart). (Also allocates space for advanced parts)
! subroutine wrtpart(part,iunit)   write part to iunit
! freepart(part)                   free memory for part
! inpart(r,part)                   (see below)
! snappart(r,part,flag,what,region)(see below)
! snapdb=.true. for debug info. .false. for silent
!
! User must also supply external subroutine
! subroutine snaperror(string)    write out errors
! If an error is encountered it will cause early termination of subroutine. It is the job
! of snaperror to set appropriate flags so that when sub returns (early) correct action is taken
! (normally we will terminate reading the input file)

! extrusion, turning, transition, neon
! sphere, box, cone (truncated), cylinder, ellipcyl, ellipsoid, torus, helix, trapezoid
! point, boundxup etc, line, arc, circle, rectangle, disk, plate, bubble

! Does not handle point, boundxup etc, but does read them in.

! ======================================================================

! snap functions have the form
!
!   function snapxxx(r,part,flag,what,region)
!
! They return the closest point to r on the surface (in global coords)
! part is the part of type Tpart
! flag=1: snap to surface, edge or corner
! flag=2: snap to edge or corner
! flag=3: snap to corner
!(in case of 2, 3 caller routine should check to see if distance too great)
! what: returns 1, 2, 3 if snapped to surf, edge corner resp (-1 if no point found)
! region: return region number of surf, edge, corner snapped to (currently not implemented)

! function inxxx(r,part)
! returns true if in the part

! The function will try edges and corners even if a surface is found
! This is done for safety and is not strictly necessary

!   function nearline(type,r1,r2,rc,rp,r3)
! takes line type, r1,r2,rc and returns the perpendicular point to rp as r3
! function returns true if r3 is in [r1,r2] and false otherwise
! NOTE 2d function

!   function insurf(data,rp)
! takes data block and returns true if 2D vector rp is in this surface

!   subroutine accum(rr,n,rp,type,reg)
! takes candidate point(s) rr(n,3) on type (1=surf, 2=edge, 3=corn), region reg
! and accumulates the point closest to rp
! if one of the points rr is closer to rp than the previous winner, it becomes the new winner
! Must be initialised by calling accst
! Winner returned by calling subroutine accres(r,what,region)
! If no points found, what=-1, region=-1

module snap

  implicit none

  private
  public snappart,inpart,getpart,freepart,wrtpart           ! main procedures
  public between,toupper,tolower,isfilled,is1d              ! auxilliary procedures
  public snapdb                                             ! debug flag (terminal only)
  public Tpart                                              ! type(Tpart)

  save

  type Tpart
     character(100) type
     real, pointer :: data(:,:) => null()
     character(1), pointer :: ltype(:) => null()         ! 'l' or 'a'
     character(2), pointer :: SSE(:) => null()           ! blank,'S','SE'
     real fab(5)
     integer region
     integer nsurf
     integer ncoat
     integer, pointer :: surftype(:) => null()  ! 1=surface region; 2=part; 3=region edge; 4=part edge
     integer, pointer :: surfto(:) => null()    ! either region or part no to snap to
     real, pointer :: surftol(:) => null()      ! default 0.9
     integer, pointer :: coatdata(:,:) => null()
     real rotate(3)
     character(1) :: rotorder(3)
     real shift(3)
     real xp(3),yp(3),zp(3)
     logical protect
  end type Tpart

  integer :: ifirst=1,mintype,minreg
  real minr(3),mindist

  real, parameter :: pi=3.141592654

  logical snapdb

  ! ######################################################################

contains

  function isfilled(part)
    logical isfilled
    type(Tpart) part

    character(*), parameter :: filled='extrusion,turning,transition,neon,'// &
         'sphere,box,cone,cylinder,ellipcyl,ellipsoid,torus,helix,trapezoid'

    if (index(filled,trim(part%type)).ne.0) then
       isfilled=.true.
    else
       isfilled=.false.
    endif
  end function isfilled

  function is1d(part)
    logical is1d
    type(Tpart) part

    if ( part%type.eq.'line'.or.part%type.eq.'arc'.or.part%type.eq.'rectangle'.or.&
         part%type.eq.'circle') then
       is1d=.true.
    else
       is1d=.false.
    endif
  end function is1d

  subroutine wrtpart(part,iunit)
    type (Tpart) part
    integer iunit,i,n

    write (iunit,'(a)') '=========================================================='
    write (iunit,'(2a)') 'type:  ',trim(part%type)

    if (associated(part%data)) then
       n=size(part%data,1)

       write (iunit,*) 'Part data:'

       if (part%type.eq.'extrusion') then
          write (iunit,'(a8,6a13,a8)') 'Type','x1','y1','x2','y2','xc','yc','SSE'
       else if (part%type.eq.'turning') then
          write (iunit,'(a8,6a13,a8)') 'Type','z1','rho1','z2','rho2','zc','rhoc','SSE'
       else if (part%type.eq.'transition') then
          write (iunit,'(8a13,a8)') 'XsDn','YsDn','XeDn','YeDn','XsUp','YsUp','XeUp','YeUp','SSE'
       else
          write (iunit,'(3a13)') 'x','y','z'
       endif

       if (part%type.eq.'extrusion'.or.part%type.eq.'turning') then
          do i=1,n
             if (part%ltype(i).eq.'l') then
                write (iunit,'(a8,4e13.5,26x,a8)') part%ltype(i),part%data(i,1:4),part%SSE(i)
             else
                write (iunit,'(a8,6e13.5,a8)') part%ltype(i),part%data(i,1:6),part%SSE(i)
             endif
          enddo
       else if (part%type.eq.'transition') then
          do i=1,n
             write (iunit,'(8e13.5,a8)') part%data(i,1:8),part%SSE(i)
          enddo
       else
          do i=1,n
             write (iunit,'(16e13.5)') part%data(i,1:16)
          enddo
       endif
    endif

    write (iunit,'(a,5e13.5)') 'fab:',part%fab
    write (iunit,'(a,i5)') 'region:',part%region
    write (iunit,'(a,3e13.5,3a2)') 'Rotate and order:',part%rotate,part%rotorder
    write (iunit,'(a,3e13.5)') 'Shift:',part%shift

    if (part%nsurf.gt.0) then
       write (iunit,*) 'Surface commands (1=region; 2=part; 3=region edge; 4=part edge)'
       write (iunit,'(3a5,a13)') '#','type','to','tol'
    endif

    do i=1,part%nsurf
       write (iunit,'(3i5,e13.5)') i,part%surftype(i),part%surfto(i),part%surftol(i)
    enddo

    if (part%ncoat.gt.0) then
       write (iunit,*) 'Coat commands'
       write (iunit,'(3a10)') '#','region','new reg'
    endif

    do i=1,part%ncoat
       write (iunit,'(3i10)') i,part%coatdata(i,:)
    enddo
    
    if (part%protect) then
       write (iunit,'(a)') 'PROTECTED'
    else
       write (iunit,'(a)') '(Not PROTECTED)'
    endif

  end subroutine wrtpart
    
  function between(phic,phi1,phi2)

! ----------------------------------------------------------------------
! return true if phic is between phi1 and phi2. 
! eg -2.5 (SW) between [3.14159 (W),-1.5707 (S)]
! ----------------------------------------------------------------------

    logical between
    real phic,phi1,phi2
    real phi21,phic1

    phi21=phi2-phi1
    phic1=phic-phi1

    if (phi21.gt.pi) then
       phi21=phi21-2*pi
    else if (phi21.lt.-pi) then
       phi21=phi21+2*pi
    endif

    if (phic1.gt.pi) then
       phic1=phic1-2*pi
    else if (phic1.lt.-pi) then
       phic1=phic1+2*pi
    endif

!    print *,'phi21,phic1',phi21,phic1

    if (phi21*phic1.gt.0.and.abs(phic1).lt.abs(phi21)) then
       between=.true.
    else
       between=.false.
    endif
  end function between

! ######################################################################
! r=  r(1)*x+  r(2)*y+  r(3)*z
!  = rp(1)*xp+rp(2)*yp+rp(3)*zp
    
  function topart(r,xp,yp,zp)
    real topart(3),r(3),xp(3),yp(3),zp(3)
    real, parameter :: x(3)=[1,0,0]
    real, parameter :: y(3)=[0,1,0]
    real, parameter :: z(3)=[0,0,1]

    topart(1)=r(1)*dot_product(x,xp) + r(2)*dot_product(y,xp) + r(3)*dot_product(z,xp)
    topart(2)=r(1)*dot_product(x,yp) + r(2)*dot_product(y,yp) + r(3)*dot_product(z,yp)
    topart(3)=r(1)*dot_product(x,zp) + r(2)*dot_product(y,zp) + r(3)*dot_product(z,zp)
  end function topart

  function toglob(rp,xp,yp,zp)
    real toglob(3),rp(3),xp(3),yp(3),zp(3)
    real, parameter :: x(3)=[1,0,0]
    real, parameter :: y(3)=[0,1,0]
    real, parameter :: z(3)=[0,0,1]

    toglob(1)=rp(1)*dot_product(xp,x)+rp(2)*dot_product(yp,x)+rp(3)*dot_product(zp,x)
    toglob(2)=rp(1)*dot_product(xp,y)+rp(2)*dot_product(yp,y)+rp(3)*dot_product(zp,y)
    toglob(3)=rp(1)*dot_product(xp,z)+rp(2)*dot_product(yp,z)+rp(3)*dot_product(zp,z)
  end function toglob

  ! ######################################################################

  function nearline(type,r1,r2,rc,rp,r3)
    
    ! ----------------------------------------------------------------------
    ! return true if rp projection on [r1,r2] touches the line between
    ! r1 and r2 (arc origin rc) and set r3 to be the projection
    ! This is all in 2D
    ! ----------------------------------------------------------------------

    logical nearline
    character(1) type
    real r1(2),r2(2),rc(2),rp(2),r3(2)

    real rad,ang1,ang2,angp,lambda
    
    if (type.eq.'l') then
       lambda=dot_product(rp-r1,r2-r1)/dot_product(r2-r1,r2-r1)
       
       if (lambda.gt.0.and.lambda.lt.1) then
          nearline=.true.
          r3=r1+lambda*(r2-r1)
       else
          nearline=.false.
       endif
    else
       rad=mag(r1-rc)

       ang1=atan2(r1(2)-rc(2),r1(1)-rc(1))
       ang2=atan2(r2(2)-rc(2),r2(1)-rc(1))
       angp=atan2(rp(2)-rc(2),rp(1)-rc(1))

!       print *,'r1=',r1
!       print *,'r2=',r2
!       print *,'rc=',rc
!       print *,'rp=',rp

!       print *,ang1,ang2,angp
       
       if (between(angp,ang1,ang2)) then
          nearline=.true.
          r3=rc+rad*[cos(angp),sin(angp)]
       else
          nearline=.false.
       endif
    endif
    
  end function nearline

  ! ######################################################################

  function insurf(type,data,rp)
    logical insurf
    character(1) :: type(:)
    real data(:,:),rp(2)

    integer npos,nneg,n,i
    real(8) xp,yp,det,mu,lambda,x1,y1,x2,y2,R,alpha,xc,yc
    real(8) ac,mua,mub,rad
    real phia,phib,ang1,ang2

    real, save :: sx=0.5,sy=0.5
    
    n=size(data,1)

    do

       npos=0
       nneg=0
    
       do i=1,n
          x1=data(i,1)
          y1=data(i,2)
          x2=data(i,3)
          y2=data(i,4)
          
          if (type(i).eq.'l') then
          
             xp=rp(1)
             yp=rp(2)
             
             det=sx*(y1-y2)-sy*(x1-x2)
             
             mu=((x1-xp)*(y1-y2)-(y1-yp)*(x1-x2))/det
             lambda=(sx*(y1-yp)-sy*(x1-xp))/det
             
             if (lambda.gt.0.and.lambda.lt.1) then
                if (mu.gt.0) then
                   npos=npos+1
                else
                   nneg=nneg+1
                endif
             endif
          else
             R=sqrt(1+sx**2/sy**2)
             
             alpha=atan(-sx/sy)
          
             xc=data(i,5)
             yc=data(i,6)

             rad=sqrt((x1-xc)**2+(y1-yc)**2)
             
             xp=rp(1)-xc         ! (xp,yp) with rc at origin
             yp=rp(2)-yc
          
             ac=acos((-sx/sy*yp+xp)/(R*rad))
             
             phia=ac+alpha
             phib=2*pi-ac+alpha
          
             mua=(rad*cos(phia)-xp)/sx
             mub=(rad*cos(phib)-xp)/sx
          
             ang1=atan2(y1-yc,x1-xc)
             ang2=atan2(y2-yc,x2-xc)
          
             if (between(phia,ang1,ang2)) then
                if (mua.gt.0) then
                   npos=npos+1
                else
                   nneg=nneg+1
                endif
             endif
          
             if (between(phib,ang1,ang2)) then
                if (mub.gt.0) then
                   npos=npos+1
                else
                   nneg=nneg+1
                endif
             endif
             
          endif
       enddo

       if (mod(npos,2).eq.0.and.mod(nneg,2).eq.0) then
          insurf=.false.
          return
       else if (mod(npos,2).eq.1.and.mod(nneg,2).eq.1) then
          insurf=.true.
          return
       else
!          print *,'insurf failed: trying again'
!          print *,npos,nneg
          call random_number(sx)
          call random_number(sy)
!          print *,sx,sy

          sx=sx*2-1
          sy=sy*2-1
       endif

    enddo

  end function insurf

! ######################################################################

  function mag(r)
    real r(:),mag
    mag=sqrt(dot_product(r,r))
  end function mag

  function normalise(r)
    real r(3),normalise(3)

    normalise=r/sqrt(r(1)**2+r(2)**2+r(3)**2)
  end function normalise

  function biggest(a,b,c)

    real a(3),b(3),c(3),biggest(3)
    real amag,bmag,cmag

    amag=sqrt(a(1)**2+a(2)**2+a(3)**2)
    bmag=sqrt(b(1)**2+b(2)**2+b(3)**2)
    cmag=sqrt(c(1)**2+c(2)**2+c(3)**2)

    if (amag.ge.bmag.and.amag.ge.cmag) then
       biggest=a
    else if (bmag.ge.amag.and.bmag.ge.cmag) then
       biggest=b
    else
       biggest=c
    endif

  end function biggest

  function cross(a,b)

    real a(3),b(3),cross(3)

    cross(1)=  a(2)*b(3)-b(2)*a(3)
    cross(2)=-(a(1)*b(3)-b(1)*a(3))
    cross(3)=  a(1)*b(2)-b(1)*a(2)

  end function cross

! ######################################################################

  subroutine accst
    ifirst=1
    mintype=-1
    minreg=-1
  end subroutine accst
  
  subroutine accum(rr,n,rp,type,reg)
    
    real rr(:,:),rp(3)
    integer type,reg,n,i
    
    real dist

    do i=1,n

       dist=mag(rr(i,:)-rp)
    
       if (ifirst.eq.1) then
          ifirst=0
          mindist=dist
          mintype=type
          minreg=reg
          minr=rr(i,:)
       else if (dist.lt.mindist) then
          mindist=dist
          mintype=type
          minreg=reg
          minr=rr(i,:)
       endif
    enddo
    
  end subroutine accum
  
  subroutine accres(r,what,region)

    real r(3)
    integer what, region

    r=minr
    what=mintype
    region=minreg
    
  end subroutine accres

! ######################################################################

  function getpart(part,iunit)

! ----------------------------------------------------------------------
! Read npart part from file iunit
! For now do not read individual surface commands. 
! If 1 or more "surface" command, snaplevel=1. 
! If 1 or more "surface edge" commands, snaplevel=2
! ----------------------------------------------------------------------

    type(Tpart) part
    integer iunit,getpart

    real x0(3),y0(3),z0(3),xp(3),yp(3),zp(3),angx,angy,angz
    logical inprt,indata,regdone,typedone,datadone,fabdone,shiftdone,rotdone,edge

    integer ios,nfab,dlen,ind,i,ifab,j,sind,cind,n,allocerr,ios1
    character rotstring*3,line*100,cdum*20,type*20

    real LxU,LxD,Ly,HH
    real, dimension(3) :: r1,r2,shift,zpp,xpp1,xpp2,xpp3,xpp,ypp

    getpart=1 ! error return

    x0=[1,0,0]
    y0=[0,1,0]
    z0=[0,0,1]
    
! ----------------------------------------------------------------------
! Defaults
! ----------------------------------------------------------------------

    inprt=.false.
    indata=.false.

    regdone=.false.
    typedone=.false.
    datadone=.false.
    fabdone=.false.
    rotdone=.false.
    shiftdone=.false.

    ind=0    ! data index
    sind=0   ! surface index
    cind=0   ! coat index

!    part%snaplevel=0
    part%rotate=0
    part%rotorder=['x','y','z']
    part%shift=0
    part%protect=.false.

! ----------------------------------------------------------------------
! Find no. of surface and coat commands and allocate
! ----------------------------------------------------------------------

    part%nsurf=nsurfaces('surface')
    if (part%nsurf.eq.-1) then
       call snaperror('End of file while looking for surface data')
       return
    endif

    if (part%nsurf.ne.0) then
       allocate (part%surftype(part%nsurf),part%surfto(part%nsurf),part%surftol(part%nsurf),stat=allocerr)
       if (allocerr.ne.0) then
          call snaperror('Could not allocate surface arrays')
          return
       endif
    endif

    part%ncoat=nsurfaces('coat')
    if (part%ncoat.eq.-1) then
       call snaperror('End of file while looking for coat data')
       return
    endif

    if (part%ncoat.ne.0) then
       allocate (part%coatdata(part%ncoat,2),stat=allocerr)
       if (allocerr.ne.0) then
          call snaperror('Could not allocate coat arrays')
          return
       endif
    endif

! ----------------------------------------------------------------------
! Major reading loop
! ----------------------------------------------------------------------

    do
       read (iunit,'(a)',iostat=ios) line

       if (ios.ne.0) then
          call snaperror('Unexpected end of file while processing part')
          return
       endif

       call tolower(line)

       if (len_trim(line).eq.0) cycle

       line=adjustl(line)

       if (line(1:1).eq.'*') cycle

! ----------------------------------------------------------------------
! Reading data for advanced parts
! ----------------------------------------------------------------------

       if (indata) then
          read (line,*) cdum

          if (cdum.eq.'end') then
             indata=.false.
          else
             ind=ind+1

             if (part%type.eq.'extrusion'.or.part%type.eq.'turning') then

                read (line,*,iostat=ios) part%ltype(ind)
             
                if (ios.ne.0) then
                   call snaperror('Could not find L or A in extrusion/turning data line')
                   return
                endif

                if (part%ltype(ind).eq.'l') then
                   read (line,*,iostat=ios) part%ltype(ind),(part%data(ind,j),j=1,4)
                else
                   read (line,*,iostat=ios) part%ltype(ind),(part%data(ind,j),j=1,6)
                endif

                if (ios.ne.0) then
                   call snaperror('Incomprehensible data line in extrusion/turning')
                   return
                endif

             else if (part%type.eq.'transition') then
                read (line,*,iostat=ios) (part%data(ind,j),j=1,8)
                if (ios.ne.0) then
                   call snaperror('Incomprehensible data line in transition')
                   return
                endif
             else
                read (line,*,iostat=ios) (part%data(ind,j),j=1,3)
                if (ios.ne.0) then
                   call snaperror('Incomprehensible data line in neon')
                   return
                endif
             endif

             if (index(line,'se').ne.0) then
                part%SSE(ind)='se'
             else if (index(line,'s').ne.0) then
                part%SSE(ind)='s'
             else
                part%SSE(ind)=' '
             endif
             
             datadone=.true.
          endif

! ----------------------------------------------------------------------
! Part and region
! ----------------------------------------------------------------------
          
       else if (line(1:4).eq.'part') then
          inprt=.true.
          
       else if (line(1:6).eq.'region'.and.inprt) then
          read (line,*,iostat=ios) cdum,part%region
          if (ios.ne.0) then
             call snaperror('Incomprehensible region number')
             return
          endif

          regdone=.true.

! ----------------------------------------------------------------------
! type. Set nfab. For advanced parts, allocate data storage and indata=.true.
! ----------------------------------------------------------------------

       else if (line(1:4).eq.'type'.and.inprt) then

          read (line,*,iostat=ios) cdum,part%type
          if (ios.ne.0) then
             call snaperror('Could not read part type')
             return
          endif


!          print *,part%type

          if (part%type.eq.'box') then
             nfab=3
          else if (part%type.eq.'cylinder') then
             nfab=2
          else if (part%type.eq.'sphere') then
             nfab=1
          else if (part%type.eq.'torus') then
             nfab=2
          else if (part%type.eq.'ellipcyl') then
             nfab=3
          else if (part%type.eq.'ellipsoid') then
             nfab=3
          else if (part%type.eq.'cone') then
             nfab=3
          else if (part%type.eq.'helix') then
             nfab=4
          else if (part%type.eq.'trapezoid') then
             nfab=4
          else if (part%type.eq.'point') then
             nfab=3
          else if (index(part%type,'bound').ne.0) then
             nfab=0
             fabdone=.true.
          else if (part%type.eq.'line') then
             nfab=1
          else if (part%type.eq.'arc') then
             nfab=3
          else if (part%type.eq.'circle') then
             nfab=1
          else if (part%type.eq.'disk') then
             nfab=1
          else if (part%type.eq.'rectangle') then
             nfab=2
          else if (part%type.eq.'plate') then
             nfab=2
          else if (part%type.eq.'bubble') then
             nfab=1
          else if (part%type.eq.'extrusion') then
             nfab=1
             dlen=datalen()
             if (dlen.eq.-1) then
                call snaperror('Unexpected end of file while looking for extrusion data')
                return
             endif

             allocate (part%data(dlen,6),part%ltype(dlen),part%SSE(dlen),stat=allocerr)
             if (allocerr.ne.0) then
                call snaperror('Could not allocate extrusion data')
                return
             endif

             indata=.true.
          else if (part%type.eq.'turning') then
             nfab=2       
             dlen=datalen()
             if (dlen.eq.-1) then
                call snaperror('Unexpected end of file while looking for turning data')
                return
             endif

             allocate (part%data(dlen,6),part%ltype(dlen),part%SSE(dlen),stat=allocerr)
             if (allocerr.ne.0) then
                call snaperror('Could not allocate turning data')
                return
             endif

             indata=.true.      
          else if (part%type.eq.'transition') then
             nfab=1
             dlen=datalen()
             if (dlen.eq.-1) then
                call snaperror('Unexpected end of file while looking for transition data')
                return
             endif

             allocate (part%data(dlen,8),part%ltype(dlen),part%SSE(dlen),stat=allocerr)
             if (allocerr.ne.0) then
                call snaperror('Could not allocate transition data')
                return
             endif

             indata=.true.      
          else if (part%type.eq.'neon') then
             nfab=1       
             dlen=datalen()
             if (dlen.eq.-1) then
                call snaperror('Unexpected end of file while looking for neon data')
                return
             endif

             allocate (part%data(dlen,16),part%ltype(dlen),part%SSE(dlen),stat=allocerr)
             if (allocerr.ne.0) then
                call snaperror('Could not allocate neon data')
                return
             endif

             indata=.true.      
          else
             call snaperror('unknown part type')
             return
          endif
          typedone=.true.

! ----------------------------------------------------------------------
! Fab
! ----------------------------------------------------------------------

       else if (line(1:3).eq.'fab'.and.inprt) then
          if (.not.typedone) then
             call snaperror('Type must be specified before fab')
             return
          endif

          if (nfab.gt.0) then
             read (line,*,iostat=ios) cdum,(part%fab(ifab),ifab=1,nfab)
             if (ios.ne.0) then
                call snaperror('Incomprehensible FAB line')
                return
             endif
             fabdone=.true.
          endif

! ----------------------------------------------------------------------
! Surface
! ----------------------------------------------------------------------

       else if (line(1:7).eq.'surface'.and.inprt) then
          sind=sind+1

          i=index(line,'edge')
          if (i.ne.0) then
             edge=.true.
             line=line(:i-1)//line(i+4:)
          else
             edge=.false.
          endif

          read (line,*,iostat=ios) cdum,cdum,cdum,cdum

          if (ios.ne.0) then    ! error
             read (line,*,iostat=ios1) cdum,type,part%surfto(sind)
             part%surftol(sind)=0.9
          else
             read (line,*,iostat=ios1) cdum,type,part%surfto(sind),part%surftol(sind)
          endif

          if (ios1.ne.0) then
             call snaperror('Incomprehensible SURFACE line')
             return
          endif

          if (type.eq.'region') then
             if (edge) then
                part%surftype(sind)=3
             else
                part%surftype(sind)=1
             endif
          else if (type.eq.'part') then
             if (edge) then
                part%surftype(sind)=4
             else
                part%surftype(sind)=2
             endif
          else
             call snaperror('Incomprehensible surface line')
             return
          endif             

       else if (line(1:4).eq.'coat'.and.inprt) then
          cind=cind+1

          read (line,*,iostat=ios) cdum,(part%coatdata(cind,j),j=1,2)
      
          if (ios.ne.0) then
             call snaperror('Incomprehensible COAT command')
             return
          endif

! ----------------------------------------------------------------------
! Rotate, shift
! ----------------------------------------------------------------------

       else if (line(1:6).eq.'rotate'.and.inprt) then

          if (rotdone) then
             call snaperror('Multiple rotate commands in part')
             return
          endif

          read (line,*,iostat=ios) cdum,(part%rotate(j),j=1,3),rotstring

          if (ios.ne.0) then
             read (line,*,iostat=ios1) cdum,(part%rotate(j),j=1,3)
             rotstring='xyz'
             if (ios1.ne.0) then
                call snaperror('Incomprehensible ROTATE command')
                return
             endif
          endif

          part%rotate=part%rotate*pi/180.0

          do j=1,3
             part%rotorder(j)=rotstring(j:j)
          enddo

          rotdone=.true.

       else if (line(1:6).eq.'shift'.and.inprt) then

          if (shiftdone) then
             call snaperror('Multiple shift commands in part')
             return
          endif

          read (line,*,iostat=ios) cdum,(part%shift(j),j=1,3)

          if (ios.ne.0) then
             call snaperror('Incomprehensible SHIFT command')
             return
          endif

          shiftdone=.true.

       else if (line(1:7).eq.'protect'.and.inprt) then
          part%protect=.true.

       else if (line(1:3).eq.'end'.and.inprt) then
          exit
       endif
    enddo

! ----------------------------------------------------------------------
! Sanity checks
! ----------------------------------------------------------------------

    if (.not.regdone) then
       call snaperror('region not set')
       return
    endif

    if (.not.typedone) then
       call snaperror('type not set')
       return
    endif

    if (index('extrusion,transition,turning,neon',trim(part%type)).ne.0.and..not.datadone) then
       call snaperror('Surface/extrusion/turning/neon data not supplied')
       return
    endif

    if (.not.fabdone) then
       call snaperror('Fab not supplied')
       return
    endif

! ----------------------------------------------------------------------
! set up extra data for neon
! ----------------------------------------------------------------------

    if (part%type.eq.'neon') then
       n=size(part%data,1)
       do i=1,n-1
          
          r1=part%data(i,1:3)
          r2=part%data(i+1,1:3)
          
          shift=0.5*(r1+r2)
          HH=mag(r2-r1)
          
          zpp=normalise(r2-r1)
          xpp1=cross(x0,zpp)
          xpp2=cross(y0,zpp)
          xpp3=cross(z0,zpp)
          
          xpp=biggest(xpp1,xpp2,xpp3)
          ypp=cross(zpp,xpp)
          
          part%data(i,4:6)=shift
          part%data(i,7:9)  =xpp
          part%data(i,10:12)=ypp
          part%data(i,13:15)=zpp
          part%data(i,16)=HH
       enddo
    endif

! ----------------------------------------------------------------------
! Trapezoid is set as an extrusion
! ----------------------------------------------------------------------

    if (part%type.eq.'trapezoid') then
       part%type='extrusion'
       allocate(part%data(4,4),part%ltype(4),part%SSE(4),stat=allocerr)
       if (allocerr.ne.0) then
          call snaperror('Could not allocate trapezoid data')
          return
       endif

       part%ltype='l'
       part%SSE='se'

       LxU=part%fab(1)
       LxD=part%fab(2)
       Ly=part%fab(3)
       part%fab(1)=part%fab(4)

       part%data(1,:)=[-LxD/2,-Ly/2, LxD/2,-Ly/2]
       part%data(2,:)=[ LxD/2,-Ly/2, LxU/2, Ly/2]
       part%data(3,:)=[ LxU/2, Ly/2,-LxU/2, Ly/2]
       part%data(4,:)=[-LxU/2, Ly/2,-LxD/2,-Ly/2]
    endif

! ----------------------------------------------------------------------
! Set xp,yp,zp for part
! ----------------------------------------------------------------------

    xp=x0
    yp=y0
    zp=z0

    angx=part%rotate(1)
    angy=part%rotate(2)
    angz=part%rotate(3)
    
    do i=1,3
       if (part%rotorder(i).eq.'x') then
          xp=rotx(xp,angx)    ! rotate xp about x0 angle angx
          yp=rotx(yp,angx)
          zp=rotx(zp,angx)
       else if (part%rotorder(i).eq.'y') then
          xp=roty(xp,angy)
          yp=roty(yp,angy)
          zp=roty(zp,angy)
       else
          xp=rotz(xp,angz)
          yp=rotz(yp,angz)
          zp=rotz(zp,angz)
       endif
    enddo

    part%xp=xp
    part%yp=yp
    part%zp=zp

    if (associated(part%SSE)) then
       part%SSE='se' ! ignore what the user says. Haha
    endif

    getpart=0    ! successful return

  contains

! ----------------------------------------------------------------------

    function datalen()
      integer datalen,ios,nbksp,i

      nbksp=0
      datalen=0
      do
         read (iunit,'(a)',iostat=ios) line
         call tolower(line)

         if (ios.ne.0) then
            datalen=-1
            return
         endif

         nbksp=nbksp+1

         line=adjustl(line)
         if (len_trim(line).eq.0.or.line(1:1).eq.'*') cycle

         if (line(1:3).eq.'end') then
            do i=1,nbksp
               backspace(iunit)
            enddo

            return
         endif

         datalen=datalen+1
      enddo

    end function datalen

! ----------------------------------------------------------------------

    function nsurfaces(s)
      character(*) s
      integer nsurfaces,nback,ios
      logical indata,inprt

      character(*), parameter :: advanced='extrusion,turning,neon,transition'

      nsurfaces=0
      nback=0
      indata=.false.
      inprt=.false.

      do
         read (iunit,'(a)',iostat=ios) line
         if (ios.ne.0) then
            nsurfaces=-1
            return
         endif

         nback=nback+1
         call tolower(line)
         line=adjustl(line)
         if (len_trim(line).eq.0.or.line(1:1).eq.'*') cycle
         
         if (line(1:4).eq.'part') then
            inprt=.true.
         else if (line(1:4).eq.'type'.and.inprt) then
            read (line,*) cdum,type
            if (index(advanced,trim(type)).ne.0) indata=.true.
         else if (line(1:3).eq.'end'.and.inprt) then
            if (indata) then
               indata=.false.
            else
               exit
            endif
         else
            read (line,*) cdum
            if (cdum.eq.s) nsurfaces=nsurfaces+1
         endif
      enddo

      do i=1,nback
         backspace(iunit)
      enddo
    end function nsurfaces

  end function getpart

! ######################################################################

  subroutine freepart(part)

    type (Tpart) part

    if (part%type.eq.'extrusion'.or.part%type.eq.'turning'.or.part%type.eq.'neon'&
         .or.part%type.eq.'transition') then
       deallocate(part%data,part%ltype,part%SSE)
    endif

    if (part%nsurf.ne.0) then
       deallocate(part%surftype,part%surfto,part%surftol)
    endif

    if (part%ncoat.ne.0) then
       deallocate(part%coatdata)
    endif

  end subroutine freepart

! ######################################################################
!     z
!     |          RH rotations about global x,y,z
!     |
!     +----- y
!    /
!   /
! x

  function rotx(rp,ang)
    real rotx(3),rp(3),ang,theta,rmag

    rotx=rp
    theta=atan2(rotx(3),rotx(2))
    rmag=sqrt(rotx(2)**2+rotx(3)**2)
  
    rotx(2)=rmag*cos(theta+ang)
    rotx(3)=rmag*sin(theta+ang)
  end function rotx

! ----------------------------------------------------------------------

  function roty(rp,ang)
    real roty(3),rp(3),ang,theta,rmag

    roty=rp
    theta=atan2(roty(3),roty(1))
    rmag=sqrt(roty(3)**2+roty(1)**2)
  
    roty(1)=rmag*cos(theta-ang)     ! NOTE negative
    roty(3)=rmag*sin(theta-ang)
  end function roty

! ----------------------------------------------------------------------

  function rotz(rp,ang)
    real rotz(3),rp(3),ang,theta,rmag
    
    rotz=rp
    theta=atan2(rotz(2),rotz(1))
    rmag=sqrt(rotz(2)**2+rotz(1)**2)
    
    rotz(1)=rmag*cos(theta+ang)
    rotz(2)=rmag*sin(theta+ang)
  end function rotz

! ######################################################################

  function snappart(r,part,flag,what,region)

    real r(3),snappart(3)
    type(Tpart) part
    integer flag,what,region

    if (part%type.eq.'box') then
       snappart=snapbox(r,part,flag,what,region)

    else if (part%type.eq.'cylinder') then
       snappart=snapcyl(r,part,flag,what,region)

    else if (part%type.eq.'cone') then
       snappart=snapcone(r,part,flag,what,region)

    else if (part%type.eq.'ellipcyl') then
       snappart=snapellipcyl(r,part,flag,what,region)

    else if (part%type.eq.'ellipsoid') then
       snappart=snapellipsoid(r,part,flag,what,region)

    else if (part%type.eq.'torus') then
       snappart=snaptorus(r,part,flag,what,region)

    else if (part%type.eq.'helix') then
       snappart=snaphelix(r,part,flag,what,region)

    else if (part%type.eq.'extrusion') then
       snappart=snapext(r,part,flag,what,region)

    else if (part%type.eq.'turning') then
       snappart=snapturn(r,part,flag,what,region)

    else if (part%type.eq.'transition') then
       snappart=snaptransition(r,part,flag,what,region)

    else if (part%type.eq.'neon') then
       snappart=snapneon(r,part,flag,what,region)

    else if (part%type.eq.'sphere') then
       snappart=snapsphere(r,part,flag,what,region)


! ----------------------------------------------------------------------

    else if (part%type.eq.'line') then
       snappart=snapline(r,part,flag,what,region)    ! flag .ge. 2

    else if (part%type.eq.'arc') then
       snappart=snaparc(r,part,flag,what,region)     ! flag .ge. 2
    
    else if (part%type.eq.'bubble') then
       snappart=snapbubble(r,part,flag,what,region)

    else if (part%type.eq.'circle') then
       snappart=snapdisk(r,part,2,what,region)       ! flag=2

    else if (part%type.eq.'disk') then
       snappart=snapdisk(r,part,flag,what,region)

    else if (part%type.eq.'rectangle') then
       snappart=snapplate(r,part,2,what,region)      ! flag=2

    else if (part%type.eq.'plate') then
       snappart=snapplate(r,part,flag,what,region)
    else
       call snaperror('Unknown part type for surface fitting')
    endif

  end function snappart

  function inpart(r,part)
    logical inpart
    real r(3)
    type (Tpart) part

    if (part%type.eq.'box') then
       inpart=inbox(r,part)

    else if (part%type.eq.'cylinder') then
       inpart=incyl(r,part)

    else if (part%type.eq.'cone') then
       inpart=incone(r,part)

    else if (part%type.eq.'ellipcyl') then
       inpart=inellipcyl(r,part)

    else if (part%type.eq.'ellipsoid') then
       inpart=inellipsoid(r,part)

    else if (part%type.eq.'torus') then
       inpart=intorus(r,part)

    else if (part%type.eq.'helix') then
       inpart=inhelix(r,part)

    else if (part%type.eq.'extrusion') then
       inpart=inext(r,part)

    else if (part%type.eq.'turning') then
       inpart=inturn(r,part)

    else if (part%type.eq.'transition') then
       inpart=intransition(r,part)

    else if (part%type.eq.'neon') then
       inpart=inneon(r,part)

    else if (part%type.eq.'sphere') then
       inpart=insphere(r,part)
    else
       call snaperror('Unknown filled part')
    endif
  end function inpart

! ######################################################################

  function inbox(r,part)

    real r(3),rp(3)
    type (Tpart) part
    logical inbox

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    if ( abs(rp(1)).lt.part%fab(1)/2.and. &
         abs(rp(2)).lt.part%fab(2)/2.and. &
         abs(rp(3)).lt.part%fab(3)/2) then
       inbox=.true.
    else
       inbox=.false.
    endif
  end function inbox

  function snapbox(r,part,flag,what,region)
    
    ! ----------------------------------------------------------------------
    ! returns nearest point to r on surface of box part.
    ! if flag is 1, returns absolute nearest point (surface, edge, corner)
    ! if flag is 2, returns nearest point on an edge or corner
    ! if flag is 3 returns nearest corner
    ! on return, what is object snapped to (1=surf,2=edge,3=corner)
    ! region is surface region number snapped to
    ! ----------------------------------------------------------------------
    
    real r(3),snapbox(3)
    type (Tpart) part
    integer flag,what,region

    real p(8,3),rp(3)
    integer ic,i,j
    real fab(3)

    fab=part%fab(1:3)/2

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    call accst

! ----------------------------------------------------------------------
! Get closest surface to rp
! ----------------------------------------------------------------------
    
    if (flag.eq.1) then

    do ic=1,3
       if (ic.eq.1) then
          i=2
          j=3
       else if (ic.eq.2) then
          i=1
          j=3
       else
          i=1
          j=2
       endif

       if ( rp(i).gt.-fab(i).and.rp(i).lt.fab(i).and. &
            rp(j).gt.-fab(j).and.rp(j).lt.fab(j)) then
          
          p(1,ic)=fab(ic)
          p(1,i)=rp(i)
          p(1,j)=rp(j)
          
          p(2,ic)=-fab(ic)
          p(2,i)=rp(i)
          p(2,j)=rp(j)

          call accum(p,2,rp,1,0)

       endif
    enddo

    endif

! ----------------------------------------------------------------------
! Get closest edge
! ----------------------------------------------------------------------

    if (flag.eq.1.or.flag.eq.2) then

    do ic=1,3

       if (ic.eq.1) then
          i=2
          j=3
       else if (ic.eq.2) then
          i=1
          j=3
       else
          i=1
          j=2
       endif

       if (rp(ic).gt.-fab(ic).and.rp(ic).lt.fab(ic)) then
          p(1,ic)=rp(ic)
          p(1,i)=fab(i)
          p(1,j)=fab(j)

          p(2,ic)=rp(ic)
          p(2,i)=fab(i)
          p(2,j)=-fab(j)

          p(3,ic)=rp(ic)
          p(3,i)=-fab(i)
          p(3,j)=fab(j)

          p(4,ic)=rp(ic)
          p(4,i)=-fab(i)
          p(4,j)=-fab(j)

          call accum(p,4,rp,2,0)
       endif
    enddo

    endif

! ----------------------------------------------------------------------
! get closest corner
! ----------------------------------------------------------------------

    p(1,:)=[ fab(1), fab(2), fab(3)]
    p(2,:)=[ fab(1), fab(2),-fab(3)]
    p(3,:)=[ fab(1),-fab(2), fab(3)]
    p(4,:)=[ fab(1),-fab(2),-fab(3)]
    p(5,:)=[-fab(1), fab(2), fab(3)]
    p(6,:)=[-fab(1), fab(2),-fab(3)]
    p(7,:)=[-fab(1),-fab(2), fab(3)]
    p(8,:)=[-fab(1),-fab(2),-fab(3)]

    call accum(p,8,rp,3,0)

    call accres(snapbox,what,region)
    snapbox=toglob(snapbox,part%xp,part%yp,part%zp)+part%shift

  end function snapbox

! ######################################################################

  function inext(r,part)
    logical inext
    real r(3),rp(3)
    type (Tpart) part

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    if (abs(rp(3)).lt.part%fab(1)/2.and.insurf(part%ltype,part%data,rp(1:2))) then
       inext=.true.
    else
       inext=.false.
    endif
  end function inext

  function snapext(r,part,flag,what,region)

    real r(3),snapext(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),p(size(part%data,1)*2,3),r1(2),r2(2),rc(2),r3(2)
    real top,bot
    integer ind,i,n

    character(1) ltype

! ----------------------------------------------------------------------

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    call accst

    top= part%fab(1)/2
    bot=-part%fab(1)/2
    
    n=size(part%data,1)

! ----------------------------------------------------------------------
! get closest surface
! data(1,:) x1,y1,x2,y2
! data(2,:) x2,y2,x3,y3
! data(3,:) x3,y3,x1,y1
! ----------------------------------------------------------------------

    if (snapdb) then
       print *,'Debug node'
       print *,flag,n
    endif

    if (flag.eq.1) then

    if (rp(3).gt.bot.and.rp(3).lt.top) then

       ind=0
       do i=1,n
          r1=part%data(i,1:2)
          r2=part%data(i,3:4)

          ltype=part%ltype(i)
          if (ltype.eq.'a') rc=part%data(i,5:6)

          if (snapdb) then
             print *, ltype,r1,r2,rc
          endif

          if (nearline(ltype,r1,r2,rc,rp(1:2),r3)) then
             ind=ind+1
             p(ind,1:2)=r3
             p(ind,3)=rp(3)
          endif
       enddo

       if (snapdb) then
          print *,'between top and bottom',ind
       endif

       if (ind.gt.0) then
          call accum(p,ind,rp,1,0)
       endif
    else
       if (snapdb) then
          print *,'NOT between top and bottom'
       endif
    endif

    if (insurf(part%ltype,part%data,rp(1:2))) then
       if (snapdb) then
          print *,'within endcaps'
       endif

       p(1,1:2)=rp(1:2)
       p(1,3)=top
       
       p(2,1:2)=rp(1:2)
       p(2,3)=bot
       
       call accum(p,2,rp,1,0)
    endif

    endif

! ----------------------------------------------------------------------
! get closest edge
! closest. return closest 3D vector in p to rp and dist is distance
! nearline: 2D vector between r1 and r2. ret TRUE if rp has perpendicular
! in [r1,r2] and, if true, return r3 on the line
! ----------------------------------------------------------------------

    if (flag.eq.1.or.flag.eq.2) then

    if (rp(3).gt.bot.and.rp(3).lt.top) then
       ind=0
       do i=1,n
          if (part%SSE(i).eq.'se') then
             ind=ind+1
             p(ind,3)=rp(3)
             p(ind,1:2)=part%data(i,1:2)
          endif
       enddo

       call accum(p,ind,rp,2,0)
    endif

    ind=0
    do i=1,n
       r1=part%data(i,1:2)
       r2=part%data(i,3:4)
       
       ltype=part%ltype(i)          
       if (ltype.eq.'a') rc=part%data(i,5:6)
       
       if (nearline(ltype,r1,r2,rc,rp(1:2),r3)) then
          ind=ind+1
          p(ind,1:2)=r3
          p(ind,3)=top
          
          ind=ind+1
          p(ind,1:2)=r3
          p(ind,3)=bot
       endif
    enddo
    
    if (ind.gt.0) then
       call accum(p,ind,rp,2,0)
    endif
    
    endif

! ----------------------------------------------------------------------
! Get nearest corner
! ----------------------------------------------------------------------

    ind=0
    do i=1,n
       if (part%SSE(i).eq.'se') then
          ind=ind+1
          p(ind,1:2)=part%data(i,1:2)
          p(ind,3)=top
          
          ind=ind+1
          p(ind,1:2)=part%data(i,1:2)
          p(ind,3)=bot
       endif
    enddo

    call accum(p,ind,rp,3,0)

    call accres(snapext,what,region)
    snapext=toglob(snapext,part%xp,part%yp,part%zp)+part%shift

  end function snapext

! ######################################################################

  function intransition(r,part)
    logical intransition
    real r(3),rp(3),HH2
    type (Tpart) part

    real data(size(part%data,1),4)
    character(1) :: ltype(size(part%data,1))

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    HH2=part%fab(1)/2

    if (abs(rp(3)).lt.HH2) then
       data=part%data(:,1:4)*0.5*(1-rp(3)/HH2)+part%data(:,5:8)*0.5*(1+rp(3)/HH2)
       ltype='l'

       if (insurf(ltype,data,rp(1:2))) then
          intransition=.true.
       else
          intransition=.false.
       endif
    else
       intransition=.false.
    endif

  end function intransition

  function snaptransition(r,part,flag,what,region)

    real r(3),snaptransition(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),p(size(part%data,1)*2,3),r1(2),r2(2),rc(2),r3(2)
    real top,bot
    integer ind,i,n

    character(1) :: ltype(size(part%data,1))
    real  datain(size(part%data,1),4)
    real datatop(size(part%data,1),4)
    real databot(size(part%data,1),4)

! ----------------------------------------------------------------------

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    call accst

    top= part%fab(1)/2
    bot=-part%fab(1)/2

    datain=part%data(:,1:4)*0.5*(1+rp(3)/bot)+part%data(:,5:8)*0.5*(1+rp(3)/top)

    databot=part%data(:,1:4)
    datatop=part%data(:,5:8)
    ltype='l'
    
    n=size(part%data,1)

! ----------------------------------------------------------------------
! get closest surface
! data(1,:) x1,y1,x2,y2
! data(2,:) x2,y2,x3,y3
! data(3,:) x3,y3,x1,y1
! ----------------------------------------------------------------------

    if (snapdb) then
       print *,'Debug node'
       print *,flag,n
    endif

    if (flag.eq.1) then

    if (rp(3).gt.bot.and.rp(3).lt.top) then

       ind=0
       do i=1,n
          r1=datain(i,1:2)
          r2=datain(i,3:4)

          if (snapdb) then
             print *,r1,r2,rc
          endif

          if (nearline('l',r1,r2,rc,rp(1:2),r3)) then
             ind=ind+1
             p(ind,1:2)=r3
             p(ind,3)=rp(3)
          endif
       enddo

       if (snapdb) then
          print *,'between top and bottom',ind
       endif

       if (ind.gt.0) then
          call accum(p,ind,rp,1,0)
       endif
    else
       if (snapdb) then
          print *,'NOT between top and bottom'
       endif
    endif


    if (insurf(ltype,datatop,rp(1:2))) then
       p(1,1:2)=rp(1:2)
       p(1,3)=top

       call accum(p,1,rp,1,0)
    endif

    if (insurf(ltype,databot,rp(1:2))) then
       p(1,1:2)=rp(1:2)
       p(1,3)=bot
       
       call accum(p,1,rp,1,0)
    endif

    endif

! ----------------------------------------------------------------------
! get closest edge
! closest. return closest 3D vector in p to rp and dist is distance
! nearline: 2D vector between r1 and r2. ret TRUE if rp has perpendicular
! in [r1,r2] and, if true, return r3 on the line
! ----------------------------------------------------------------------

    if (flag.eq.1.or.flag.eq.2) then

    if (rp(3).gt.bot.and.rp(3).lt.top) then
       ind=0
       do i=1,n
          if (part%SSE(i).eq.'se') then
             ind=ind+1
             p(ind,3)=rp(3)
             p(ind,1:2)=datain(i,1:2)
          endif
       enddo

       call accum(p,ind,rp,2,0)
    endif

    ind=0
    do i=1,n
       r1=datatop(i,1:2)
       r2=datatop(i,3:4)
       
       if (nearline('l',r1,r2,rc,rp(1:2),r3)) then
          ind=ind+1
          p(ind,1:2)=r3
          p(ind,3)=top
       endif
    enddo
    
    if (ind.gt.0) then
       call accum(p,ind,rp,2,0)
    endif

    ind=0
    do i=1,n
       r1=databot(i,1:2)
       r2=databot(i,3:4)
       
       if (nearline('l',r1,r2,rc,rp(1:2),r3)) then
          ind=ind+1
          p(ind,1:2)=r3
          p(ind,3)=bot
       endif
    enddo
    
    if (ind.gt.0) then
       call accum(p,ind,rp,2,0)
    endif
    
    endif

! ----------------------------------------------------------------------
! Get nearest corner
! ----------------------------------------------------------------------

    ind=0
    do i=1,n
       if (part%SSE(i).eq.'se') then
          ind=ind+1
          p(ind,1:2)=datatop(i,1:2)
          p(ind,3)=top
          
          ind=ind+1
          p(ind,1:2)=databot(i,1:2)
          p(ind,3)=bot
       endif
    enddo

    call accum(p,ind,rp,3,0)

    call accres(snaptransition,what,region)
    snaptransition=toglob(snaptransition,part%xp,part%yp,part%zp)+part%shift

  end function snaptransition

! ######################################################################
     
  function incyl(r,part)

    logical incyl
    real r(3),rp(3),rho
    type (Tpart) part

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    rho=sqrt(rp(1)**2+rp(2)**2)

    if (abs(rp(3)).lt.part%fab(2)/2.and.rho.lt.part%fab(1)) then
       incyl=.true.
    else
       incyl=.false.
    endif
  end function incyl

  function snapcyl(r,part,flag,what,region)

    real r(3),snapcyl(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),top,bot,rad,rho,p(2,3)

! ----------------------------------------------------------------------

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    top= part%fab(2)/2
    bot=-part%fab(2)/2

    rad=part%fab(1)
    rho=sqrt(rp(1)**2+rp(2)**2)

    call accst

! ----------------------------------------------------------------------
! Surface
! ----------------------------------------------------------------------

    if (flag.eq.1) then
       
       if (rp(3).gt.bot.and.rp(3).lt.top) then
          p(1,:)=[rp(1)*rad/rho,rp(2)*rad/rho,rp(3)]
          call accum(p,1,rp,1,0)
       endif

       if (rho.lt.rad) then
          p(1,:)=[rp(1),rp(2),top]
          p(2,:)=[rp(1),rp(2),bot]
          
          call accum(p,2,rp,1,0)
       endif
    endif
    
! ----------------------------------------------------------------------
! top bottom circular edge
! ----------------------------------------------------------------------

    p(1,:)=[rp(1)*rad/rho,rp(2)*rad/rho,top]
    p(2,:)=[rp(1)*rad/rho,rp(2)*rad/rho,bot]

    call accum(p,2,rp,2,0)

    call accres(snapcyl,what,region)
    snapcyl=toglob(snapcyl,part%xp,part%yp,part%zp)+part%shift

  end function snapcyl

! ######################################################################
     
  function inellipcyl(r,part)

    logical inellipcyl
    real r(3),rp(3),rho2,a2,theta,Rx,Ry,HH
    type (Tpart) part

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    rho2=rp(1)**2+rp(2)**2

    Rx=part%fab(1)
    Ry=part%fab(2)
    HH=part%fab(3)

    theta=atan2(rp(2),rp(1))

    a2=Rx**2*cos(theta)**2+Ry**2*sin(theta)**2

    if (abs(rp(3)).lt.HH/2.and.rho2.lt.a2) then
       inellipcyl=.true.
    else
       inellipcyl=.false.
    endif
  end function inellipcyl

  function snapellipcyl(r,part,flag,what,region)

    real r(3),snapellipcyl(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),top,bot,rad,rho,p(2,3),theta,Rx,Ry,HH

! ----------------------------------------------------------------------

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    Rx=part%fab(1)
    Ry=part%fab(2)
    HH=part%fab(3)

    top= HH/2
    bot=-HH/2

    rho=sqrt(rp(1)**2+rp(2)**2)
    theta=atan2(rp(2),rp(1))

    rad=sqrt(Rx**2*cos(theta)**2+Ry**2*sin(theta)**2)

    call accst

! ----------------------------------------------------------------------
! Surface
! ----------------------------------------------------------------------

    if (flag.eq.1) then
       
       if (rp(3).gt.bot.and.rp(3).lt.top) then
          p(1,:)=[rp(1)*rad/rho,rp(2)*rad/rho,rp(3)]
          call accum(p,1,rp,1,0)
       endif

       if (rho.lt.rad) then
          p(1,:)=[rp(1),rp(2),top]
          p(2,:)=[rp(1),rp(2),bot]
          
          call accum(p,2,rp,1,0)
       endif
    endif
    
! ----------------------------------------------------------------------
! top bottom circular edge
! ----------------------------------------------------------------------

    p(1,:)=[rp(1)*rad/rho,rp(2)*rad/rho,top]
    p(2,:)=[rp(1)*rad/rho,rp(2)*rad/rho,bot]

    call accum(p,2,rp,2,0)

    call accres(snapellipcyl,what,region)
    snapellipcyl=toglob(snapellipcyl,part%xp,part%yp,part%zp)+part%shift

  end function snapellipcyl

! ######################################################################

  function intorus(r,part)
    
    logical intorus
    real r(3),rp(3),rho,RR,rrad,dist
    type (Tpart) part
    
    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    
    rho=sqrt(rp(1)**2+rp(2)**2)
    
    RR=part%fab(1)
    rrad=part%fab(2)
    
    dist=sqrt((rho-RR)**2+rp(3)**2)
    
    if (dist.lt.rrad) then
       intorus=.true.
    else
       intorus=.false.
    endif
  end function intorus

! ######################################################################
  
  function snaptorus(r,part,flag,what,region)
    
    real r(3),snaptorus(3)
    type(Tpart) part
    integer flag,what,region
    
    real rp(3),rho,z,RR,rrad,offset(2),rho2,z2,dist
    
    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    rho=sqrt(rp(1)**2+rp(2)**2)
    z=rp(3)
    
    RR=part%fab(1)
    rrad=part%fab(2)
    
    offset=[rho-RR,z]
    
    dist=sqrt((rho-RR)**2+z**2)
    
    offset=offset*rrad/dist
    
    rho2=RR+offset(1)
    z2=offset(2)
    
    snaptorus=[rp(1)*rho2/rho, rp(2)*rho2/rho, z2]
    snaptorus=toglob(snaptorus,part%xp,part%yp,part%zp)+part%shift

    what=1
    region=0
  end function snaptorus

! ######################################################################

  function inellipsoid(r,part)
    
    logical inellipsoid
    real r(3),rp(3),Rx,Ry,Rz,phi,theta,rr2,rad2
    type(Tpart) part

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    Rx=part%fab(1)
    Ry=part%fab(2)
    Rz=part%fab(3)

    phi=atan2(rp(2)*Rx,rp(1)*Ry)
    theta=acos(rp(3)/Rz)

    rr2=rp(1)**2+rp(2)**2+rp(3)**2

    rad2=Rx**2*sin(theta)**2*cos(phi)**2+Ry**2*sin(theta)**2*sin(phi)**2+Rz**2*cos(theta)**2

    if (rr2.lt.rad2) then
       inellipsoid=.true.
    else
       inellipsoid=.false.
    endif
  end function inellipsoid

  function snapellipsoid(r,part,flag,what,region)

    real r(3),snapellipsoid(3)
    type(Tpart) part
    integer flag,what,region

    real rp(3),Rx,Ry,Rz,phi,theta,rr,rad

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    Rx=part%fab(1)
    Ry=part%fab(2)
    Rz=part%fab(3)

    phi=atan2(rp(2)*Rx,rp(1)*Ry)
    theta=acos(rp(3)/Rz)

    rr=sqrt(rp(1)**2+rp(2)**2+rp(3)**2)

    rad=sqrt(Rx**2*sin(theta)**2*cos(phi)**2+Ry**2*sin(theta)**2*sin(phi)**2+Rz**2*cos(theta)**2)

    snapellipsoid=rp*rad/rr
    snapellipsoid=toglob(snapellipsoid,part%xp,part%yp,part%zp)+part%shift

    what=1
    region=0

  end function snapellipsoid

! ######################################################################

  function inhelix(r,part)
    
    logical inhelix
    real r(3),rp(3),rho,RR,rrad,dist,HH,Hw,zwind,theta
    type (Tpart) part
    
    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    
    rho=sqrt(rp(1)**2+rp(2)**2)
    theta=atan2(rp(2),rp(1))
    if (theta.lt.0) theta=theta+2*pi
    
    RR=part%fab(1)
    rrad=part%fab(2)
    HH=part%fab(3)
    Hw=part%fab(4)

    zwind=Hw*theta/(2*pi)-HH/2

    inhelix=.false.
    do
       if (zwind.gt.HH/2) exit

       dist=sqrt((rho-RR)**2+(rp(3)-zwind)**2)
       if (dist.lt.rrad) then
          inhelix=.true.
          return
       endif

       zwind=zwind+Hw
    enddo

  end function inhelix
  
  function snaphelix(r,part,flag,what,region)
    
    real r(3),snaphelix(3)
    type(Tpart) part
    integer flag,what,region,ind,n,allocerr
    
    real rp(3),rho,z,RR,rrad,offset(2),rho2,z2,dist,zwind
    real HH,Hw,theta,theta2

    real, allocatable :: p(:,:)
    character(100) dum

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    call accst

    rho=sqrt(rp(1)**2+rp(2)**2)

    theta=atan2(rp(2),rp(1))
    if (theta.lt.0) theta=theta+2*pi

    z=rp(3)

    RR=part%fab(1)
    rrad=part%fab(2)
    HH=part%fab(3)
    Hw=part%fab(4)

    n=int(HH/Hw)+1
    allocate (p(2*n,3),stat=allocerr)

    if (allocerr.ne.0) then
       call snaperror('Could not allocate helix array')
       return
    endif

! ----------------------------------------------------------------------
! Winding surface
! ----------------------------------------------------------------------

    if (flag.eq.1) then
    
    zwind=Hw*theta/(2*pi)-HH/2

    ind=0
    do
       if (zwind.gt.HH/2) exit

       offset=[rho-RR,z-zwind]
    
       dist=sqrt((rho-RR)**2+(z-zwind)**2)
    
       offset=offset*rrad/dist
    
       rho2=RR +offset(1)
       z2=zwind+offset(2)
    
       ind=ind+1
       p(ind,:)=[rp(1)*rho2/rho,rp(2)*rho2/rho,z2]

       zwind=zwind+Hw
    enddo

    call accum(p,ind,rp,1,0)

! ----------------------------------------------------------------------
! End caps
! ----------------------------------------------------------------------

    theta2=2*pi*HH/Hw

    dist=sqrt((rho-RR)**2+(z-HH/2)**2)
    if (dist.lt.rrad) then
       p(1,:)=[rho*cos(theta2),rho*sin(theta2),z]
       call accum(p,1,rp,1,66)
    endif

    dist=sqrt((rho-RR)**2+(z+HH/2)**2)
    if (dist.lt.rrad) then
!       write (dum,*) rho,0.0,z,rp
!       call snaperror('End cap'//dum)
       p(1,:)=[rho,0.0,z]
       call accum(p,1,rp,1,99)
    endif

    endif

! ----------------------------------------------------------------------
! End cap edges
! ----------------------------------------------------------------------

    offset=[rho-RR,z-HH/2]
    dist=sqrt((rho-RR)**2+(z-HH/2)**2)

    offset=offset*rrad/dist
    
    rho2=RR+offset(1)
    z2=HH/2+offset(2)

    p(1,:)=[rho2*cos(theta2),rho2*sin(theta2),z2]

    call accum(p,1,rp,2,11)

! ----------------------------------------------------------------------

    offset=[rho-RR,z+HH/2]
    dist=sqrt((rho-RR)**2+(z+HH/2)**2)

    offset=offset*rrad/dist
    
    rho2=RR+offset(1)
    z2=-HH/2+offset(2)

    p(1,:)=[rho2,0.0,z2]

    call accum(p,1,rp,2,22)

    call accres(snaphelix,what,region)

!!$    if (region.eq.99) call snaperror('99 is the winner')
!!$    if (region.eq.66) call snaperror('66 is the winner')
!!$    if (region.eq.99) call snaperror('11 is the winner')
!!$    if (region.eq.66) call snaperror('22 is the winner')

    snaphelix=toglob(snaphelix,part%xp,part%yp,part%zp)+part%shift

    deallocate (p)

  end function snaphelix
    
! ######################################################################

  function inturn(r,part)

    logical inturn
    type (Tpart) part
    real r(3),rp(3),phi1,phi2,rho,phi,z

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    phi1=part%fab(1)*pi/180
    phi2=part%fab(2)*pi/180

    rho=sqrt(rp(1)**2+rp(2)**2)
    phi=atan2(rp(2),rp(1))           ! returns in [-pi,pi]
    z=rp(3)

    if (phi.lt.0) phi=phi+2*pi       ! in [0,2*pi]

    if (phi.gt.phi1.and.phi.lt.phi2.and.insurf(part%ltype,part%data,[z,rho])) then
       inturn=.true.
    else
       inturn=.false.
    endif
  end function inturn

! ######################################################################

  function snapturn(r,part,flag,what,region)

! ----------------------------------------------------------------------
! Turing part. Need to deal with case where 360 deg turning is used
! in this case, there are no end faces and no edges
! ----------------------------------------------------------------------

    real r(3),snapturn(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),phi1,phi2,phi,r1(2),r2(2),rc(2),rho1,z1,rho,z,r3(2)
    integer ind,i,n
    real p(2*size(part%data,1),3)
    character(1) ltype

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    phi1=part%fab(1)*pi/180
    phi2=part%fab(2)*pi/180

    rho=sqrt(rp(1)**2+rp(2)**2)
    phi=atan2(rp(2),rp(1))
    z=rp(3)

    if (phi.lt.0) phi=phi+2*pi       ! in [0,2*pi]

    call accst

    n=size(part%data,1)

! ----------------------------------------------------------------------
! Surfaces
! ----------------------------------------------------------------------

    if (flag.eq.1) then

! if phi=[phi1,phi2]
       if (phi.gt.phi1.and.phi.lt.phi2) then

          ind=0
          do i=1,n
             r1=part%data(i,1:2)
             r2=part%data(i,3:4)
             
             ltype=part%ltype(i)
             if (ltype.eq.'a') rc=part%data(i,5:6)
             
             if (nearline(ltype,r1,r2,rc,[z,rho],r3)) then
                ind=ind+1
                p(ind,1:2)=r3(2)*[cos(phi),sin(phi)]
                p(ind,3)=r3(1)
             endif
          enddo
       
          if (ind.gt.0) then
             call accum(p,ind,rp,1,0)
          endif
       endif

       if (insurf(part%ltype,part%data,[z,rho])) then
          p(1,1:2)=rho*[cos(phi1),sin(phi1)]
          p(1,3)=z
          
          p(2,1:2)=rho*[cos(phi2),sin(phi2)]
          p(2,3)=z
          
          call accum(p,2,rp,1,0)
       endif
       
    endif

! ----------------------------------------------------------------------
! Edges
! ----------------------------------------------------------------------

    if (flag.eq.1.or.flag.eq.2) then
       
! if phi=[phi1,phi2]
       if (phi.gt.phi1.and.phi.lt.phi2) then
          do i=1,n
             z1=part%data(i,1)
             rho1=part%data(i,2)
             
             p(i,:)=[rho1*cos(phi),rho1*sin(phi),z1]
          enddo
          
          call accum(p,n,rp,2,0)
       endif

       ind=0
       do i=1,n
             
          r1=part%data(i,1:2)
          r2=part%data(i,3:4)
          ltype=part%ltype(i)          
          if (ltype.eq.'a') rc=part%data(i,5:6)
          
          if (nearline(ltype,r1,r2,rc,[z,rho],r3)) then
             ind=ind+1
             p(ind,1:2)=r3(2)*[cos(phi1),sin(phi1)]
             p(ind,3)=r3(1)
             
             ind=ind+1
             p(ind,1:2)=r3(2)*[cos(phi2),sin(phi2)]
             p(ind,3)=r3(1)
          endif
       enddo
          
       if (ind.gt.0) then
          call accum(p,ind,rp,2,0)
       endif
       
    endif

! ----------------------------------------------------------------------
! Corners
! ----------------------------------------------------------------------

    ind=0
    do i=1,n
       z1  =part%data(i,1)
       rho1=part%data(i,2)

       ind=ind+1
       p(ind,:)=[rho1*cos(phi1),rho1*sin(phi1),z1]

       ind=ind+1
       p(ind,:)=[rho1*cos(phi2),rho1*sin(phi2),z1]
    enddo

    call accum(p,ind,rp,3,0)

    call accres(snapturn,what,region)
    snapturn=toglob(snapturn,part%xp,part%yp,part%zp)+part%shift

  end function snapturn

! ######################################################################

  function insphere(r,part)

    logical insphere
    type (Tpart) part
    real r(3),rp(3)

    rp=r-part%shift

    if (rp(1)**2+rp(2)**2+rp(3)**2.lt.part%fab(1)**2) then
       insphere=.true.
    else
       insphere=.false.
    endif
  end function insphere

  function snapsphere(r,part,flag,what,region)

    real snapsphere(3),r(3),rp(3),rr
    type (Tpart) part
    integer flag,what,region

    rp=r-part%shift

    rr=sqrt(rp(1)**2+rp(2)**2+rp(3)**2)

    snapsphere=rp*part%fab(1)/rr+part%shift
    what=1
    region=0
  end function snapsphere

! ######################################################################

  function incone(r,part)

    logical incone
    real r(3),rp(3),rho,RR,HH,Hz,rad,z
    type (Tpart) part

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    rho=sqrt(rp(1)**2+rp(2)**2)
    z=rp(3)

    RR=part%fab(1)
    HH=part%fab(2)
    Hz=part%fab(3)

    rad=RR-z/HH*RR    ! radius at current height

    if (z.gt.0.and.z.lt.Hz.and.rho.lt.rad) then
       incone=.true.
    else
       incone=.false.
    endif

  end function incone

  function snapcone(r,part,flag,what,region)

    real r(3),snapcone(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),rho,p(2,3),RR,HH,Hz,z,rvec(2),lambda,radtop

! ----------------------------------------------------------------------

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    RR=part%fab(1)
    HH=part%fab(2)
    Hz=part%fab(3)

    radtop=RR*(1-Hz/HH)

    rho=sqrt(rp(1)**2+rp(2)**2)
    z=rp(3)

    call accst

! ----------------------------------------------------------------------
! Curved surface and top and bottom
! ----------------------------------------------------------------------

    if (flag.eq.1) then

       lambda=((RR-rho)*RR+z*HH)/(RR**2+HH**2)

       rvec=[RR,0.0]+lambda*[-RR,HH]

       if (rvec(2).gt.0.and.rvec(2).lt.Hz.and.rho.ne.0.0) then
          p(1,:)=[rp(1)*rvec(1)/rho, rp(2)*rvec(1)/rho, rvec(2)]
          call accum(p,1,rp,1,0)
       endif

       if (rho.lt.RR) then
          p(1,:)=[rp(1),rp(2),0.0]
          call accum(p,1,rp,1,0)
       endif

       if (rho.lt.radtop) then
          p(1,:)=[rp(1),rp(2),Hz]
          call accum(p,1,rp,1,0)
       endif
    endif
    
! ----------------------------------------------------------------------
! top bottom circular edge
! ----------------------------------------------------------------------

    p(1,:)=[rp(1)*RR/rho,rp(2)*RR/rho,0.0]
    p(2,:)=[rp(1)*radtop/rho,rp(2)*radtop/rho,Hz]

    call accum(p,2,rp,2,0)

    call accres(snapcone,what,region)
    snapcone=toglob(snapcone,part%xp,part%yp,part%zp)+part%shift

  end function snapcone

! ######################################################################

  function inneon(r,part)

    logical inneon
    real r(3),rp(3),rr2,shift(3),xpp(3),ypp(3),zpp(3),HH,rpp(3)
    real rho2,dist12,dist22
    type(Tpart) part
    integer n,i

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
 
    n=size(part%data,1)

    rr2=part%fab(1)**2

    inneon=.false.

    do i=1,n-1
       shift=part%data(i,4:6)
       xpp=part%data(i,7:9)
       ypp=part%data(i,10:12)
       zpp=part%data(i,13:15)
       HH=part%data(i,16)

       rpp=topart(rp-shift,xpp,ypp,zpp)

       rho2=rpp(1)**2+rpp(2)**2
       dist12=rpp(1)**2+rpp(2)**2+(rpp(3)-HH/2)**2
       dist22=rpp(1)**2+rpp(2)**2+(rpp(3)+HH/2)**2

       if ((abs(rpp(3)).lt.HH/2.and.rho2.lt.rr2).or.dist12.lt.rr2.or.dist22.lt.rr2) then
          inneon=.true.
          return
       endif
    enddo

  end function inneon
          
  function snapneon(r,part,flag,what,region)

    real r(3),snapneon(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),rr,xpp(3),ypp(3),zpp(3),HH,rpp(3),rho,test(3)
    integer n,i
    real p(2,3),dist,offset(3),shift(3)

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    n=size(part%data,1)

    rr=part%fab(1)

    call accst

    do i=1,n-1
       shift=part%data(i,4:6)
       xpp=part%data(i,7:9)
       ypp=part%data(i,10:12)
       zpp=part%data(i,13:15)
       HH=part%data(i,16)

       rpp=topart(rp-shift,xpp,ypp,zpp)
       rho=sqrt(rpp(1)**2+rpp(2)**2)

       if (abs(rpp(3)).lt.HH/2) then
          test=[rpp(1)*rr/rho,rpp(2)*rr/rho,rpp(3)]
          p(1,:)=toglob(test,xpp,ypp,zpp)+shift
          call accum(p,1,rp,1,0)
       endif

       offset=rpp-[0.0,0.0,HH/2]
       dist=mag(offset)
       offset=offset*rr/dist
       test=offset+[0.0,0.0,HH/2]

       p(1,:)=toglob(test,xpp,ypp,zpp)+shift

       offset=rpp-[0.0,0.0,-HH/2]
       dist=mag(offset)
       offset=offset*rr/dist
       test=offset+[0.0,0.0,-HH/2]

       p(2,:)=toglob(test,xpp,ypp,zpp)+shift

       call accum(p,2,rp,1,0)
    enddo

    call accres(snapneon,what,region)
    snapneon=toglob(snapneon,part%xp,part%yp,part%zp)+part%shift

  end function snapneon    

! ######################################################################

  function snapplate(r,part,flag,what,region)

    real r(3),snapplate(3)
    type (Tpart) part
    integer flag,what,region

    real Lx,Ly,rp(3)
    real p(4,3)

! call with flag=2 for rectangle

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    Lx=part%fab(1)
    Ly=part%fab(2)

    call accst

    if (flag.eq.1) then

       if (abs(rp(1)).lt.Lx/2.and.abs(rp(2)).lt.Ly/2) then
          p(1,:)=[rp(1),rp(2),0.0]
          call accum(p,1,rp,1,0)
       endif

    endif

    if (flag.eq.1.or.flag.eq.2) then

       if (abs(rp(1)).lt.Lx/2) then
          p(1,:)=[rp(1), Ly/2,0.0]
          p(2,:)=[rp(1),-Ly/2,0.0]
          call accum(p,2,rp,2,0)
       endif
       
       if (abs(rp(2)).lt.Ly/2) then
          p(1,:)=[ Lx/2,rp(2),0.0]
          p(2,:)=[-Lx/2,rp(2),0.0]
          call accum(p,2,rp,2,0)
       endif
    endif

    p(1,:)=[ Lx/2, Ly/2,0.0]
    p(2,:)=[-Lx/2, Ly/2,0.0]
    p(3,:)=[ Lx/2,-Ly/2,0.0]
    p(4,:)=[-Lx/2,-Ly/2,0.0]

    call accum(p,4,rp,3,0)

    call accres(snapplate,what,region)

    snapplate=toglob(snapplate,part%xp,part%yp,part%zp)+part%shift

  end function snapplate

! ######################################################################

  function snapdisk(r,part,flag,what,region)

    real r(3),snapdisk(3)
    type (Tpart) part
    integer flag,what,region

    real RR,dist,rp(3)
    real p(1,3)

! call with flag=2 for circle

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    RR=part%fab(1)

    dist=sqrt(rp(1)**2+rp(2)**2)

    call accst

    if (flag.eq.1) then

       if (dist.lt.RR) then
          p(1,:)=[rp(1),rp(2),0.0]
          call accum(p,1,rp,1,0)
       endif
    endif

    p(1,:)=[rp(1)*RR/dist,rp(2)*RR/dist,0.0]

    call accum(p,1,rp,2,0)

    call accres(snapdisk,what,region)

    snapdisk=toglob(snapdisk,part%xp,part%yp,part%zp)+part%shift

  end function snapdisk

! ######################################################################

  function snapbubble(r,part,flag,what,region)
    
    real r(3),snapbubble(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),RR,dist

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    dist=sqrt(rp(1)**2+rp(2)**2+rp(3)**2)
    RR=part%fab(1)

    snapbubble=rp*RR/dist
    snapbubble=toglob(snapbubble,part%xp,part%yp,part%zp)+part%shift
    
    what=1
    region=0

  end function snapbubble

! ######################################################################

  function snapline(r,part,flag,what,region)

    real r(3),snapline(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),LL
    real p(2,3)

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)
    LL=part%fab(1)
    call accst

    if (flag.eq.2) then
       if (abs(rp(3)).lt.LL/2) then
          p(1,:)=[0.0,0.0,rp(3)]
          call accum(p,1,rp,2,0)
       endif
    endif

    p(1,:)=[0.0,0.0, LL/2]
    p(2,:)=[0.0,0.0,-LL/2]

    call accum(p,2,rp,3,0)
    call accres(snapline,what,region)

    snapline=toglob(snapline,part%xp,part%yp,part%zp)+part%shift
  end function snapline

! ######################################################################

  function snaparc(r,part,flag,what,region)

    real r(3),snaparc(3)
    type (Tpart) part
    integer flag,what,region

    real rp(3),RR,theta1,theta2,theta,rho
    real p(2,3)

    rp=topart(r-part%shift, part%xp, part%yp, part%zp)

    RR=part%fab(1)
    theta1=part%fab(2)*pi/180
    theta2=part%fab(3)*pi/180

    theta=atan2(rp(2),rp(1))
    rho=sqrt(rp(1)**2+rp(2)**2)

    call accst

    if (flag.eq.2) then
       if (between(theta,theta1,theta2)) then
          p(1,:)=[rp(1),rp(2),0.0]*RR/rho
          call accum(p,1,rp,2,0)
       endif
    endif

    p(1,:)=[RR*cos(theta1),RR*sin(theta1),0.0]
    p(2,:)=[RR*cos(theta2),RR*sin(theta2),0.0]

    call accum(p,2,rp,3,0)

    call accres(snaparc,what,region)

    snaparc=toglob(snaparc,part%xp,part%yp,part%zp)+part%shift
  end function snaparc

! ######################################################################

  subroutine toupper(string)

    character(*) string
    character(26), parameter :: lower='abcdefghijklmnopqrstuvwxyz'
    character(26), parameter :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    character(1) c
    integer i,j

    do i=1,len_trim(string)
       c=string(i:i)
       j=index(lower,c)
       if (j.ne.0) string(i:i)=upper(j:j)
    enddo
  end subroutine toupper

  subroutine tolower(string)

    character(*) string
    character(26), parameter :: lower='abcdefghijklmnopqrstuvwxyz'
    character(26), parameter :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    character(1) c
    integer i,j

    do i=1,len_trim(string)
       c=string(i:i)
       j=index(upper,c)
       if (j.ne.0) string(i:i)=lower(j:j)
    enddo
  end subroutine tolower

end module snap
