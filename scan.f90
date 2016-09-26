! Part of Zinc FE package. Author: John Blackburn

module scan

! ----------------------------------------------------------------------
! Set u and ftol and variables in "common" BEFORE calling linescan, planescan
! Note that one of the args of line/planescan is type(Texpr)
! so this type is made available as well
! ----------------------------------------------------------------------

  use evaluate, only : Tlist

  implicit none

  save
  private
  public u,ftol,linescan,planescan,surfint,volint,Texpr,findu,fdchk,replaceBrace

  double precision, allocatable :: u(:,:,:,:)
  double precision ftol

  type Texpr
     double precision expr
     character(1000) sexpr
     type(Tlist) texpr
     integer iexpr
     
     integer nbrak
     integer, allocatable :: nperbrak(:)         ! no. terms in each bracket
     
     character(1000), allocatable :: sbrak(:,:)  ! string for each bracket term
     double precision, allocatable :: brak(:,:)  ! value for each bracket term
     integer, allocatable :: ibrak(:,:)          ! index for each bracket term
     integer, allocatable :: rbrak(:,:)          ! region for each bracket term
     type(Tlist), allocatable :: tbrak(:,:)      ! token list for each bracket term
  end type Texpr

contains
  
! ######################################################################

subroutine getu(i,j,k,x,y,z,ur,dur,xi,eta,mu,ifail)

! ----------------------------------------------------------------------
!     For a given point r=(x,y,z) in a given element specified by its
!     low indices i,j,k, calculate the value of u_ii and its derivatives
!     dur(1,2)=du1/dx2 etc. 
!     If point not in box containing element, ifail=2
!     If point is in box but not element return ifail=1
!     success: ifail=0
!     shares block fcnb with sub fcn
!     Uses u,ftol from common block
! ----------------------------------------------------------------------

  use fcnb
  use shape
  use geom
  use common, only : nvar,imax,jmax,kmax,rnode

  integer i,j,k,ifail
  double precision x,y,z
  double precision xi,eta,mu,ur(nvar),dur(nvar,3)      ! dimensioned from common

  integer, parameter :: lwa=3*(3+5)+3 ! n*(m+5)+m for iopt=1 or 2 (dnls1e)
  integer iw(3),itab(8,4,3),itabsrt(8,3),icomb(2,3)
  
  double precision rtab(8,3),wa(lwa),xvec(3),fvec(3)
  double precision ul(8,nvar),dN(3)          ! ul auto array
  double precision w(100)
  double precision jac(3,3),invjac(3,3),jacdet
  
  integer ifirst,l,inode,jnode,knode,m,n,nprint,lw,ii,jj,info
  double precision x1,y1,z1,xmin,ymin,zmin,xmax,ymax,zmax,tol

!  external fcn,fcn2,fcnnag
 
! ----------------------------------------------------------------------
!     Get element and its node u values
!     put copy of point on fcnb common
! ----------------------------------------------------------------------

  icomb(1,1)=i
  icomb(1,2)=j
  icomb(1,3)=k
  
  icomb(2,1)=i+1
  icomb(2,2)=j+1
  icomb(2,3)=k+1
  
  call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)

! ----------------------------------------------------------------------
!     First check if point x,y,z is inside box containing element
!     If not, return will ifail=2
! ----------------------------------------------------------------------

  ifirst=1
  do l=1,8
     x1=rtab(l,1)
     y1=rtab(l,2)
     z1=rtab(l,3)
     
     if (ifirst.eq.1) then
        xmax=x1
        ymax=y1
        zmax=z1
        
        xmin=x1
        ymin=y1
        zmin=z1
        
        ifirst=0
     else
        xmax=max(x1,xmax)
        ymax=max(y1,ymax)
        zmax=max(z1,zmax)
        
        xmin=min(x1,xmin)
        ymin=min(y1,ymin)
        zmin=min(z1,zmin)
     endif
  enddo
  
  if (.not.(x.le.xmax.and.x.ge.xmin.and. &
            y.le.ymax.and.y.ge.ymin.and. &
            z.le.zmax.and.z.ge.zmin)) then
     ifail=2
     return
  endif

! ----------------------------------------------------------------------
!     Prepare element in sorted order, get u at each local node (if allocated)
! ----------------------------------------------------------------------

  call sort(itab,rtab,itabsrt,rtabsrt) ! rtabsrt passed to fcn

  if (allocated(u)) then
  
     do l=1,8
        inode=itabsrt(l,1)
        jnode=itabsrt(l,2)
        knode=itabsrt(l,3)
        
        do ii=1,nvar
           ul(l,ii)=u(inode,jnode,knode,ii)
        enddo
     enddo

  endif

  xx=x                      ! passed to fcn
  yy=y
  zz=z

! ----------------------------------------------------------------------
!     Initialise xvec, call dnls1e and recover solution xi,eta,mu
!     fcn just provides fve! (iopt=1), fcn2 also jacobian (iopt=2)
!     fcnnag for NAG library call
!     Uncomment one of these
!     Recommend use dnls1e with fcn (iopt=1) as others did not work
!     correctly
! ----------------------------------------------------------------------

  xvec(1)=0
  xvec(2)=0
  xvec(3)=0

  m=3
  n=3
  tol=1e-6
  nprint=-1
  lw=100
  
  call dnls1e(fcn, 1,m,n,xvec,fvec,tol,nprint,info,iw,wa,lwa)
!     call dnls1e(fcn2,2,m,n,xvec,fvec,tol,nprint,info,iw,wa,lwa)
!     call E04FYF(m,n,fcnnag,xvec,fsumsq,w,lw,iuser,user,ifail)

  xi=xvec(1)
  eta=xvec(2)
  mu=xvec(3)
  
! ----------------------------------------------------------------------
!     If local coords in range then call jacobian (for benefit of getdN)
!     Hence prepare position and derivatives. ftol gives margin for error
! ----------------------------------------------------------------------

  if ( xi .ge.-1-ftol.and.xi .le.+1+ftol.and.&
       eta.ge.-1-ftol.and.eta.le.+1+ftol.and.&
       mu .ge.-1-ftol.and.mu .le.+1+ftol) then
     ifail=0
  else
     ifail=1
     return
  endif

  call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet)
  
  if (allocated(u)) then

     do ii=1,nvar
        ur(ii)=0
        do jj=1,3
           dur(ii,jj)=0
        enddo
     enddo
     
     do l=1,8
        call getdN(l,xi,eta,mu,invjac,dN)
        
        do ii=1,nvar
           ur(ii)=ur(ii)+ul(l,ii)*Nl(l,xi,eta,mu)
        enddo
        
        do ii=1,nvar
           do jj=1,3
              dur(ii,jj)=dur(ii,jj)+ul(l,ii)*dN(jj)
           enddo
        enddo
     enddo

  endif

end subroutine getu

! ######################################################################

subroutine fcn(iflag,m,n,xvec,fvec,fjac,ldfjac)

! ----------------------------------------------------------------------
!     set fvec as diff between x1 for xvec=(xi,eta,mu) and xx target etc
!     suitable for dnls1e with IOPT=1, ie, onlt function passed not jacobian
! ----------------------------------------------------------------------

  use fcnb
  use shape

  double precision xvec(n),fvec(m),fjac
  integer iflag,m,n,ldfjac,l

  double precision xi,eta,mu,x1,y1,z1
  
  xi=xvec(1)
  eta=xvec(2)
  mu=xvec(3)
  
  x1=0
  y1=0
  z1=0

  do l=1,8
     x1=x1+rtabsrt(l,1)*Nl(l,xi,eta,mu)
     y1=y1+rtabsrt(l,2)*Nl(l,xi,eta,mu)
     z1=z1+rtabsrt(l,3)*Nl(l,xi,eta,mu)
  enddo
  
  fvec(1)=x1-xx
  fvec(2)=y1-yy
  fvec(3)=z1-zz
  
end subroutine fcn

! ######################################################################

subroutine fcnnag(m,n,xvec,fvec,iuser,user)

! ----------------------------------------------------------------------
!     set fvec as difference between x1 for xvec=(xi,eta,mu) and xx target
!     Suitable for nag function E04FYF
! ----------------------------------------------------------------------

  use fcnb
  use shape
  
  integer m,n,iuser,l
  double precision xvec(n),fvec(m),user

  double precision xi,eta,mu,x1,y1,z1

    xi=xvec(1)
  eta=xvec(2)
  mu=xvec(3)
  
  x1=0
  y1=0
  z1=0
  
  do l=1,8
     x1=x1+rtabsrt(l,1)*Nl(l,xi,eta,mu)
     y1=y1+rtabsrt(l,2)*Nl(l,xi,eta,mu)
     z1=z1+rtabsrt(l,3)*Nl(l,xi,eta,mu)
  enddo
  
  fvec(1)=x1-xx
  fvec(2)=y1-yy
  fvec(3)=z1-zz
  
end subroutine fcnnag

! ######################################################################

subroutine fcn2(iflag,m,n,xvec,fvec,fjac,ldfjac)

! ----------------------------------------------------------------------
!     set fvec as difference between x1 for xvec=(xi,eta,mu) and xx target
!     This version also prepares jacobian at specified point
!     Suitable for dnls1e IOPT=2
!     calls functions Nl, dNdxi, dNdeta, dNdmu
! ----------------------------------------------------------------------

  use fcnb
  use shape

  double precision xvec(n),fvec(m),fjac(ldfjac,n)
  integer iflag,m,n,ldfjac,l,i,j

  double precision xi,eta,mu,x1,y1,z1

  xi=xvec(1)
  eta=xvec(2)
  mu=xvec(3)
  
  if (iflag.eq.1) then
     x1=0
     y1=0
     z1=0
     
     do l=1,8
        x1=x1+rtabsrt(l,1)*Nl(l,xi,eta,mu)
        y1=y1+rtabsrt(l,2)*Nl(l,xi,eta,mu)
        z1=z1+rtabsrt(l,3)*Nl(l,xi,eta,mu)
     enddo
     
     fvec(1)=x1-xx
     fvec(2)=y1-yy
     fvec(3)=z1-zz
     
  else if (iflag.eq.2) then
     do i=1,3
        do j=1,3
           fjac(i,j)=0
        enddo
     enddo
     
     do i=1,3
        do l=1,8
           fjac(i,1)=fjac(i,1)+rtabsrt(l,i)*dNdxi(l,xi,eta,mu)
           fjac(i,2)=fjac(i,2)+rtabsrt(l,i)*dNdeta(l,xi,eta,mu)
           fjac(i,3)=fjac(i,3)+rtabsrt(l,i)*dNdmu(l,xi,eta,mu)
        enddo
     enddo
  endif

end subroutine fcn2
      
! ######################################################################

subroutine linescan(iunit,N,x1,y1,z1,x2,y2,z2,expr)

! ----------------------------------------------------------------------
!     Scan between points and write out to unit iunit. N intervals
!     a little bit smart in that it tries the current element first
!     Uses ijkmax,nvar from common
! ----------------------------------------------------------------------

  use evaluate, only : Tlist
  use common, only : nvar,imax,jmax,kmax,iregup

  integer iunit,N
  double precision x1,y1,z1,x2,y2,z2
  type(Texpr) expr

  double precision ur(nvar),dur(nvar,3) ! auto arrays

  double precision xi,eta,mu,dx,dy,dz,x,y,z,pval,dist,nx,ny,nz
  integer iold,jold,kold,is,irege,ios,ifail

  nx=0; ny=0; nz=0
  
  dx=(x2-x1)/N
  dy=(y2-y1)/N
  dz=(z2-z1)/N
  
  iold=imax/2
  jold=jmax/2
  kold=kmax/2
  
! ----------------------------------------------------------------------
!     Linescan. iold, jold, kold is last result. We will try this first
! ----------------------------------------------------------------------

  do is=0,N
     
     x=x1+is*dx
     y=y1+is*dy
     z=z1+is*dz
     
     call findu(iold,jold,kold,x,y,z,ur,dur,xi,eta,mu,ifail)
     
! ----------------------------------------------------------------------
!     If success, ifail=0 and we jumped to 1. If failure, the loop
!     termination and ifail=1,2. In that case, set ur,dur,irege=0
! ----------------------------------------------------------------------

     if (ifail.eq.0) then  
        irege=iregup(iold,jold,kold)
        pval=getval(x,y,z,nx,ny,nz,ur,dur,expr,irege)
     else
        irege=0
        pval=0
     endif

     dist=sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)

     write (iunit,'(5e13.5,5i6)',iostat=ios) dist,x,y,z,pval,ifail,iold,jold,kold,irege
     
     if (ios.ne.0) then
        print *,'Error writing to linescan.out'
        call exit(1)
     endif

  enddo

end subroutine linescan

! ######################################################################

subroutine planescan(iunit,inorm,pnorm,a1,a2,b1,b2,Na,Nb,expr)

! ----------------------------------------------------------------------
!     Scan plane normal to inorm, 1=x, 2=y, 3=z, at pnorm
!     number of intervals is Na, Nb between [a1,a2] etc
!     eg, inorm=3, plane at z=pnorm, scan x=[a1,a2], y=[b1,b2]
!     uses only imax etc from common
! ----------------------------------------------------------------------

  use evaluate, only : Tlist
  use common, only : nvar,imax,jmax,kmax,iregup

  integer iunit,inorm,Na,Nb
  double precision pnorm,a1,a2,b1,b2,nx,ny,nz
  type(Texpr) expr

  double precision ur(nvar),dur(nvar,3),r(3) ! auto arrays

  double precision xi,eta,mu,x,y,z,pval
  double precision da,db,a,b
  integer iold,jold,kold,irege,ios
  integer iac,ibc,ia,ib,ifail

  nx=0; ny=0; nz=0

  if (inorm.eq.1) then
     iac=2
     ibc=3
  else if (inorm.eq.2) then
     iac=1
     ibc=3
  else
     iac=1
     ibc=2
  endif
  
  da=(a2-a1)/Na
  db=(b2-b1)/Nb
  
  iold=imax/2
  jold=jmax/2
  kold=kmax/2
  
  do ia=0,Na
     do ib=0,Nb
        a=a1+ia*da
        b=b1+ib*db
        
        r(iac)=a
        r(ibc)=b
        r(inorm)=pnorm
        
        x=r(1)
        y=r(2)
        z=r(3)
        
        call findu(iold,jold,kold,x,y,z,ur,dur,xi,eta,mu,ifail)
        
        if (ifail.eq.0) then  
           irege=iregup(iold,jold,kold)
           pval=getval(x,y,z,nx,ny,nz,ur,dur,expr,irege)
        else
           irege=0
           pval=0
        endif

        write (iunit,'(4e13.5,5i6)',iostat=ios) x,y,z,pval,ifail,iold,jold,kold,irege

        if (ios.ne.0) then
           print *,'Error writing to planescan.out'
           call exit(1)
        endif

     enddo
     write (iunit,*)
  enddo

end subroutine planescan

! ######################################################################

  function surfint(ifrom,ito,expr)

! ----------------------------------------------------------------------
! Eg: surfint 2-1 = eps1*(Vx*nx+Vy*ny+Vz*nz)
! itabsrt, rtabsrt is face
! itabsrt1,2,el rtabsrt1,2,el are for element (same for irege12)
! centre12 are normally element centres, facecentre is obvious
! ul is for element, ulface is for face
! ----------------------------------------------------------------------

    use common, only : imax,jmax,kmax,rnode,iregup,nvar,gauss,wt,ng_xi,ng_mu,ng_eta
    use geom
    use integrate
    use matrices, only : lookind
    use util
    use shape
    use iofile, only : prterr

    double precision surfint
    type(Texpr) expr
    integer ifrom,ito

    integer idir,itop,jtop,ktop,i,j,k,irege1,irege2
    integer itabsrt(4,3),itabsrt1(8,3),itabsrt2(8,3),icomb(2,3),itabsrtel(8,3)
    integer inode,jnode,knode,itab(8,4,3)
    integer ixi,ieta,icoord,ii,jj,l,facelst(4)

    double precision ul(8,nvar),ulface(4,nvar),ur(nvar),dur(nvar,3)      ! auto arrays
    double precision rtabsrt(4,3),rtab(8,3),rtabsrt1(8,3),rtabsrt2(8,3),rtabsrtel(8,3)
    double precision facecentre(3),centre1(3),centre2(3)
    double precision vec1(3),vec2(3),norm(3),norm1(3),norm2(3),norm3(3),norm4(3),test(3),eta,xi
    double precision drdxi(3),drdeta(3),trailer,coord,nx,ny,nz
    double precision jac(3,3),invjac(3,3),jacdet,dN(3),x,y,z

!    print *,'surfint',ifrom,ito

    surfint=0
    
    do idir=1,3
       
       if (idir.eq.1) then
          itop=imax-1
          jtop=jmax-1
          ktop=kmax
       else if (idir.eq.2) then
          itop=imax-1
          ktop=kmax-1
          jtop=jmax
       else
          jtop=jmax-1
          ktop=kmax-1
          itop=imax
       endif

! ----------------------------------------------------------------------
! Loop over all faces. get face geometry in itabsrt, rtabsrt
! get the region number and element setup in elements either side of the face
! elements stored in itabsrt12, rtabsrt12. Centres in centre12
! face centre is facecentre
! ----------------------------------------------------------------------

       do i=0,itop
       do j=0,jtop
       do k=0,ktop
          
          if (idir.eq.1) then
             rtabsrt(1,:)=rnode(i,j,k,:)
             rtabsrt(2,:)=rnode(i+1,j,k,:)
             rtabsrt(3,:)=rnode(i+1,j+1,k,:)
             rtabsrt(4,:)=rnode(i,j+1,k,:)
             
             itabsrt(1,:)=[i,j,k]
             itabsrt(2,:)=[i+1,j,k]
             itabsrt(3,:)=[i+1,j+1,k]
             itabsrt(4,:)=[i,j+1,k]
             
             facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4

             if (k.eq.kmax) then
                irege1=lookind('ZMAX')
                centre1=facecentre
             else
                irege1=iregup(i,j,k)
                centre1=centre(i,j,k)
                icomb(1,:)=[i,j,k]
                icomb(2,:)=[i+1,j+1,k+1]
                call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
                call sort(itab,rtab,itabsrt1,rtabsrt1)
             endif
             
             if (k.eq.0) then
                irege2=lookind('ZMIN')
                centre2=facecentre
             else
                irege2=iregup(i,j,k-1)
                centre2=centre(i,j,k-1)
                icomb(1,:)=[i,j,k-1]
                icomb(2,:)=[i+1,j+1,k]
                call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
                call sort(itab,rtab,itabsrt2,rtabsrt2)
             endif

! ----------------------------------------------------------------------
             
          else if (idir.eq.2) then
             rtabsrt(1,:)=rnode(i,j,k,:)
             rtabsrt(2,:)=rnode(i+1,j,k,:)
             rtabsrt(3,:)=rnode(i+1,j,k+1,:)
             rtabsrt(4,:)=rnode(i,j,k+1,:)
             
             itabsrt(1,:)=[i,j,k]
             itabsrt(2,:)=[i+1,j,k]
             itabsrt(3,:)=[i+1,j,k+1]
             itabsrt(4,:)=[i,j,k+1]
             
             facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4

             if (j.eq.jmax) then
                irege1=lookind('YMAX')
                centre1=facecentre
             else
                irege1=iregup(i,j,k)
                centre1=centre(i,j,k)
                icomb(1,:)=[i,j,k]
                icomb(2,:)=[i+1,j+1,k+1]
                call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
                call sort(itab,rtab,itabsrt1,rtabsrt1)
             endif
             
             if (j.eq.0) then
                irege2=lookind('YMIN')
                centre2=facecentre
             else
                irege2=iregup(i,j-1,k)
                centre2=centre(i,j-1,k)
                icomb(1,:)=[i,j-1,k]
                icomb(2,:)=[i+1,j,k+1]
                call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
                call sort(itab,rtab,itabsrt2,rtabsrt2)
             endif

! ----------------------------------------------------------------------

          else
             rtabsrt(1,:)=rnode(i,j,k,:)
             rtabsrt(2,:)=rnode(i,j+1,k,:)
             rtabsrt(3,:)=rnode(i,j+1,k+1,:)
             rtabsrt(4,:)=rnode(i,j,k+1,:)
             
             itabsrt(1,:)=[i,j,k]
             itabsrt(2,:)=[i,j+1,k]
             itabsrt(3,:)=[i,j+1,k+1]
             itabsrt(4,:)=[i,j,k+1]

             facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4
             
             if (i.eq.imax) then
                irege1=lookind('XMAX')
                centre1=facecentre
             else
                irege1=iregup(i,j,k)
                centre1=centre(i,j,k)
                icomb(1,:)=[i,j,k]
                icomb(2,:)=[i+1,j+1,k+1]
                call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
                call sort(itab,rtab,itabsrt1,rtabsrt1)
             endif
             
             if (i.eq.0) then
                irege2=lookind('XMIN')
                centre2=facecentre
             else
                irege2=iregup(i-1,j,k)
                centre2=centre(i-1,j,k)
                icomb(1,:)=[i-1,j,k]
                icomb(2,:)=[i,j+1,k+1]
                call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
                call sort(itab,rtab,itabsrt2,rtabsrt2)
             endif
          endif

! ----------------------------------------------------------------------
! If this is a face to integrate, find the unit normal vector
! discover icoord and coord indicating which face of the "to" element
! is actuall the current face. icoord = 1,2,3 => xi,eta,mu
! also get ul for relevant element
! ----------------------------------------------------------------------

          if ( (irege1.eq.ifrom.and.irege2.eq.ito).or. &
               (irege2.eq.ifrom.and.irege1.eq.ito)) then

!             print *,'found integration face:',irege1,irege2,idir,i,j,k

             vec1=rtabsrt(2,:)-rtabsrt(1,:)
             vec2=rtabsrt(4,:)-rtabsrt(1,:)
             norm1=cross(vec1,vec2)

             vec1=rtabsrt(3,:)-rtabsrt(2,:)
             vec2=rtabsrt(1,:)-rtabsrt(2,:)
             norm2=cross(vec1,vec2)

             vec1=rtabsrt(4,:)-rtabsrt(3,:)
             vec2=rtabsrt(2,:)-rtabsrt(3,:)
             norm3=cross(vec1,vec2)

             vec1=rtabsrt(1,:)-rtabsrt(4,:)
             vec2=rtabsrt(3,:)-rtabsrt(4,:)
             norm4=cross(vec1,vec2)
             
             if (dot_product(norm1,norm2).lt.0.or.dot_product(norm1,norm3).lt.0.or.dot_product(norm1,norm4).lt.0.or.&
                 dot_product(norm2,norm3).lt.0.or.dot_product(norm2,norm4).lt.0.or.dot_product(norm3,norm4).lt.0) then
                call prterr('incompatible normal vectors')
             endif
                
             norm=(norm1+norm2+norm3+norm4)/4
             norm=norm/mag(norm)
     
             if (irege1.eq.ifrom.and.irege2.eq.ito) then
                test=centre2-centre1
                call lockface(itabsrt2,itabsrt,icoord,coord,facelst)

                do l=1,8
                   inode=itabsrt2(l,1)
                   jnode=itabsrt2(l,2)
                   knode=itabsrt2(l,3)
                   
                   do ii=1,nvar
                      ul(l,ii)=u(inode,jnode,knode,ii)
                   enddo
                enddo

                rtabsrtel=rtabsrt2
                itabsrtel=itabsrt2

             else
                test=centre1-centre2
                call lockface(itabsrt1,itabsrt,icoord,coord,facelst)

                do l=1,8
                   inode=itabsrt1(l,1)
                   jnode=itabsrt1(l,2)
                   knode=itabsrt1(l,3)
                   
                   do ii=1,nvar
                      ul(l,ii)=u(inode,jnode,knode,ii)
                   enddo
                enddo

                rtabsrtel=rtabsrt1
                itabsrtel=itabsrt1
             endif

             if (dot_product(norm,test).lt.0) norm=-norm

! ----------------------------------------------------------------------
! Get the solution at the face nodes
! ----------------------------------------------------------------------
        
             do l=1,4
                inode=itabsrtel(facelst(l),1)
                jnode=itabsrtel(facelst(l),2)
                knode=itabsrtel(facelst(l),3)
                
                do ii=1,nvar
                   ulface(l,ii)=u(inode,jnode,knode,ii)
                enddo
             enddo
             
!             print *,'ulface',ulface,ng

! ----------------------------------------------------------------------
! Integrate
! ----------------------------------------------------------------------
        
             do ixi=1,ng_xi
             do ieta=1,ng_eta
!                print *,'ixi,ieta',ixi,ieta
                xi=gauss(ng_xi,ixi)
                eta=gauss(ng_eta,ieta)
                
                drdxi=0
                drdeta=0
                
                do l=1,4
                   drdxi(:)= drdxi(:)+ rtabsrtel(facelst(l),:)*dN2d_dxi(l,xi,eta)
                   drdeta(:)=drdeta(:)+rtabsrtel(facelst(l),:)*dN2d_deta(l,xi,eta)
                enddo

                trailer=mag(cross(drdxi,drdeta))*wt(ng_xi,ixi)*wt(ng_eta,ieta)
                
                ur=0
                do ii=1,nvar
                   do l=1,4
                      ur(ii)=ur(ii)+ulface(l,ii)*N2d(l,xi,eta)
                   enddo
                enddo

                if (icoord.eq.1) then
                   call jacobian(coord,xi,eta,rtabsrtel,jac,invjac,jacdet)
                else if (icoord.eq.2) then
                   call jacobian(xi,coord,eta,rtabsrtel,jac,invjac,jacdet)
                else
                   call jacobian(xi,eta,coord,rtabsrtel,jac,invjac,jacdet)
                endif

                dur(:,:)=0
            
                do l=1,8
                   if (icoord.eq.1) then
                      call getdN(l,coord,xi,eta,invjac,dN)
                   else if (icoord.eq.2) then
                      call getdN(l,xi,coord,eta,invjac,dN)
                   else
                      call getdN(l,xi,eta,coord,invjac,dN)
                   endif
                   
                   do ii=1,nvar
                      do jj=1,3
                         dur(ii,jj)=dur(ii,jj)+ul(l,ii)*dN(jj)
                      enddo
                   enddo
                enddo

                x=facecentre(1)
                y=facecentre(2)
                z=facecentre(3)

                nx=norm(1)
                ny=norm(2)
                nz=norm(3)

!                print *,'integrate: ',xi,eta,trailer
                surfint=surfint+getval(x,y,z,nx,ny,nz,ur,dur,expr,ito)*trailer

             enddo
             enddo

!             stop
          endif

       enddo
       enddo
       enddo
    enddo

!    print *,'return=',surfint

  end function surfint

! ######################################################################

  function volint(region,expr)

    use common
    use geom
    use shape
    use indexq, only : indQ

    integer region
    type(Texpr) expr
    double precision volint

    integer i,j,k,irege,icomb(2,3),itabsrt(8,3),itab(8,4,3)
    double precision rtab(8,3),rtabsrt(8,3)
    double precision ul(8,nvar),nx,ny,nz,x,y,z,ur(nvar),dur(nvar,3)

    integer ixi,ieta,imu,ii,jj,inode,jnode,knode,l
    double precision xi,eta,mu,trailer,dN(3),Nl1
    double precision jac(3,3),invjac(3,3),jacdet

    nx=0
    ny=0
    nz=0

    volint=0

    do i=0,imax-1
    do j=0,jmax-1
    do k=0,kmax-1

! ----------------------------------------------------------------------
! Get region number
! ----------------------------------------------------------------------

       irege=iregup(i,j,k)

       if (irege.eq.region) then
     
! ----------------------------------------------------------------------
! Get the sorted element nodes. Get local u values in ul(l,ii)
! ----------------------------------------------------------------------

          icomb(1,1)=i
          icomb(1,2)=j
          icomb(1,3)=k
     
          icomb(2,1)=i+1
          icomb(2,2)=j+1
          icomb(2,3)=k+1
     
          call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
          call sort(itab,rtab,itabsrt,rtabsrt)
          
          x=sum(rtabsrt(:,1))/8
          y=sum(rtabsrt(:,2))/8
          z=sum(rtabsrt(:,3))/8
          
          do l=1,8
             inode=itabsrt(l,1)
             jnode=itabsrt(l,2)
             knode=itabsrt(l,3)
             
             do ii=1,nvar
                ul(l,ii)=u(inode,jnode,knode,ii)
             enddo
          enddo
          
! ----------------------------------------------------------------------
! Element integral loops
! Integrate over element, get du(ii,jj), ur(ii)
! ----------------------------------------------------------------------

          do ixi=1,ng_xi
             do ieta=1,ng_eta
                do imu=1,ng_mu

                   xi=gauss(ng_xi,ixi)
                   eta=gauss(ng_eta,ieta)
                   mu=gauss(ng_mu,imu)

                   call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet)
                   trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)*wt(ng_mu,imu)
            
                   dur(:,:)=0
            
                   do l=1,8
                      call getdN(l,xi,eta,mu,invjac,dN)
                      do ii=1,nvar
                         do jj=1,3
                            dur(ii,jj)=dur(ii,jj)+ul(l,ii)*dN(jj)
                         enddo
                      enddo
                   enddo

                   ur(:)=0
                   
                   do l=1,8
                      Nl1=Nl(l,xi,eta,mu)
                      do ii=1,nvar
                         ur(ii)=ur(ii)+ul(l,ii)*Nl1
                      enddo
                   enddo
                   
                   volint=volint+getval(x,y,z,nx,ny,nz,ur,dur,expr,region)*trailer

                enddo
             enddo
          enddo

       endif

    enddo
    enddo
    enddo

  end function volint

! ######################################################################

subroutine findu(iold,jold,kold,x,y,z,ur,dur,xi,eta,mu,ifail)

! ----------------------------------------------------------------------
! find ur, dur for point (x,y,z). iold,jold,kold is initial guess for element
! (bottom left corner node)
! and will be reset to the true value of the element
! ifail=0 success, ifail=1 failure
! ----------------------------------------------------------------------

  use common, only : nvar,imax,jmax,kmax

  integer iold,jold,kold,ifail
  double precision x,y,z,ur(nvar),dur(nvar,3),xi,eta,mu

  integer ijkmax,irad,i1,i2,j1,j2,k1,k2,i,j,k

  logical iedge,jedge,kedge

  ijkmax=max(imax,jmax,kmax)

  do irad=0,ijkmax
     i1=max(iold-irad,0)
     i2=min(iold+irad,imax-1)

     j1=max(jold-irad,0)
     j2=min(jold+irad,jmax-1)

     k1=max(kold-irad,0)
     k2=min(kold+irad,kmax-1)

     do i=i1,i2
        do j=j1,j2
           do k=k1,k2

              iedge=i.eq.i1.or.i.eq.i2
              jedge=j.eq.j1.or.j.eq.j2
              kedge=k.eq.k1.or.k.eq.k2

              if (iedge.or.jedge.or.kedge) then

                 call getu(i,j,k,x,y,z,ur,dur,xi,eta,mu,ifail)

                 if (ifail.eq.0) then
                    iold=i
                    jold=j
                    kold=k
                    ifail=0
                    return
                 endif
              endif

           enddo
        enddo
     enddo

  enddo

  ifail=1

end subroutine findu

! ######################################################################

function getval(x,y,z,nx,ny,nz,ur,dur,expr,irege)

  use common, only : us,dus,nvar
  use evaluate
  use iofile, only : prterr
  use matrices, only : scanfun

  double precision getval,x,y,z,nx,ny,nz,ur(nvar),dur(nvar,3)
  type(Texpr) expr
  integer irege

  double precision res,val
  integer ivar,i,j,ival
  logical found
  character(2) temp

! ----------------------------------------------------------------------
! set right values of x,y,z, nx,ny,nz, ur,dur on variable stack
! ----------------------------------------------------------------------

  call defparam('X',x)
  call defparam('Y',y)
  call defparam('Z',z)

  call defparam('NX',nx)
  call defparam('NY',ny)
  call defparam('NZ',nz)
     
  do ivar=1,nvar
     call defparam(us(ivar),ur(ivar))
     do j=1,3
        call defparam(dus(ivar,j),dur(ivar,j))
     enddo
  enddo

! ----------------------------------------------------------------------
! Go through brackets and set regexp01, regexp02 etc on stack
! with appropriate value refering to index number (3=token is not allowed)
! ----------------------------------------------------------------------

  do i=1,expr%nbrak

     found=.false.
     do j=1,expr%nperbrak(i)
        if (expr%rbrak(i,j).eq.irege) then

           ival=expr%ibrak(i,j)

           if (ival.lt.0) then
              val=ur(abs(ival))
           else if (ival.eq.0) then
              val=expr%brak(i,j)
           else if (ival.ge.1.and.ival.le.2) then
              call evaltoken(expr%tbrak(i,j),val)
           else
              call prterr('Illegal {} expression type')
           endif

           write (temp,'(i2.2)') i
           call defparam('regexp'//temp,val)
           found=.true.
           exit
        endif
     enddo
     if (.not.found) call prterr('While scanning, could not find needed region in {}')
  enddo

! ----------------------------------------------------------------------
! no brackets, evaluate expression based on index (incl. token=3)
! with brackets, use evaltoken with regexp01 etc set above (ignore %iexpr)
! Note that scanfun may be an actual function or a procedure pointer
! depending on matrices
! ----------------------------------------------------------------------

  if (expr%nbrak.eq.0) then
     ival=expr%iexpr

     if (ival.lt.0) then
        getval=ur(abs(ival))
     else if (ival.eq.0) then
        getval=expr%expr
     else if (ival.eq.1.or.ival.eq.2) then
        call evaltoken(expr%texpr,res)
        getval=res
     else if (ival.eq.3) then
        getval=scanfun(expr%sexpr,x,y,z,nx,ny,nz,ur,dur,nvar)
     else
        call prterr('Illegal ival in getval')
     endif
  else
     call evaltoken(expr%texpr,res)
     getval=res
  endif

end function getval

! ######################################################################

subroutine fdchk(iunit,inorm,ipnorm,ii)

! ----------------------------------------------------------------------
! Check equation ii, plane normal to inorm (1=i, 2=j, 3=k), (i,j,k)=ipnorm
! Eg FDCHK 3 5 1, check plane k=5, equation 1
! ----------------------------------------------------------------------

  use common, only : us,imax,jmax,kmax,nvar,rnode

  integer iunit,inorm,ipnorm,ii
  integer nmax,mmax,m,n,jj,kk,ll
  double precision Cr,Cl,dxr,dxl,dx,Ct(nvar),At(nvar),Ft,ur,ul,uru,urd,ulu,uld,uc
  integer, dimension(3) :: ix,ixr,ixl,ixru,ixrd,ixlu,ixld

  ! testing
  integer i,j,k
  double precision c11,c66,c55,c12,c13,e31,e15,dy,dz,Cux,Cuy,Cuz,CV

  write (iunit,'(3a5,100a13)') '#   i','j','k',('C.'//trim(us(kk)),kk=1,nvar), & 
       ('A.'//trim(us(kk)),kk=1,nvar),'F'

  if (inorm.eq.1) then
     mmax=jmax-1
     nmax=kmax-1
  elseif (inorm.eq.2) then
     mmax=imax-1
     nmax=kmax-1
  else
     mmax=imax-1
     nmax=jmax-1
  endif

  do m=1,mmax
  do n=1,nmax

     if (inorm.eq.1) then
        ix(1)=ipnorm
        ix(2)=m
        ix(3)=n
     elseif (inorm.eq.2) then
        ix(1)=m
        ix(2)=ipnorm
        ix(3)=n
     else
        ix(1)=m
        ix(2)=n
        ix(3)=ipnorm
     endif

     Ct=0
     At=0

     do jj=1,3
     do kk=1,nvar
     do ll=1,3

        call compass(jj,ll,ix, ixr,ixl,ixru,ixrd,ixlu,ixld)

        if (jj.eq.ll) then
           Cr=getMav("C",ix,ixr,ii,jj,kk,ll)       ! av of elements connecting ix and ixr
           Cl=getMav("C",ix,ixl,ii,jj,kk,ll)
           
           ur=u(ixr(1),ixr(2),ixr(3),kk)
           ul=u(ixl(1),ixl(2),ixl(3),kk)
           uc=u(ix(1),ix(2),ix(3),kk)
           
           dxr=rnode(ixr(1),ixr(2),ixr(3),jj)-rnode( ix(1), ix(2), ix(3),jj)
           dxl=rnode( ix(1), ix(2), ix(3),jj)-rnode(ixl(1),ixl(2),ixl(3),jj)

           dx=(dxr+dxl)/2

           if (dxr<0.or.dxl<0) print *,'Warning: dxr or dxl<0'
!           print *,'jj=ll',dxr,dxl,dx
           
           Ct(kk)=Ct(kk)+(Cr*(ur-uc)/dxr-Cl*(uc-ul)/dxl)/dx
        else
         
           Cr=getMav("C",ix,ixr,ii,jj,kk,ll)         ! av of all elements connecting ix and ixr
           Cl=getMav("C",ix,ixl,ii,jj,kk,ll)
           
           uru=u(ixru(1),ixru(2),ixru(3),kk)
           urd=u(ixrd(1),ixrd(2),ixrd(3),kk)
           
           ulu=u(ixlu(1),ixlu(2),ixlu(3),kk)
           uld=u(ixld(1),ixld(2),ixld(3),kk)
           
           dxr=rnode(ixru(1),ixru(2),ixru(3),ll)-rnode(ixrd(1),ixrd(2),ixrd(3),ll)
           dxl=rnode(ixlu(1),ixlu(2),ixlu(3),ll)-rnode(ixld(1),ixld(2),ixld(3),ll)
           
           dx=rnode(ixr(1),ixr(2),ixr(3),jj)-rnode(ixl(1),ixl(2),ixl(3),jj)
           
           if (dxr<0.or.dxl<0.or.dx<0) print *,'Warning: dx(r,l)<0'
!           print *,'jj <> ll',dxr,dxl,dx

           Ct(kk)=Ct(kk)+(Cr*(uru-urd)/dxr-Cl*(ulu-uld)/dxl)/dx
        endif

      enddo
      enddo
      enddo

      do kk=1,nvar
         At(kk)=getMav("A",ix,ix,ii,0,kk,0)*u(ix(1),ix(2),ix(3),kk)
      enddo

      Ft=getMav("F",ix,ix,ii,0,0,0)

      i=ix(1)
      j=ix(2)
      k=ix(3)

      c11=12.6e10
      c66=2.347e10
      c55=2.3e10

      c12=7.95e10
      c13=8.41e10

      e31=-6.5
      e15=17

      dx=0.00125/2; dy=0.00125/2; dz=0.00125/2

      Cux= (c11*(u(i+1,j,k,1)-u(i,j,k,1))/dx-c11*(u(i,j,k,1)-u(i-1,j,k,1))/dx)/dx &    ! c1111
          +(c66*(u(i,j+1,k,1)-u(i,j,k,1))/dy-c66*(u(i,j,k,1)-u(i,j-1,k,1))/dy)/dy &    ! c1212
          +(c55*(u(i,j,k+1,1)-u(i,j,k,1))/dz-c55*(u(i,j,k,1)-u(i,j,k-1,1))/dz)/dz      ! c1313

!      Cuy= (c12*(u(i+1,j+1,k,2)-u(i+1,j-1,k,2))/dy-c12*(u(i-1,j+1,k,2)-u(i-1,j-1,k,2))/dy)/dx & ! c1122
!          +(c66*(u(i+1,j+1,k,2)-u(i-1,j+1,k,2))/dx-c66*(u(i+1,j-1,k,2)-u(i-1,j-1,k,2))/dx)/dy   ! c1221

!      Cuz= (c13*(u(i+1,j,k+1,3)-u(i+1,j,k-1,3))/dz-c13*(u(i-1,j,k+1,3)-u(i-1,j,k-1,3))/dz)/dx & ! c1133
!          +(c55*(u(i+1,j,k+1,3)-u(i-1,j,k+1,3))/dx-c55*(u(i+1,j,k-1,3)-u(i-1,j,k-1,3))/dx)/dz   ! c1331

!      CV=  (e31*(u(i+1,j,k+1,4)-u(i+1,j,k-1,4))/dz-e31*(u(i-1,j,k+1,4)-u(i-1,j,k-1,4))/dz)/dx & ! e311
!          +(e15*(u(i+1,j,k+1,4)-u(i-1,j,k+1,4))/dx-e15*(u(i+1,j,k-1,4)-u(i-1,j,k-1,4))/dx)/dz   ! e113

      Cuy= (c12*(u(i+1,j+1,k,2)-u(i+1,j  ,k,2))/dy-c12*(u(i  ,j+1,k,2)-u(i  ,j  ,k,2))/dy)/dx & ! c1122
          +(c66*(u(i+1,j+1,k,2)-u(i  ,j+1,k,2))/dx-c66*(u(i+1,j  ,k,2)-u(i  ,j  ,k,2))/dx)/dy   ! c1221

      Cuz= (c13*(u(i+1,j,k+1,3)-u(i+1,j,k  ,3))/dz-c13*(u(i  ,j,k+1,3)-u(i  ,j,k  ,3))/dz)/dx & ! c1133
          +(c55*(u(i+1,j,k+1,3)-u(i  ,j,k+1,3))/dx-c55*(u(i+1,j,k  ,3)-u(i  ,j,k  ,3))/dx)/dz   ! c1331

      CV=  (e31*(u(i+1,j,k+1,4)-u(i+1,j,k  ,4))/dz-e31*(u(i  ,j,k+1,4)-u(i  ,j,k  ,4))/dz)/dx & ! e311
          +(e15*(u(i+1,j,k+1,4)-u(i  ,j,k+1,4))/dx-e15*(u(i+1,j,k  ,4)-u(i  ,j,k  ,4))/dx)/dz   ! e113


!      Cuy=Cuy/4
!      Cuz=Cuz/4
!      CV=CV/4

      write (iunit,'(3i5,100e13.5)') ix,(Ct(kk),kk=1,nvar),(At(kk),kk=1,nvar),Ft,Cux,Cuy,Cuz,CV

   enddo
   write (iunit,*)
   enddo

end subroutine fdchk

! ######################################################################

subroutine compass(m,n,ix,ixr,ixl,ixru,ixrd,ixlu,ixld)

! ----------------------------------------------------------------------
!     Taking m to be r-l and n to be u-d, return ix vectors
!     for centre, r, l, ru, rd, lu, ld. 
! ----------------------------------------------------------------------

  integer m,n
  integer, dimension(3) :: ix,ixr,ixl,ixru,ixrd,ixlu,ixld
  integer i

! ----------------------------------------------------------------------
!     r,u
! ----------------------------------------------------------------------

  do i=1,3
     ixr(i)=ix(i)
     ixl(i)=ix(i)
  enddo
  
  ixr(m)=ix(m)+1
  ixl(m)=ix(m)-1
  
  if (m.eq.n) return

! ----------------------------------------------------------------------
!     ru,rd,lu,ld
! ----------------------------------------------------------------------
  
  do i=1,3
     ixru(i)=ix(i)
     ixrd(i)=ix(i)
     ixlu(i)=ix(i)
     ixld(i)=ix(i)
  enddo

  ixru(m)=ix(m)+1
  ixru(n)=ix(n)+1
  
  ixrd(m)=ix(m)+1
  ixrd(n)=ix(n)-1
  
  ixlu(m)=ix(m)-1
  ixlu(n)=ix(n)+1
  
  ixld(m)=ix(m)-1
  ixld(n)=ix(n)-1
  
end subroutine compass

! ######################################################################

function getMav(M,ix1,ix2,ii,jj,kk,ll)

! ----------------------------------------------------------------------
! Given vector nodes ix1 and ix2 (each triplets i,j,k)
! Find average M=(C,A,F) of all elements that connect to both nodes
! For now we will just support linear systems so that C, A, F
! dont depend on ur, dur
! ----------------------------------------------------------------------

  use common, only : nvar,iregup,imax,jmax,kmax
  use matrices, only : CCval, AAval, FFval

  integer ix1(3),ix2(3),ii,jj,kk,ll,iCCx,iaax
  double precision getMav
  character(1) M

  double precision ur(nvar),dur(nvar,3),rc(3)
  integer icomb(2,3),iep,ie,je,ke,ic,jc,kc,i1,j1,k1,irege,nfound
  logical found

  integer ies(8),jes(8),kes(8)

  data ies/ 1, 1, 1, 1,-1,-1,-1,-1/
  data jes/ 1, 1,-1,-1, 1, 1,-1,-1/
  data kes/ 1,-1, 1,-1, 1,-1, 1,-1/

  ur=0       ! not used (yet)
  dur=0
  rc=0

  getMav=0
  nfound=0
  do iep=1,8
     ie=ies(iep)
     je=jes(iep)
     ke=kes(iep)

     icomb(1,1)=ix1(1)
     icomb(1,2)=ix1(2)
     icomb(1,3)=ix1(3)
     
     icomb(2,1)=ix1(1)+ie
     icomb(2,2)=ix1(2)+je
     icomb(2,3)=ix1(3)+ke

     if ( (icomb(2,1).ge.0.and.icomb(2,1).le.imax).and. &
          (icomb(2,2).ge.0.and.icomb(2,2).le.jmax).and. &
          (icomb(2,3).ge.0.and.icomb(2,3).le.kmax)) then
     
        found=.false.
        do ic=1,2
           do jc=1,2
              do kc=1,2
                 if (icomb(ic,1)==ix2(1).and.icomb(jc,2)==ix2(2).and.icomb(kc,3)==ix2(3)) found=.true.
              enddo
           enddo
        enddo

!     print *,ix1
!     print *,icomb(1,:)
!     print *,icomb(2,:)
!     print *,ix2,found
     
        if (found) then
           i1=ix1(1)+(ie-1)/2
           j1=ix1(2)+(je-1)/2
           k1=ix1(3)+(ke-1)/2
           
           irege=iregup(i1,j1,k1)

           if (irege<1 .or. irege>2) then
              print *,'irege out of range'
              stop
           endif

           if (M=='C') then
              getMav=getMav+CCval(irege,ii,jj,kk,ll,ur,dur,rc,iCCx)   ! Need ur, dur and rc
           elseif (M=='A') then
              getMav=getMav+AAval(irege,ii,kk,ur,dur,rc,iaax)
           elseif (M=="F") then
              getMav=getMav+FFVal(irege,ii,ur,dur,rc)
           else
              print *,'Error, getMav, bad M'
              stop
           endif
           
!           if (M=="C".and.irege==2) print *,ix1,'|',i1,j1,k1,irege,getMav
           
           nfound=nfound+1
        endif

     endif
  enddo    ! element loop

  if (nfound==0) then
     print *,'Error, element not found'
     stop
  endif

  getMav=getMav/nfound

!  print *,'getMav=',getMav
!  stop

end function getMav

! ######################################################################

function replaceBrace(string)

  character(*) string
  character(len(string)) replaceBrace
  integer i

  do i=1,len(string)
     if (string(i:i)=='{') then
        replaceBrace(i:i)='['
     else if (string(i:i)=='}') then
        replaceBrace(i:i)=']'
     else
        replaceBrace(i:i)=string(i:i)
     endif
  enddo
end function replaceBrace

end module scan
