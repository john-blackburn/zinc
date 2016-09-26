! Part of Zinc FE package. Author: John Blackburn

module geom

! subroutine getel
! subroutine sort
! function equal
! subroutine jacobian
! subroutine getdN

implicit none

private
public getel,sort,jacobian,getdN,centre

integer, private :: itabcopy(8,4,3)

contains

subroutine getel(icomb,rnode,imax,jmax,kmax,itab,rtab)

! ----------------------------------------------------------------------
!     Given combination array icomb, calculate itab,rtab to specify element
!     These can be passed to sort
!     Eg: 
!     icomb(1,1)=i
!     icomb(1,2)=j
!     icomb(1,3)=k
!     
!     icomb(2,1)=i+1
!     icomb(2,2)=j+1
!     icomb(2,3)=k+1
! ----------------------------------------------------------------------

  use common, only : NDIM
  use iofile, only : prterr

  integer imax,jmax,kmax
  integer icomb(2,NDIM),itab(2**NDIM,NDIM+1,NDIM)
  double precision rnode(0:imax,0:jmax,0:kmax,NDIM),rtab(2**NDIM,NDIM)
  
  integer i,j,k,ind,ic,jc,kc,ip,ip2,idiff,icol
  integer inode,jnode,knode,inode2,jnode2,knode2
  logical found

!  do i=1,8
!     do j=1,4
!        do k=1,3
!           itab(i,j,k)=-1
!        enddo
!     enddo
!  enddo

  itab=-1

! ----------------------------------------------------------------------
! 3D
! ----------------------------------------------------------------------
  
  if (NDIM==3) then
     ind=0
     do ic=1,2
        do jc=1,2
           do kc=1,2
              ind=ind+1
              itab(ind,1,1)=icomb(ic,1)
              itab(ind,1,2)=icomb(jc,2)
              itab(ind,1,3)=icomb(kc,3)
           enddo
        enddo
     enddo

     do ip=1,8
        inode=itab(ip,1,1)
        jnode=itab(ip,1,2)
        knode=itab(ip,1,3)
         
        do ip2=1,8
           inode2=itab(ip2,1,1)
           jnode2=itab(ip2,1,2)
           knode2=itab(ip2,1,3)
           idiff=abs(inode2-inode)+abs(jnode2-jnode)+abs(knode2-knode)
           if (idiff.eq.1) then

              found=.false.
              do icol=2,4
                 if (itab(ip,icol,1).eq.-1) then
                    itab(ip,icol,1)=inode2
                    itab(ip,icol,2)=jnode2
                    itab(ip,icol,3)=knode2
                    found=.true.
                    exit
                    !                 goto 1
                 endif
              enddo
              if (.not.found) call prterr('error getel 3d')
              !1          continue
           endif
        enddo
     enddo

     do ip=1,8
        inode=itab(ip,1,1)
        jnode=itab(ip,1,2)
        knode=itab(ip,1,3)
        
        rtab(ip,1)=rnode(inode,jnode,knode,1)
        rtab(ip,2)=rnode(inode,jnode,knode,2)
        rtab(ip,3)=rnode(inode,jnode,knode,3)
     enddo

! ----------------------------------------------------------------------
! 2D
! ----------------------------------------------------------------------

  else if (NDIM==2) then
     ind=0
     do ic=1,2
        do jc=1,2
           ind=ind+1
           itab(ind,1,1)=icomb(ic,1)
           itab(ind,1,2)=icomb(jc,2)
        enddo
     enddo

     do ip=1,4
        inode=itab(ip,1,1)
        jnode=itab(ip,1,2)
         
        do ip2=1,4
           inode2=itab(ip2,1,1)
           jnode2=itab(ip2,1,2)
           idiff=abs(inode2-inode)+abs(jnode2-jnode)

           if (idiff.eq.1) then
              found=.false.
              do icol=2,3
                 if (itab(ip,icol,1).eq.-1) then
                    itab(ip,icol,1)=inode2
                    itab(ip,icol,2)=jnode2
                    found=.true.
                    exit
                 endif
              enddo
              if (.not.found) call prterr('error getel 2d')
           endif
        enddo
     enddo

     do ip=1,4
        inode=itab(ip,1,1)
        jnode=itab(ip,1,2)
        
        rtab(ip,1)=rnode(inode,jnode,0,1)
        rtab(ip,2)=rnode(inode,jnode,0,2)
     enddo

! ----------------------------------------------------------------------
! 1D
! ----------------------------------------------------------------------

  else
     ind=0
     do ic=1,2
        ind=ind+1
        itab(ind,1,1)=icomb(ic,1)
     enddo

     do ip=1,2
        inode=itab(ip,1,1)
         
        do ip2=1,2
           inode2=itab(ip2,1,1)

           idiff=abs(inode2-inode)

           if (idiff.eq.1) then
              found=.false.
              do icol=2,2
                 if (itab(ip,icol,1).eq.-1) then
                    itab(ip,icol,1)=inode2
                    found=.true.
                    exit
                 endif
              enddo
              if (.not.found) call prterr('error getel 1d')
           endif
        enddo
     enddo

     do ip=1,2
        inode=itab(ip,1,1)
        
        rtab(ip,1)=rnode(inode,0,0,1)
     enddo

  endif

end subroutine getel

! ######################################################################

subroutine sort(itab,rtab,itabsrt,rtabsrt)

! ----------------------------------------------------------------------
!     Given a table of 8 numbered nodes and connection data (itab,rtab)
!     sort the nodes into order according to algorithm given in 
!     Uni colorado AFEM for hexahedra (chapter 18)
! ----------------------------------------------------------------------

  use common, only : NDIM
  use iofile, only : prterr

  integer itab(2**NDIM,NDIM+1,NDIM),itabsrt(2**NDIM,NDIM)
  double precision rtab(2**NDIM,NDIM),rtabsrt(2**NDIM,NDIM)

  integer list(2**NDIM,NDIM),listop(8,3),list2(2**NDIM,NDIM),listop2(8,3)
  double precision face(3),faceop(3),s(3),t(3),r12(NDIM),r23(NDIM)
  double precision rlist(2**NDIM,NDIM),rlistop(8,3)
  integer next(3),itest(3)

  integer i,j,k,i1,j1,ipos,ip,ifail,ind,idiff
  double precision dp
  logical found,done(4)

! ----------------------------------------------------------------------
! 1D: node 2 should be to the right of node 1 rtabsrt(2,1) > rtabsrt(1,1)
! Note : => 1:1
! ----------------------------------------------------------------------

  if (NDIM==1) then
     if (rtab(2,1)>rtab(1,1)) then
        itabsrt(1,:)=itab(1,1,:)
        itabsrt(2,:)=itab(2,1,:)
        
        rtabsrt(1,:)=rtab(1,:)
        rtabsrt(2,:)=rtab(2,:)
     else
        itabsrt(1,:)=itab(2,1,:)
        itabsrt(2,:)=itab(1,1,:)
        
        rtabsrt(1,:)=rtab(2,:)
        rtabsrt(2,:)=rtab(1,:)
     endif

! ----------------------------------------------------------------------
! 2D: Use only itab(j,1,:). [i|r]tabsrt nodes are anticlockwise
! ----------------------------------------------------------------------

  else if (NDIM==2) then
     done=.false.

     list(1,:)=itab(1,1,:)
     done(1)=.true.

     do i=2,4
        found=.false.
        do j=1,4
           idiff=abs(itab(j,1,1)-list(i-1,1))+abs(itab(j,1,2)-list(i-1,2))
           if (.not.done(j).and.idiff==1) then
              list(i,:)=itab(j,1,:)
              done(j)=.true.
              found=.true.
              exit
           endif
        enddo
        if (.not.found) call prterr('Sort 2d error')
     enddo

     do i=1,4  ! run over itab
        do k=1,4
           if (all(itab(i,1,:)==list(k,:))) then
              rlist(k,:)=rtab(i,:)
           endif
        enddo
     enddo

     r12=(rlist(2,:)-rlist(1,:)+rlist(3,:)-rlist(4,:))/2
     r23=(rlist(3,:)-rlist(2,:)+rlist(4,:)-rlist(1,:))/2

     if (r12(1)*r23(2)-r23(1)*r12(2) < 0) then
        list2(1,:)=list(1,:)
        list2(2,:)=list(4,:)
        list2(3,:)=list(3,:)
        list2(4,:)=list(2,:)

        list=list2
     endif

     itabsrt=list

     do i=1,4
        found=.false.
        do j=1,4
           if (all(itab(j,1,:)==list(i,:))) then
              rtabsrt(i,:)=rtab(j,:)
              found=.true.
              exit
           endif
        enddo
        if (.not.found) call prterr('Sort 2d error 2')
     enddo

  else if (NDIM==3) then

! ----------------------------------------------------------------------
!     Make copy of itab for use of function equal
!     itabcopy is global in this module
! ----------------------------------------------------------------------

  do i=1,8
     do j=1,4
        do k=1,3
           itabcopy(i,j,k)=itab(i,j,k)
        enddo
     enddo
  enddo

! ----------------------------------------------------------------------
!     First list member is first node in table
!     second is its first connection
!     third is second's connection which does not go back on itself
! ----------------------------------------------------------------------

  list(1,1)=itab(1,1,1)
  list(1,2)=itab(1,1,2)
  list(1,3)=itab(1,1,3)
  
  next(1)=itab(1,2,1)
  next(2)=itab(1,2,2)
  next(3)=itab(1,2,3)
  
  do i=1,8
     if (equal(i,1,next(1),next(2),next(3))) then
        list(2,1)=itab(i,1,1)
        list(2,2)=itab(i,1,2)
        list(2,3)=itab(i,1,3)
        if (.not.equal(i,2,list(1,1),list(1,2),list(1,3))) then
           next(1)=itab(i,2,1)
           next(2)=itab(i,2,2)
           next(3)=itab(i,2,3)
        else
           next(1)=itab(i,3,1)
           next(2)=itab(i,3,2)
           next(3)=itab(i,3,3)
        endif
        goto 1
     endif
  enddo

1 continue

! ----------------------------------------------------------------------
!     4th list member is connected to 3rd member and also back to first
!     member
! ----------------------------------------------------------------------

  do i=1,8
     if (equal(i,1,next(1),next(2),next(3))) then
        list(3,1)=itab(i,1,1)
        list(3,2)=itab(i,1,2)
        list(3,3)=itab(i,1,3)
        do j=2,4
           if (.not.equal(i,j,list(2,1),list(2,2),list(2,3))) then
              itest(1)=itab(i,j,1)
              itest(2)=itab(i,j,2)
              itest(3)=itab(i,j,3)
              do i1=1,8
                 if (equal(i1,1,itest(1),itest(2),itest(3))) then
                    do j1=2,4
                       if (equal(i1,j1,list(1,1),list(1,2),list(1,3))) then
                          list(4,1)=itab(i1,1,1)
                          list(4,2)=itab(i1,1,2)
                          list(4,3)=itab(i1,1,3)
                       endif
                    enddo
                 endif
              enddo
           endif
        enddo
     endif
  enddo

! ----------------------------------------------------------------------
!     prepare opposite list, node is opposite if it connects to node
!     but is not in list
! ----------------------------------------------------------------------
      
  do i=1,4                  ! scan list
     ipos=-1
     do ip=1,8
        if (equal(ip,1,list(i,1),list(i,2),list(i,3))) ipos=ip
     enddo
     
     do j=2,4               
        itest(1)=itab(ipos,j,1) ! scan connection from listed node
        itest(2)=itab(ipos,j,2) 
        itest(3)=itab(ipos,j,3) 
        ifail=0
        do k=1,4            ! scan list
           if (itest(1).eq.list(k,1).and.&
               itest(2).eq.list(k,2).and.&
               itest(3).eq.list(k,3)) ifail=1
        enddo
        if (ifail.eq.0) then
           listop(i,1)=itest(1)
           listop(i,2)=itest(2)
           listop(i,3)=itest(3)
        endif
     enddo
  enddo

! ----------------------------------------------------------------------
!     Prepare position vectors of list'ed nodes, rlist, rlistop
! ----------------------------------------------------------------------

  do i=1,8
     itest(1)=itab(i,1,1)
     itest(2)=itab(i,1,2)
     itest(3)=itab(i,1,3)
     
     do k=1,4
        if ( itest(1).eq.list(k,1).and.&
             itest(2).eq.list(k,2).and.&
             itest(3).eq.list(k,3)) then
           do j=1,3
              rlist(k,j)=rtab(i,j)
           enddo
        endif
     enddo

     do k=1,4
        if ( itest(1).eq.listop(k,1).and.&
             itest(2).eq.listop(k,2).and.&
             itest(3).eq.listop(k,3)) then
           do j=1,3
              rlistop(k,j)=rtab(i,j)
           enddo
        endif
     enddo
  enddo

! ----------------------------------------------------------------------
!     prepare vectors r12=(r12+r43)/2, r23=(r23+r14)/2
!     s=r12 x r23 should be parallel to t
! ----------------------------------------------------------------------

  do i=1,3
     r12(i)=(rlist(2,i)-rlist(1,i)+rlist(3,i)-rlist(4,i))/2
     r23(i)=(rlist(3,i)-rlist(2,i)+rlist(4,i)-rlist(1,i))/2
  enddo

  s(1)=r12(2)*r23(3)-r23(2)*r12(3)
  s(2)=-(r12(1)*r23(3)-r23(1)*r12(3))
  s(3)=r12(1)*r23(2)-r23(1)*r12(2)

! ----------------------------------------------------------------------
!     Prepare centroid of list face and opposite faces, t is difference
!     from face to opposite
! ----------------------------------------------------------------------

  do i=1,3
     face(i)=0
     faceop(i)=0
     do j=1,4
        face(i)=face(i)+rlist(j,i)
        faceop(i)=faceop(i)+rlistop(j,i)
     enddo
     face(i)=face(i)/4
     faceop(i)=faceop(i)/4
  enddo

  do i=1,3
     t(i)=faceop(i)-face(i)
  enddo

! ----------------------------------------------------------------------
!     dp=s.t, if positive ok, else reverse both lists
! ----------------------------------------------------------------------

  dp=0
  do i=1,3
     dp=dp+s(i)*t(i)
  enddo

  if (dp.lt.0) then
     do j=1,3
        list2(1,j)=list(1,j)
        list2(2,j)=list(4,j)
        list2(3,j)=list(3,j)
        list2(4,j)=list(2,j)

        listop2(1,j)=listop(1,j)
        listop2(2,j)=listop(4,j)
        listop2(3,j)=listop(3,j)
        listop2(4,j)=listop(2,j)
     enddo

     do i=1,4
        do j=1,3
           list(i,j)=list2(i,j)
           listop(i,j)=listop2(i,j)
        enddo
     enddo
  endif

! ----------------------------------------------------------------------
!     append listop to list
! ----------------------------------------------------------------------

  do i=1,4
     do j=1,3
        list(i+4,j)=listop(i,j)
     enddo
  enddo
  
! ----------------------------------------------------------------------
!     prepare itabsrt which is sorted in local index order
! ----------------------------------------------------------------------

  do i=1,8
     do j=1,3
        itabsrt(i,j)=list(i,j)
     enddo
  enddo

  ind=0
  do i=1,8
     do j=1,8
        if (equal(j,1,list(i,1),list(i,2),list(i,3))) then
           ind=ind+1
           do k=1,3
              rtabsrt(ind,k)=rtab(j,k)
           enddo
        endif
     enddo
  enddo

  endif
  
end subroutine sort
      
! ######################################################################

logical function equal(itabr,itabc,i,j,k)

! ----------------------------------------------------------------------
!     Compare itab(itabr,itabc,1-3) with i,j,k resp
! ----------------------------------------------------------------------

  implicit none
  integer itabr,itabc,i,j,k
  
  equal=itabcopy(itabr,itabc,1).eq.i.and.&
        itabcopy(itabr,itabc,2).eq.j.and.&
        itabcopy(itabr,itabc,3).eq.k

end function equal

! ######################################################################

subroutine jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet)

! ----------------------------------------------------------------------
!     Calculate jacobian matrix its inv and det at point xi,eta,mu.
!     returns jac,invjac,jacdet through common (does not read from common)
! ----------------------------------------------------------------------

  use common, only : NDIM
  use shape
  use iofile, only : prterr
  
  double precision xi,eta,mu
  double precision rtabsrt(2**NDIM,NDIM),jac(NDIM,NDIM),invjac(NDIM,NDIM),jacdet
  integer mlo(3),mhi(3)
  double precision Cof(3,3)
  integer i,j,l,ilo,ihi,jlo,jhi

  save mlo,mhi
  
  data mlo/2,1,1/
  data mhi/3,3,2/

  do j=1,NDIM
     jac(1,j)=0
     do l=1,2**NDIM
        jac(1,j)=jac(1,j)+rtabsrt(l,j)*dNdxi(l,xi,eta,mu)
     enddo
  enddo

  if (NDIM>=2) then
     do j=1,NDIM
        jac(2,j)=0
        do l=1,2**NDIM
           jac(2,j)=jac(2,j)+rtabsrt(l,j)*dNdeta(l,xi,eta,mu)
        enddo
     enddo
  endif

  if (NDIM==3) then
     do j=1,NDIM
        jac(3,j)=0
        do l=1,2**NDIM
           jac(3,j)=jac(3,j)+rtabsrt(l,j)*dNdmu(l,xi,eta,mu)
        enddo
     enddo
  endif

! ----------------------------------------------------------------------
!     Calculate cofactor matrix. (-1)**(i+j)=pm(i,j)
! ----------------------------------------------------------------------

  if (NDIM==3) then

     do i=1,3
        do j=1,3
           ilo=mlo(i)
           ihi=mhi(i)
           
           jlo=mlo(j)
           jhi=mhi(j)
           
           Cof(i,j)=(jac(ilo,jlo)*jac(ihi,jhi)-jac(ilo,jhi)*jac(ihi,jlo)) &
                *(-1)**(i+j)
        enddo
     enddo

     jacdet=jac(1,1)*Cof(1,1)+jac(1,2)*Cof(1,2)+jac(1,3)*Cof(1,3)
     if (jacdet.lt.0) call prterr('Error, negative Jacobian')

     do i=1,3
        do j=1,3
           invjac(i,j)=Cof(j,i)/jacdet
        enddo
     enddo
  else if (NDIM==2) then

     jacdet=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)

     if (jacdet<0) call prterr('Error, negative Jacobian 2d')

     invjac(1,1)= jac(2,2)/jacdet
     invjac(2,2)= jac(1,1)/jacdet
     invjac(1,2)=-jac(1,2)/jacdet
     invjac(2,1)=-jac(2,1)/jacdet
     
  else
     invjac(1,1)=1/jac(1,1)
     jacdet=jac(1,1)

     if (jacdet<0) call prterr('Error, negative Jacobian 1d')
  endif

end subroutine jacobian

! ######################################################################

subroutine getdN(l,xi,eta,mu,invjac,dN)

! ----------------------------------------------------------------------
!     Calculate vector of dN=(dNl/dx,dNl/dy,dNl/dz) at point (xi,eta,mu)
! ----------------------------------------------------------------------

  use common, only : NDIM
  use shape

  integer l,i,j
  double precision xi,eta,mu,invjac(NDIM,NDIM),dN(NDIM)
  double precision vec(NDIM)
  
               vec(1)=dNdxi(l,xi,eta,mu)
  if (NDIM>=2) vec(2)=dNdeta(l,xi,eta,mu)
  if (NDIM==3) vec(3)=dNdmu(l,xi,eta,mu)
  
  do i=1,NDIM
     dN(i)=0
     do j=1,NDIM
        dN(i)=dN(i)+invjac(i,j)*vec(j)
     enddo
  enddo

end subroutine getdN

! ######################################################################

function centre(i,j,k)
  
  use common, only : NDIM
  use common, only : rnode
  
  integer i,j,k
  double precision, dimension(NDIM) :: r1,r2,r3,r4,r5,r6,r7,r8,centre
  
  if (NDIM==3) then
     r1=rnode(i,  j,  k,:)
     r2=rnode(i+1,j,  k,:)
     r3=rnode(i+1,j+1,k,:)
     r4=rnode(i,  j+1,k,:)
     
     r5=rnode(i,  j,  k+1,:)
     r6=rnode(i+1,j,  k+1,:)
     r7=rnode(i+1,j+1,k+1,:)
     r8=rnode(i,  j+1,k+1,:)
  
     centre=(r1+r2+r3+r4+r5+r6+r7+r8)/8
  else if (NDIM==2) then
     r1=rnode(i,  j,  0,:)
     r2=rnode(i+1,j,  0,:)
     r3=rnode(i+1,j+1,0,:)
     r4=rnode(i,  j+1,0,:)
     centre=(r1+r2+r3+r4)/4
  else
     r1=rnode(i,  0,  0,:)
     r2=rnode(i+1,0,  0,:)
     centre=(r1+r2)/2
  endif

end function centre

end module geom
