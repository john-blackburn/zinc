module cuboid

  implicit none

  private
  public getQcuboid,getQcuboid2,getRcuboid,getJacobianCuboid

  integer, parameter :: xbar(8)=[-1, 1, 1,-1,  -1, 1, 1,-1]
  integer, parameter :: ybar(8)=[-1,-1, 1, 1,  -1,-1, 1, 1]
  integer, parameter :: zbar(8)=[-1,-1,-1,-1,   1, 1, 1, 1]

  double precision g_rtabsrt(8,3)
  double precision g_xmin,g_xmax, g_ymin,g_ymax, g_zmin,g_zmax

  interface Ix
     module procedure Ix0,Ix1,Ix2,Ix3
  end interface Ix

  interface Iy
     module procedure Iy0,Iy1,Iy2,Iy3
  end interface Iy

  interface Iz
     module procedure Iz0,Iz1,Iz2,Iz3
  end interface Iz

contains

  subroutine getQcuboid(rmFixed)

    use common
    use matrices
    use indexq
    use iofile

    logical rmFixed,inin,jnin,knin,iein,jein,kein,found
    integer i,j,k,in,jn,kn,iep,ie,je,ke,iCCx
    integer L1,L2,ip,inode,jnode,knode,i1,j1,k1,irege,itabsrt(8,3)
    double precision xmax,ymax,zmax,xmin,ymin,zmin,dx,dy,dz,rtabsrt(8,3),xL1,yL1,zL1,xL2,yL2,zL2
    double precision rc(3),elint(3,3),tot,ur(nvar),dur(nvar,3)
    integer ii,kk,ialpha,ibeta,jj,ll,lenQst,L

    integer, parameter :: ies(8)=[ 1,-1, 1,-1, 1,-1, 1,-1]
    integer, parameter :: jes(8)=[ 1, 1,-1,-1, 1, 1,-1,-1]
    integer, parameter :: kes(8)=[ 1, 1, 1, 1,-1,-1,-1,-1]
    integer, parameter :: x=1, y=2, z=3

    if (rmFixed) call prterr('getQcuboid: removeFixed not supported')

    print *,'Calculating Q matrix (assume all elements are cuboids)'
    lenQ=0

    do k=0,kmax
    write (*,*) 'k-plane:',k,'/',kmax
    do j=0,jmax
    do i=0,imax

       do in=i-1,i+1
       do jn=j-1,j+1
       do kn=k-1,k+1

          inin=(in.ge.0.and.in.le.imax)
          jnin=(jn.ge.0.and.jn.le.jmax)
          knin=(kn.ge.0.and.kn.le.kmax)

          if (inin.and.jnin.and.knin) then
        
             lenQst=lenQ+1
         
             do iep=1,8
                ie=ies(iep)
                je=jes(iep)
                ke=kes(iep)
           
                iein=(i+ie.ge.0.and.i+ie.le.imax)
                jein=(j+je.ge.0.and.j+je.le.jmax)
                kein=(k+ke.ge.0.and.k+ke.le.kmax)
               
                if (iein.and.jein.and.kein) then

                   i1=i+(ie-1)/2
                   j1=j+(je-1)/2
                   k1=k+(ke-1)/2            ! for 2D, k=0, ke=1 => k1=0
                   irege=iregup(i1,j1,k1)
                   
                   itabsrt(1,:)=[i1,    j1,  k1]
                   itabsrt(2,:)=[i1+1,  j1,  k1]
                   itabsrt(3,:)=[i1+1,j1+1,  k1]
                   itabsrt(4,:)=[i1  ,j1+1,  k1]
                   
                   itabsrt(5,:)=[i1,    j1,k1+1]
                   itabsrt(6,:)=[i1+1,  j1,k1+1]
                   itabsrt(7,:)=[i1+1,j1+1,k1+1]
                   itabsrt(8,:)=[i1  ,j1+1,k1+1]

                   found=.false.
                   do ip=1,8
                      if (itabsrt(ip,1).eq.in.and.itabsrt(ip,2).eq.jn.and.itabsrt(ip,3).eq.kn) found=.true.
                   enddo

                   if (found) then

                      L1=0; L2=0
                      do ip=1,8
                         inode=itabsrt(ip,1)
                         jnode=itabsrt(ip,2)
                         knode=itabsrt(ip,3)
                         if (inode.eq.i.and.jnode.eq.j.and.knode.eq.k) L1=ip
                         if (inode.eq.in.and.jnode.eq.jn.and.knode.eq.kn) L2=ip
                      enddo
                      
                      if (L1==0.or.L2==0) call prterr('Error in getQ, bad L1,L2')
                      
                      xmax=rnode(i1+1,j1,k1,1); xmin=rnode(i1,j1,k1,1)
                      ymax=rnode(i1,j1+1,k1,2); ymin=rnode(i1,j1,k1,2)
                      zmax=rnode(i1,j1,k1+1,3); zmin=rnode(i1,j1,k1,3)

                      dx=xmax-xmin
                      dy=ymax-ymin
                      dz=zmax-zmin
                        
                      do L=1,8
                         rtabsrt(L,:)=rnode(itabsrt(L,1),itabsrt(L,2),itabsrt(L,3),:)
                      enddo
                      
                      rc=0
                      do L=1,8
                         rc=rc+rtabsrt(L,:)
                      enddo
                      rc=rc/8
                      
                      call getelinfo(itabsrt,rtabsrt,ur,dur)

                      xL1=rtabsrt(L1,1)       ! rnode(i,j,k,1)
                      yL1=rtabsrt(L1,2)
                      zL1=rtabsrt(L1,3)

                      xL2=rtabsrt(L2,1)       ! rnode(in,jn,kn,1)
                      yL2=rtabsrt(L2,2)
                      zL2=rtabsrt(L2,3)

!                      print *,xL1,yL1,zL1,'?=',rnode(i,j,k,:)
!                      print *,xL2,yL2,zL2,'?=',rnode(in,jn,kn,:)
!                      stop

!   double precision g_xmin,g_xmax, g_ymin,g_ymax, g_zmin,g_zmax

                   elint(x,x)=xbar(L1)*xbar(L2)/dx &
                        *dintf(ymin, ymax, ybar(L1)/dy, yL1, ybar(L2)/dy, yL2) &
                        *dintf(zmin, zmax, zbar(L1)/dz, zL1, zbar(L2)/dz, zL2)

                   elint(y,y)=ybar(L1)*ybar(L2)/dy &
                        *dintf(xmin, xmax, xbar(L1)/dx, xL1, xbar(L2)/dx, xL2) &
                        *dintf(zmin, zmax, zbar(L1)/dz, zL1, zbar(L2)/dz, zL2)

                   elint(z,z)=zbar(L1)*zbar(L2)/dz &
                        *dintf(xmin, xmax, xbar(L1)/dx, xL1, xbar(L2)/dx, xL2) &
                        *dintf(ymin, ymax, ybar(L1)/dy, yL1, ybar(L2)/dy, yL2)

!----------------------------------------------------------------------

                   elint(x,y)=xbar(L1)*ybar(L2)/(dx*dy) &
                        *(xmax+xbar(L2)*(xmax-xL2)**2/(2*dx)-xmin-xbar(L2)*(xmin-xL2)**2/(2*dx)) &
                        *(ymax+ybar(L1)*(ymax-yL1)**2/(2*dy)-ymin-ybar(L1)*(ymin-yL1)**2/(2*dy)) &
                        *dintf(zmin, zmax, zbar(l1)/dz, zL1, zbar(L2)/dz, zL2)

                   elint(y,x)=ybar(L1)*xbar(L2)/(dy*dx) &
                        *(ymax+ybar(L2)*(ymax-yL2)**2/(2*dy)-ymin-ybar(L2)*(ymin-yL2)**2/(2*dy)) &
                        *(xmax+xbar(L1)*(xmax-xL1)**2/(2*dx)-xmin-xbar(L1)*(xmin-xL1)**2/(2*dx)) &
                        *dintf(zmin, zmax, zbar(l1)/dz, zL1, zbar(L2)/dz, zL2)

!----------------------------------------------------------------------

                   elint(x,z)=xbar(L1)*zbar(L2)/(dx*dz) &
                        *(xmax+xbar(L2)*(xmax-xL2)**2/(2*dx)-xmin-xbar(L2)*(xmin-xL2)**2/(2*dx)) &
                        *(zmax+zbar(L1)*(zmax-zL1)**2/(2*dz)-zmin-zbar(L1)*(zmin-zL1)**2/(2*dz)) &
                        *dintf(ymin, ymax, ybar(L1)/dy, yL1, ybar(L2)/dy, yL2)

                   elint(z,x)=zbar(L1)*xbar(L2)/(dz*dx) &
                        *(zmax+zbar(L2)*(zmax-zL2)**2/(2*dz)-zmin-zbar(L2)*(zmin-zL2)**2/(2*dz)) &
                        *(xmax+xbar(L1)*(xmax-xL1)**2/(2*dx)-xmin-xbar(L1)*(xmin-xL1)**2/(2*dx)) &
                        *dintf(ymin, ymax, ybar(L1)/dy, yL1, ybar(L2)/dy, yL2)

!----------------------------------------------------------------------

                   elint(y,z)=ybar(L1)*zbar(L2)/(dy*dz) &
                        *(ymax+ybar(L2)*(ymax-yL2)**2/(2*dy)-ymin-ybar(L2)*(ymin-yL2)**2/(2*dy)) &
                        *(zmax+zbar(L1)*(zmax-zL1)**2/(2*dz)-zmin-zbar(L1)*(zmin-zL1)**2/(2*dz)) &
                        *dintf(xmin, xmax, xbar(L1)/dx, xL1, xbar(L2)/dx, xL2)

                   elint(z,y)=zbar(L1)*ybar(L2)/(dz*dy) &
                        *(zmax+zbar(L2)*(zmax-zL2)**2/(2*dz)-zmin-zbar(L2)*(zmin-zL2)**2/(2*dz)) &
                        *(ymax+ybar(L1)*(ymax-yL1)**2/(2*dy)-ymin-ybar(L1)*(ymin-yL1)**2/(2*dy)) &
                        *dintf(xmin, xmax, xbar(L1)/dx, xL1, xbar(L2)/dx, xL2)

                      do ii=1,nvar
                         do kk=1,nvar
                            ialpha=indQ(i,j,k,ii)
                            ibeta =indQ(in,jn,kn,kk)

                            tot=0
                            do jj=1,3
                               do ll=1,3
                                  tot=tot+CCval(irege,ii,jj,kk,ll,ur,dur,rc,iCCx)*elint(jj,ll)
                               enddo
                            enddo

!                         tot=tot+aaval(irege,ii,kk,ur,dur,rc,iaax)*elintNN

                            if (tot /= 0.or.newton=="YES") then
                               found=.false.
                               do ip=lenQst,lenQ   ! prev lenQst,lenQ
                                  if (iQ(ip) == ialpha.and.jQ(ip) == ibeta) then
                                     Qval(ip)=Qval(ip)+tot
                                     found=.true.
                                     exit
                                  endif
                               enddo
                            
                               if (.not.found) then                 
                                  lenQ=lenQ+1
                                  iQ(lenQ)=ialpha
                                  jQ(lenQ)=ibeta
                                  Qval(lenQ)=tot
                               endif
                            endif

                         enddo    ! var number ii
                      enddo

                   endif  ! found=true
                endif     ! iein
             enddo        ! iep
          endif           ! inin
       enddo
       enddo
       enddo

    enddo
    enddo
    enddo
  end subroutine getQcuboid

! ######################################################################

  subroutine getQcuboid2(rmFixed)

    use common
    use matrices
    use indexq
    use iofile

    logical rmFixed
    double precision tot,rtabsrt(8,3)
    double precision xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz,rc(3),ur(nvar),dur(nvar,3)
    double precision elint(3,3),xL1,yL1,zL1,xL2,yL2,zL2
    integer, parameter :: x=1, y=2, z=3
    integer i,j,k,irege,L,L1,L2,ii,kk,jj,ll,ip,ialpha,ibeta,iCCx
    integer itabsrt(8,3),iL1,jL1,kL1,iL2,jL2,kL2,endprevrow,lenQst
    logical found

    if (rmFixed) call prterr('getQcuboid: removeFixed not supported')

    print *,'Calculating Q-matrix (assume all elements are cuboids)'

    lenQ=0
    lenQst=1

    do i=0,imax-1
       write (*,*) 'i-plane:',i,'/',imax

       do j=0,jmax-1
          do k=0,kmax-1

             print *,i,j,k,lenQ,lenQst
             
             irege=iregup(i,j,k)

             itabsrt(1,:)=[i,    j,  k]
             itabsrt(2,:)=[i+1,  j,  k]
             itabsrt(3,:)=[i+1,j+1,  k]
             itabsrt(4,:)=[i  ,j+1,  k]

             itabsrt(5,:)=[i,    j,k+1]
             itabsrt(6,:)=[i+1,  j,k+1]
             itabsrt(7,:)=[i+1,j+1,k+1]
             itabsrt(8,:)=[i  ,j+1,k+1]

             xmax=rnode(i+1,j,k,1); xmin=rnode(i,j,k,1)
             ymax=rnode(i,j+1,k,2); ymin=rnode(i,j,k,2)
             zmax=rnode(i,j,k+1,3); zmin=rnode(i,j,k,3)

             dx=xmax-xmin
             dy=ymax-ymin
             dz=zmax-zmin
             
             do L=1,8
                rtabsrt(L,:)=rnode(itabsrt(L,1),itabsrt(L,2),itabsrt(L,3),:)
             enddo

             rc=0
             do L=1,8
                rc=rc+rtabsrt(L,:)
             enddo
             rc=rc/8

             call getelinfo(itabsrt,rtabsrt,ur,dur)

             do L1=1,8
                do L2=1,8

                   iL1=itabsrt(L1,1)
                   jL1=itabsrt(L1,2)
                   kL1=itabsrt(L1,3)

                   iL2=itabsrt(L2,1)
                   jL2=itabsrt(L2,2)
                   kL2=itabsrt(L2,3)

                   xL1=rtabsrt(L1,1)
                   yL1=rtabsrt(L1,2)
                   zL1=rtabsrt(L1,3)

                   xL2=rtabsrt(L2,1)
                   yL2=rtabsrt(L2,2)
                   zL2=rtabsrt(L2,3)

                   elint(x,x)=xbar(L1)*xbar(L2)/dx &
                        *dintf(ymin, ymax, ybar(L1)/dy, yL1, ybar(L2)/dy, yL2) &
                        *dintf(zmin, zmax, zbar(L1)/dz, zL1, zbar(L2)/dz, zL2)

                   elint(y,y)=ybar(L1)*ybar(L2)/dy &
                        *dintf(xmin, xmax, xbar(L1)/dx, xL1, xbar(L2)/dx, xL2) &
                        *dintf(zmin, zmax, zbar(L1)/dz, zL1, zbar(L2)/dz, zL2)

                   elint(z,z)=zbar(L1)*zbar(L2)/dz &
                        *dintf(xmin, xmax, xbar(L1)/dx, xL1, xbar(L2)/dx, xL2) &
                        *dintf(ymin, ymax, ybar(L1)/dy, yL1, ybar(L2)/dy, yL2)

!----------------------------------------------------------------------

                   elint(x,y)=xbar(L1)*ybar(L2)/(dx*dy) &
                        *(xmax+xbar(L2)*(xmax-xL2)**2/(2*dx)-xmin-xbar(L2)*(xmin-xL2)**2/(2*dx)) &
                        *(ymax+ybar(L1)*(ymax-yL1)**2/(2*dy)-ymin-ybar(L1)*(ymin-yL1)**2/(2*dy)) &
                        *dintf(zmin, zmax, zbar(l1)/dz, zL1, zbar(L2)/dz, zL2)

                   elint(y,x)=ybar(L1)*xbar(L2)/(dy*dx) &
                        *(ymax+ybar(L2)*(ymax-yL2)**2/(2*dy)-ymin-ybar(L2)*(ymin-yL2)**2/(2*dy)) &
                        *(xmax+xbar(L1)*(xmax-xL1)**2/(2*dx)-xmin-xbar(L1)*(xmin-xL1)**2/(2*dx)) &
                        *dintf(zmin, zmax, zbar(l1)/dz, zL1, zbar(L2)/dz, zL2)

!----------------------------------------------------------------------

                   elint(x,z)=xbar(L1)*zbar(L2)/(dx*dz) &
                        *(xmax+xbar(L2)*(xmax-xL2)**2/(2*dx)-xmin-xbar(L2)*(xmin-xL2)**2/(2*dx)) &
                        *(zmax+zbar(L1)*(zmax-zL1)**2/(2*dz)-zmin-zbar(L1)*(zmin-zL1)**2/(2*dz)) &
                        *dintf(ymin, ymax, ybar(L1)/dy, yL1, ybar(L2)/dy, yL2)

                   elint(z,x)=zbar(L1)*xbar(L2)/(dz*dx) &
                        *(zmax+zbar(L2)*(zmax-zL2)**2/(2*dz)-zmin-zbar(L2)*(zmin-zL2)**2/(2*dz)) &
                        *(xmax+xbar(L1)*(xmax-xL1)**2/(2*dx)-xmin-xbar(L1)*(xmin-xL1)**2/(2*dx)) &
                        *dintf(ymin, ymax, ybar(L1)/dy, yL1, ybar(L2)/dy, yL2)

!----------------------------------------------------------------------

                   elint(y,z)=ybar(L1)*zbar(L2)/(dy*dz) &
                        *(ymax+ybar(L2)*(ymax-yL2)**2/(2*dy)-ymin-ybar(L2)*(ymin-yL2)**2/(2*dy)) &
                        *(zmax+zbar(L1)*(zmax-zL1)**2/(2*dz)-zmin-zbar(L1)*(zmin-zL1)**2/(2*dz)) &
                        *dintf(xmin, xmax, xbar(L1)/dx, xL1, xbar(L2)/dx, xL2)

                   elint(z,y)=zbar(L1)*ybar(L2)/(dz*dy) &
                        *(zmax+zbar(L2)*(zmax-zL2)**2/(2*dz)-zmin-zbar(L2)*(zmin-zL2)**2/(2*dz)) &
                        *(ymax+ybar(L1)*(ymax-yL1)**2/(2*dy)-ymin-ybar(L1)*(ymin-yL1)**2/(2*dy)) &
                        *dintf(xmin, xmax, xbar(L1)/dx, xL1, xbar(L2)/dx, xL2)

                   do ii=1,nvar
                      do kk=1,nvar
                         ialpha=indQ(iL1,jL1,kL1,ii)
                         ibeta =indQ(iL2,jL2,kL2,kk)

                         tot=0
                         do jj=1,3
                            do ll=1,3
                               tot=tot+CCval(irege,ii,jj,kk,ll,ur,dur,rc,iCCx)*elint(jj,ll)
                            enddo
                         enddo

!                         print *,"tot=",L1,L2,ii,kk,tot

!                         tot=tot+aaval(irege,ii,kk,ur,dur,rc,iaax)*elintNN

                         if (tot /= 0.or.newton=="YES") then
                            found=.false.
                            do ip=lenQst,lenQ   ! prev lenQst,lenQ
                               if (iQ(ip) == ialpha.and.jQ(ip) == ibeta) then
                                  Qval(ip)=Qval(ip)+tot
                                  found=.true.
                                  exit
                               endif
                            enddo
                            
                            if (.not.found) then                 
                               lenQ=lenQ+1
                               iQ(lenQ)=ialpha
                               jQ(lenQ)=ibeta
                               Qval(lenQ)=tot
                            endif
                         endif

                      enddo    ! var number
                   enddo

                enddo          ! local numbers L1,L2
             enddo

          enddo                ! i,j,k
       enddo

       if (i /= 0) then
          lenQst=endprevrow+1
       endif

       endprevrow=lenQ

    enddo

  end subroutine getQcuboid2

! ######################################################################

  subroutine getRcuboid(append)

! ----------------------------------------------------------------------
! Form the R vector. Uses rnode, ff, iff, u (non-linear)
! ----------------------------------------------------------------------

    use common
    use geom
    use shape
    use indexQ
    use matrices
    use util
    use iofile, only : prterr
    
    logical append
    
    integer itabsrt(2**NDIM,NDIM)
    
    double precision rtabsrt(2**NDIM,NDIM),elintN
    double precision ur(nvar),dur(nvar,NDIM),rc(NDIM)
    double precision xmax,ymax,zmax,xmin,ymin,zmin,dx,dy,dz,xL1,yL1,zL1
    
    integer i,j,k,iL1,jL1,kL1,irege,ii,L1,ialpha,L

! ----------------------------------------------------------------------
! Loop over nodes. getQ must have been called already so RR has non-zero
! entries by now
! ----------------------------------------------------------------------

    print *,'Calculating R matrix'

    if (.not.append) RR=0            ! start from scratch

    do k=0,kmax-1
    do j=0,jmax-1
    do i=0,imax-1

       irege=iregup(i,j,k)

       itabsrt(1,:)=[i,    j,  k]
       itabsrt(2,:)=[i+1,  j,  k]
       itabsrt(3,:)=[i+1,j+1,  k]
       itabsrt(4,:)=[i  ,j+1,  k]
       
       itabsrt(5,:)=[i,    j,k+1]
       itabsrt(6,:)=[i+1,  j,k+1]
       itabsrt(7,:)=[i+1,j+1,k+1]
       itabsrt(8,:)=[i  ,j+1,k+1]
       
       xmax=rnode(i+1,j,k,1); xmin=rnode(i,j,k,1)
       ymax=rnode(i,j+1,k,2); ymin=rnode(i,j,k,2)
       zmax=rnode(i,j,k+1,3); zmin=rnode(i,j,k,3)
       
       dx=xmax-xmin
       dy=ymax-ymin
       dz=zmax-zmin
        
       do L=1,8
          rtabsrt(L,:)=rnode(itabsrt(L,1),itabsrt(L,2),itabsrt(L,3),:)
       enddo

       rc=0
       do L=1,8
          rc=rc+rtabsrt(L,:)
       enddo
       rc=rc/8
       
       call getelinfo(itabsrt,rtabsrt,ur,dur)

       do L1=1,8

          iL1=itabsrt(L1,1)
          jL1=itabsrt(L1,2)
          kL1=itabsrt(L1,3)
          
          xL1=rtabsrt(L1,1)
          yL1=rtabsrt(L1,2)
          zL1=rtabsrt(L1,3)
          
! ----------------------------------------------------------------------
! Calculate elintN = int N_alpha^(e) dV, alpha -> l1
! for given node and element
! ----------------------------------------------------------------------

          elintN=(xmax+xbar(L1)*(xmax-xL1)**2/(2*dx)-xmin-xbar(L1)*(xmin-xL1)**2/(2*dx)) &
                *(ymax+ybar(L1)*(ymax-yL1)**2/(2*dy)-ymin-ybar(L1)*(ymin-yL1)**2/(2*dy)) &
                *(zmax+zbar(L1)*(zmax-zL1)**2/(2*dz)-zmin-zbar(L1)*(zmin-zL1)**2/(2*dz))
          
          do ii=1,nvar
             ialpha=indQ(iL1,jL1,kL1,ii)
             RR(ialpha)=RR(ialpha)+ffval(irege,ii,ur,dur,rc)*elintN
          enddo
       enddo

    enddo
    enddo
    enddo

  end subroutine getRcuboid

! ######################################################################

  function intf(s,a,b,c,d)
    double precision intf,s,a,b,c,d
    intf=s+a*(s-b)**2/2+c*(s-d)**2/2+a*c*(s**3/3-s**2/2*(b+d)+s*b*d)
  end function intf

! ######################################################################

  function dintf(smin,smax,a,b,c,d)
    double precision dintf,smin,smax,a,b,c,d
    dintf=intf(smax,a,b,c,d)-intf(smin,a,b,c,d)
  end function dintf

! ######################################################################

  subroutine getJacobianCuboid
    call getJacobianCuboidC
    call getJacobianCuboidF
  end subroutine getJacobianCuboid

! ######################################################################

  subroutine getJacobianCuboidC

! ----------------------------------------------------------------------
! Add components to Qval to make it into Jacobian
! Works only for cuboid elements and accounts only for nonlinear C
! ----------------------------------------------------------------------

    use common
    use matrices
    use indexq

    integer, parameter :: x=1,y=2,z=3
    integer i,j,k,irege,itabsrt(8,3),L,la,lb,lg,ip
    integer ia,ja,ka,ib,jb,kb,ig,jg,kg,ii,jj,kk,ll,nn,irow,icol
    double precision xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz,ur(nvar),dur(nvar,3)
    double precision elint2(3,3),elint3(3,3,3),elint(nvar,nvar,nvar),tot
    double precision rtabsrt(8,3),dCdu(nvar),dCddu(nvar,3),xc,yc,zc
    character(EXPRLEN) token
    logical found

    print *,'Preparing Jacobian (getJacobianCuboidC)'

    do i=0,imax-1
       print *,'i-plane: ',i,'/',imax-1
    do j=0,jmax-1
    do k=0,kmax-1

       irege=iregup(i,j,k)

       itabsrt(1,:)=[i,    j,  k]
       itabsrt(2,:)=[i+1,  j,  k]
       itabsrt(3,:)=[i+1,j+1,  k]
       itabsrt(4,:)=[i  ,j+1,  k]

       itabsrt(5,:)=[i,    j,k+1]
       itabsrt(6,:)=[i+1,  j,k+1]
       itabsrt(7,:)=[i+1,j+1,k+1]
       itabsrt(8,:)=[i  ,j+1,k+1]
       
       xmax=rnode(i+1,j,k,1); xmin=rnode(i,j,k,1)
       ymax=rnode(i,j+1,k,2); ymin=rnode(i,j,k,2)
       zmax=rnode(i,j,k+1,3); zmin=rnode(i,j,k,3)

       dx=xmax-xmin; xc=(xmax+xmin)/2
       dy=ymax-ymin; yc=(ymax+ymin)/2
       dz=zmax-zmin; zc=(zmax+zmin)/2
             
       do L=1,8
          rtabsrt(L,:)=rnode(itabsrt(L,1),itabsrt(L,2),itabsrt(L,3),:)
       enddo
       
       call getelinfo(itabsrt,rtabsrt,ur,dur)    ! get ur, dur at centre

       ! Global variables for Ix, Iy, Iz functions
       g_rtabsrt=rtabsrt
       g_xmax=xmax; g_xmin=xmin
       g_ymax=ymax; g_ymin=ymin
       g_zmax=zmax; g_zmin=zmin

       do la=1,8
       do lb=1,8
       do lg=1,8

          ia=itabsrt(la,1)
          ja=itabsrt(la,2)
          ka=itabsrt(la,3)

          ib=itabsrt(lb,1)
          jb=itabsrt(lb,2)
          kb=itabsrt(lb,3)
          
          ig=itabsrt(lg,1)
          jg=itabsrt(lg,2)
          kg=itabsrt(lg,3)

        !        a b
        ! elint2( , )=  Ix(la,lb,lg)*Iy(la,lb,lg)*Iz(la,lb,lg)
          elint2(x,x)=  Ix(      lg)*Iy(la,lb,lg)*Iz(la,lb,lg) *xbar(la)/dx *xbar(lb)/dx
          elint2(x,y)=  Ix(   lb,lg)*Iy(la,   lg)*Iz(la,lb,lg) *xbar(la)/dx *ybar(lb)/dy
          elint2(x,z)=  Ix(   lb,lg)*Iy(la,lb,lg)*Iz(la,   lg) *xbar(la)/dx *zbar(lb)/dz
    
          elint2(y,x)=  Ix(la,   lg)*Iy(   lb,lg)*Iz(la,lb,lg) *ybar(la)/dy *xbar(lb)/dx
          elint2(y,y)=  Ix(la,lb,lg)*Iy(      lg)*Iz(la,lb,lg) *ybar(la)/dy *ybar(lb)/dy
          elint2(y,z)=  Ix(la,lb,lg)*Iy(   lb,lg)*Iz(la,   lg) *ybar(la)/dy *zbar(lb)/dz
          
          elint2(z,x)=  Ix(la,   lg)*Iy(la,lb,lg)*Iz(   lb,lg) *zbar(la)/dz *xbar(lb)/dx
          elint2(z,y)=  Ix(la,lb,lg)*Iy(la,   lg)*Iz(   lb,lg) *zbar(la)/dz *ybar(lb)/dy
          elint2(z,z)=  Ix(la,lb,lg)*Iy(la,lb,lg)*Iz(      lg) *zbar(la)/dz *zbar(lb)/dz

        !        a b g
        ! elint3( , , )=Ix(la,lb,lg)*Iy(la,lb,lg)*Iz(la,lb,lg)
          elint3(x,x,x)=Ix(        )*Iy(la,lb,lg)*Iz(la,lb,lg) *xbar(la)/dx *xbar(lb)/dx *xbar(lg)/dx
          elint3(x,x,y)=Ix(      lg)*Iy(la,lb   )*Iz(la,lb,lg) *xbar(la)/dx *xbar(lb)/dx *ybar(lg)/dy
          elint3(x,x,z)=Ix(      lg)*Iy(la,lb,lg)*Iz(la,lb   ) *xbar(la)/dx *xbar(lb)/dx *zbar(lg)/dz
          
          elint3(x,y,x)=Ix(   lb   )*Iy(la,   lg)*Iz(la,lb,lg) *xbar(la)/dx *ybar(lb)/dy *xbar(lg)/dx
          elint3(x,y,y)=Ix(   lb,lg)*Iy(la      )*Iz(la,lb,lg) *xbar(la)/dx *ybar(lb)/dy *ybar(lg)/dy
          elint3(x,y,z)=Ix(   lb,lg)*Iy(la,   lg)*Iz(la,lb   ) *xbar(la)/dx *ybar(lb)/dy *zbar(lg)/dz
          
          elint3(x,z,x)=Ix(   lb   )*Iy(la,lb,lg)*Iz(la,   lg) *xbar(la)/dx *zbar(lb)/dz *xbar(lg)/dx
          elint3(x,z,y)=Ix(   lb,lg)*Iy(la,lb   )*Iz(la,   lg) *xbar(la)/dx *zbar(lb)/dz *ybar(lg)/dy
          elint3(x,z,z)=Ix(   lb,lg)*Iy(la,lb,lg)*Iz(la      ) *xbar(la)/dx *zbar(lb)/dz *zbar(lg)/dz

! ----------------------------------------------------------------------

          elint3(y,x,x)=Ix(la      )*Iy(   lb,lg)*Iz(la,lb,lg) *ybar(la)/dy *xbar(lb)/dx *xbar(lg)/dx
          elint3(y,x,y)=Ix(la   ,lg)*Iy(   lb   )*Iz(la,lb,lg) *ybar(la)/dy *xbar(lb)/dx *ybar(lg)/dy
          elint3(y,x,z)=Ix(la   ,lg)*Iy(   lb,lg)*Iz(la,lb   ) *ybar(la)/dy *xbar(lb)/dx *zbar(lg)/dz
          
          elint3(y,y,x)=Ix(la,lb   )*Iy(      lg)*Iz(la,lb,lg) *ybar(la)/dy *ybar(lb)/dy *xbar(lg)/dx
          elint3(y,y,y)=Ix(la,lb,lg)*Iy(        )*Iz(la,lb,lg) *ybar(la)/dy *ybar(lb)/dy *ybar(lg)/dy
          elint3(y,y,z)=Ix(la,lb,lg)*Iy(      lg)*Iz(la,lb   ) *ybar(la)/dy *ybar(lb)/dy *zbar(lg)/dz
          
          elint3(y,z,x)=Ix(la,lb   )*Iy(   lb,lg)*Iz(la   ,lg) *ybar(la)/dy *zbar(lb)/dz *xbar(lg)/dx
          elint3(y,z,y)=Ix(la,lb,lg)*Iy(   lb   )*Iz(la   ,lg) *ybar(la)/dy *zbar(lb)/dz *ybar(lg)/dy
          elint3(y,z,z)=Ix(la,lb,lg)*Iy(   lb,lg)*Iz(la      ) *ybar(la)/dy *zbar(lb)/dz *zbar(lg)/dz
          
! ----------------------------------------------------------------------

          elint3(z,x,x)=Ix(la      )*Iy(la,lb,lg)*Iz(   lb,lg) *zbar(la)/dz *xbar(lb)/dx *xbar(lg)/dx
          elint3(z,x,y)=Ix(la   ,lg)*Iy(la,lb   )*Iz(   lb,lg) *zbar(la)/dz *xbar(lb)/dx *ybar(lg)/dy
          elint3(z,x,z)=Ix(la   ,lg)*Iy(la,lb,lg)*Iz(   lb   ) *zbar(la)/dz *xbar(lb)/dx *zbar(lg)/dz
    
          elint3(z,y,x)=Ix(la,lb   )*Iy(la   ,lg)*Iz(   lb,lg) *zbar(la)/dz *ybar(lb)/dy *xbar(lg)/dx
          elint3(z,y,y)=Ix(la,lb,lg)*Iy(la      )*Iz(   lb,lg) *zbar(la)/dz *ybar(lb)/dy *ybar(lg)/dy
          elint3(z,y,z)=Ix(la,lb,lg)*Iy(la   ,lg)*Iz(   lb   ) *zbar(la)/dz *ybar(lb)/dy *zbar(lg)/dz
    
          elint3(z,z,x)=Ix(la,lb   )*Iy(la,lb,lg)*Iz(      lg) *zbar(la)/dz *zbar(lb)/dz *xbar(lg)/dx
          elint3(z,z,y)=Ix(la,lb,lg)*Iy(la,lb   )*Iz(      lg) *zbar(la)/dz *zbar(lb)/dz *ybar(lg)/dy
          elint3(z,z,z)=Ix(la,lb,lg)*Iy(la,lb,lg)*Iz(        ) *zbar(la)/dz *zbar(lb)/dz *zbar(lg)/dz

          elint=0
          do ii=1,nvar
          do kk=1,nvar
          do nn=1,nvar
             do jj=1,3
             do ll=1,3

                if (iCC(irege,ii,jj,kk,ll)==3) then
                   token=sCC(irege,ii,jj,kk,ll)

                   call dCfun_du (token,xc,yc,zc,ur,dur,nvar,istep, &
                        ireg,iregup,rnode,vec,imax,jmax,kmax,dCdu)

                   call dCfun_ddu(token,xc,yc,zc,ur,dur,nvar,istep, &
                        ireg,iregup,rnode,vec,imax,jmax,kmax,dCddu)

!                   dCdu=dCfun_du(token,ur,dur)
!                   dCddu=dCfun_ddu(token,ur,dur)
                   
                   elint(ii,kk,nn)=elint(ii,kk,nn)+dCdu(nn)*elint2(jj,ll) &
                        +dCddu(nn,1)*elint3(jj,ll,1) &
                        +dCddu(nn,2)*elint3(jj,ll,2) &
                        +dCddu(nn,3)*elint3(jj,ll,3)
                endif

             enddo
             enddo
          enddo
          enddo
          enddo

          do ii=1,nvar
          do nn=1,nvar

             tot=0
             do kk=1,nvar
                tot=tot+vec(indQ(ib,jb,kb,kk))*elint(ii,kk,nn)
             enddo
             
             irow=indQ(ia,ja,ka,ii)
             icol=indQ(ig,jg,kg,nn)
           
!             call addtoQ(irow,icol,tot)

             if (iunk(irow) == 0) cycle   ! no contribution to jacobian in fixed node

             if (irowst(irow) == 1.and.irowed(irow) == 0) then
                print *,'getJacobian: Found empty row, fail!'
                stop
             endif

             found=.false.
             do ip=irowst(irow),irowed(irow)
                if (jQ(ip)==icol) then
                   found=.true.
                   Qval(ip)=Qval(ip)+tot
                   exit
                endif
             enddo

             if (.not.found) then
                print *,'getJacobianCuboidC, element not found',irow,icol,tot
                stop
             endif

          enddo
          enddo

       enddo    ! la,lb,lg
       enddo
       enddo

    enddo       ! i,j,k
    enddo
    enddo

  end subroutine getJacobianCuboidC

! ######################################################################

  subroutine getJacobianCuboidF

    use common
    use matrices
    use indexq

    integer i,j,k,irege,itabsrt(8,3),L,la,lg,ip
    integer ia,ja,ka,ig,jg,kg,ii,nn,irow,icol
    double precision xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz,ur(nvar),dur(nvar,3)
    double precision elint2,elint3(3),tot,xc,yc,zc
    double precision rtabsrt(8,3),dfdu(nvar),dfddu(nvar,3)
    character(EXPRLEN) token
    logical found

    print *,'Preparing Jacobian (getJacobianCuboidF)'

    do i=0,imax-1
       print *,'i-plane: ',i,'/',imax-1
    do j=0,jmax-1
    do k=0,kmax-1

       irege=iregup(i,j,k)

       itabsrt(1,:)=[i,    j,  k]
       itabsrt(2,:)=[i+1,  j,  k]
       itabsrt(3,:)=[i+1,j+1,  k]
       itabsrt(4,:)=[i  ,j+1,  k]

       itabsrt(5,:)=[i,    j,k+1]
       itabsrt(6,:)=[i+1,  j,k+1]
       itabsrt(7,:)=[i+1,j+1,k+1]
       itabsrt(8,:)=[i  ,j+1,k+1]
       
       xmax=rnode(i+1,j,k,1); xmin=rnode(i,j,k,1)
       ymax=rnode(i,j+1,k,2); ymin=rnode(i,j,k,2)
       zmax=rnode(i,j,k+1,3); zmin=rnode(i,j,k,3)

       dx=xmax-xmin; xc=(xmax+xmin)/2
       dy=ymax-ymin; yc=(ymax+ymin)/2
       dz=zmax-zmin; zc=(zmax+zmin)/2
             
       do L=1,8
          rtabsrt(L,:)=rnode(itabsrt(L,1),itabsrt(L,2),itabsrt(L,3),:)
       enddo
       
       call getelinfo(itabsrt,rtabsrt,ur,dur)    ! get ur, dur at centre

       ! Global variables for Ix, Iy, Iz functions
       g_rtabsrt=rtabsrt
       g_xmax=xmax; g_xmin=xmin
       g_ymax=ymax; g_ymin=ymin
       g_zmax=zmax; g_zmin=zmin

       do la=1,8
       do lg=1,8

          ia=itabsrt(la,1)
          ja=itabsrt(la,2)
          ka=itabsrt(la,3)

          ig=itabsrt(lg,1)
          jg=itabsrt(lg,2)
          kg=itabsrt(lg,3)

          elint2=Ix(la,lg)*Iy(la,lg)*Iz(la,lg)

          elint3(1)=Ix(la   )*Iy(la,lg)*Iz(la,lg)*xbar(lg)/dx
          elint3(2)=Ix(la,lg)*Iy(la   )*Iz(la,lg)*ybar(lg)/dy
          elint3(3)=Ix(la,lg)*Iy(la,lg)*Iz(la   )*zbar(lg)/dz

          do ii=1,nvar
          do nn=1,nvar

             if (iff(irege,ii)==3) then
                token=sff(irege,ii)

!                dfdu=dffun_du(token,ur,dur)
!                dfddu=dffun_ddu(token,ur,dur)

                call dffun_du (token,xc,yc,zc,ur,dur,nvar,istep, &
                     ireg,iregup,rnode,vec,imax,jmax,kmax,dfdu)

                call dffun_ddu(token,xc,yc,zc,ur,dur,nvar,istep, &
                     ireg,iregup,rnode,vec,imax,jmax,kmax,dfddu)

                tot=-dfdu(nn)*elint2 &
                    -dfddu(nn,1)*elint3(1) &
                    -dfddu(nn,2)*elint3(2) &
                    -dfddu(nn,3)*elint3(3)

                irow=indQ(ia,ja,ka,ii)
                icol=indQ(ig,jg,kg,nn)
             
                if (iunk(irow)==0) cycle       ! very important

                if (irowst(irow) == 1.and.irowed(irow) == 0) then
                   print *,'getJacobian: Found empty row, fail!'
                   stop
                endif

                found=.false.
                do ip=irowst(irow),irowed(irow)
                   if (jQ(ip)==icol) then
                      found=.true.
                      Qval(ip)=Qval(ip)+tot
                      exit
                   endif
                enddo
          
                if (.not.found) then
                   print *,'getJacobianCuboidF, element not found',irow,icol,tot
                   stop
                endif
             endif

          enddo  ! nn=1,nvar
          enddo

       enddo     ! la,lg
       enddo

    enddo
    enddo
    enddo

  end subroutine getJacobianCuboidF

! ######################################################################

function Ix0()
  double precision Ix0
  Ix0=g_xmax-g_xmin
end function Ix0

function Ix1(l)
  integer l
  double precision Ix1,xl,dx

  xl=g_rtabsrt(l,1)
  dx=g_xmax-g_xmin

  Ix1=g_xmax+xbar(l)/dx*(g_xmax-xl)**2/2 - (g_xmin+xbar(l)/dx*(g_xmin-xl)**2/2)
end function Ix1

function Ix2(l1,l2)
  integer l1,l2
  double precision Ix2,xl1,xl2,dx

  xl1=g_rtabsrt(l1,1)
  xl2=g_rtabsrt(l2,1)
  dx=g_xmax-g_xmin

  Ix2=  Ix2i(g_xmax,xbar(l1),xbar(l2),xl1,xl2,dx) &
       -Ix2i(g_xmin,xbar(l1),xbar(l2),xl1,xl2,dx)
end function Ix2

function Ix3(l1,l2,l3)
  integer l1,l2,l3
  double precision Ix3,xl1,xl2,xl3,dx

  xl1=g_rtabsrt(l1,1)
  xl2=g_rtabsrt(l2,1)
  xl3=g_rtabsrt(l3,1)
  dx=g_xmax-g_xmin

  Ix3=  Ix3i(g_xmax,xbar(l1),xbar(l2),xbar(l3),xl1,xl2,xl3,dx) &
       -Ix3i(g_xmin,xbar(l1),xbar(l2),xbar(l3),xl1,xl2,xl3,dx)
end function Ix3

! ######################################################################

function Iy0()
  double precision Iy0
  Iy0=g_ymax-g_ymin
end function Iy0

function Iy1(l)
  integer l
  double precision Iy1,yl,dy

  yl=g_rtabsrt(l,2)
  dy=g_ymax-g_ymin

  Iy1=g_ymax+ybar(l)/dy*(g_ymax-yl)**2/2 - (g_ymin+ybar(l)/dy*(g_ymin-yl)**2/2)
end function Iy1

function Iy2(l1,l2)
  integer l1,l2
  double precision Iy2,yl1,yl2,dy

  yl1=g_rtabsrt(l1,2)
  yl2=g_rtabsrt(l2,2)
  dy=g_ymax-g_ymin

  Iy2=  Ix2i(g_ymax,ybar(l1),ybar(l2),yl1,yl2,dy) &
       -Ix2i(g_ymin,ybar(l1),ybar(l2),yl1,yl2,dy)
end function Iy2

function Iy3(l1,l2,l3)
  integer l1,l2,l3
  double precision Iy3,yl1,yl2,yl3,dy

  yl1=g_rtabsrt(l1,2)
  yl2=g_rtabsrt(l2,2)
  yl3=g_rtabsrt(l3,2)
  dy=g_ymax-g_ymin

  Iy3=  Ix3i(g_ymax,ybar(l1),ybar(l2),ybar(l3),yl1,yl2,yl3,dy) &
       -Ix3i(g_ymin,ybar(l1),ybar(l2),ybar(l3),yl1,yl2,yl3,dy)
end function Iy3

! ######################################################################

function Iz0()
  double precision Iz0
  Iz0=g_zmax-g_zmin
end function Iz0

function Iz1(l)
  integer l
  double precision Iz1,zl,dz

  zl=g_rtabsrt(l,3)
  dz=g_zmax-g_zmin

  Iz1=g_zmax+zbar(l)/dz*(g_zmax-zl)**2/2 - (g_zmin+zbar(l)/dz*(g_zmin-zl)**2/2)
end function Iz1

function Iz2(l1,l2)
  integer l1,l2
  double precision Iz2,zl1,zl2,dz

  zl1=g_rtabsrt(l1,3)
  zl2=g_rtabsrt(l2,3)
  dz=g_zmax-g_zmin

  Iz2=  Ix2i(g_zmax,zbar(l1),zbar(l2),zl1,zl2,dz) &
       -Ix2i(g_zmin,zbar(l1),zbar(l2),zl1,zl2,dz)
end function Iz2

function Iz3(l1,l2,l3)
  integer l1,l2,l3
  double precision Iz3,zl1,zl2,zl3,dz

  zl1=g_rtabsrt(l1,3)
  zl2=g_rtabsrt(l2,3)
  zl3=g_rtabsrt(l3,3)
  dz=g_zmax-g_zmin

  Iz3=  Ix3i(g_zmax,zbar(l1),zbar(l2),zbar(l3),zl1,zl2,zl3,dz) &
       -Ix3i(g_zmin,zbar(l1),zbar(l2),zbar(l3),zl1,zl2,zl3,dz)
end function Iz3

! ######################################################################

function Ix2i(x,xbarl1,xbarl2,xl1,xl2,dx)
  double precision Ix2i,x,xl1,xl2,dx
  integer xbarl1,xbarl2
  Ix2i=x + xbarl1*(x-xl1)**2/(2*dx) + xbarl2*(x-xl2)**2/(2*dx) + xbarl1*xbarl2/dx**2*(x**3/3-x**2/2*(xl1+xl2)+xl1*xl2*x)
end function Ix2i

function Ix3i(x,xbarl1,xbarl2,xbarl3,xl1,xl2,xl3,dx)
  double precision Ix3i,x,xl1,xl2,xl3,dx
  integer xbarl1,xbarl2,xbarl3

  Ix3i=x+xbarl1*(x-xl1)**2/(2*dx)+xbarl2*(x-xl2)**2/(2*dx)+xbarl3*(x-xl3)**2/(2*dx) &
       +xbarl1*xbarl2/dx**2*(x**3/3-x**2/2*(xl1+xl2)+xl1*xl2*x) &
       +xbarl1*xbarl3/dx**2*(x**3/3-x**2/2*(xl1+xl3)+xl1*xl3*x) &
       +xbarl2*xbarl3/dx**2*(x**3/3-x**2/2*(xl2+xl3)+xl2*xl3*x) &
       +xbarl1*xbarl2*xbarl3/dx**3*(x**4/4-x**3/3*(xl1+xl2+xl3)+x**2/2*(xl1*xl2+xl1*xl3+xl2*xl3)-xl1*xl2*xl3*x)
end function Ix3i

! ######################################################################

subroutine runjactest(iunit,Nres0,vec0,delta_vec)

! ----------------------------------------------------------------------
! For testing, calculate the Jacobian by perturbation
! perturb one dof at a time based on vec0 and recalculate Q, R (must be cuboids)
! Hence calculate new negative residual and compare with Nres0
! to form Jacobian directly. delta_vec is perturbation
! Write results out to iunit
! ----------------------------------------------------------------------

  use common, only : Qval,iQ,jQ,lenQ,ndof,irowst,irowed,imax,jmax,kmax,nvar,RR,iunk
  use indexq
  use heap

  integer iunit 
  double precision Nres0(:),vec0(:),delta_vec(:)

  double precision Jval(nvar*ndof*27), Qval0(nvar*ndof*27)
  double precision Nres(ndof),vec(ndof),val
  integer i,j,k,ii,ialpha,ibeta,ip,i1,j1,k1,ii1
  logical found

  Qval0=Qval ! stores integrated Jacobian
  vec0=vec
  Nres0=Nres
  Jval=0     ! use iQ, jQ indexing

  do i=0,imax                  ! loop over columns
  do j=0,jmax
  do k=0,kmax
  do ii=1,nvar
     ialpha=indQ(i,j,k,ii)
                 
     vec=vec0
     vec(ialpha)=vec0(ialpha)+delta_vec(ii)
     
     call getQcuboid(.false.)
     
     RR=0                     ! getR does not zero
     call getRCuboid(.false.)
     
     if (.not.sorted(iQ,lenQ)) call heapsort(iQ,lenQ)    ! heapsort calls swap
     call archive
     
     do ip=1,lenQ               ! replace equations corresponding to fixed DOFs
        if (iunk(iQ(ip)) == 0) then
           if (iQ(ip) == jQ(ip)) then
              Qval(ip)=1
           else
              Qval(ip)=0
           endif
        endif
     enddo
     
     do ip=1,ndof                ! to be equal to u=const
        if (iunk(ip) == 0) then
           RR(ip)=vec0(ip)
        endif
     enddo
     
     ! new residual
     do i1=1,ndof
        Nres(i1)=RR(i1)
        do ip=irowst(i1),irowed(i1)
           Nres(i1)=Nres(i1)-Qval(ip)*vec(jQ(ip))
        enddo
     enddo
     
     do i1=0,imax                   ! loop over rows for each column
     do j1=0,jmax
     do k1=0,kmax
     do ii1=1,nvar
        ibeta=indQ(i1,j1,k1,ii1)
        
        if (Nres(ibeta) == Nres0(ibeta)) cycle
        val=-(Nres(ibeta)-Nres0(ibeta))/delta_vec(ii1)   
        
        !          J(ibeta,ialpha)=-(Nres(ibeta)-Nres0(ibeta))/delta_vec(ii1)
        found=.false.
        do ip=irowst(ibeta),irowed(ibeta)
           if (jQ(ip) == ialpha) then
              found=.true.
              Jval(ip)=Jval(ip)+val      ! Jacobian stored in Q
              exit
           endif
        enddo
        
        if (.not.found) then
           print *,'Newton jactest, element not found'
           print *,'ibeta=',ibeta,':',i1,j1,k1
           print *,'ialpha=',ialpha,':',i,j,k
           stop
        endif
        
     enddo
     enddo
     enddo
     enddo

  enddo
  enddo
  enddo
  enddo

  write (iunit,*) 'Archive:'
  do i=1,ndof
     write (iunit,*) i,irowst(i),irowed(i)
  enddo
  
  write (iunit,*) 'J: (int and perturbation)'
  write (iunit,'(2a10,2a13)') 'iQ','jQ','Jint','Jpert'
  
  do i=1,min(lenQ,5000)
     write (iunit,'(2i10,2e13.5)') iQ(i),jQ(i),Qval0(i),Jval(i)
  enddo
  
end subroutine runjactest

end module cuboid
