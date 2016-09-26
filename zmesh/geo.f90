! ######################################################################

module geo

implicit none

contains

  subroutine shift(yshift,zshift)
    use global

    real yshift,zshift
    integer i,j,k,ic,jmx

    if (fastdraw) then
       do i=1,npoly
          if (preg(i).gt.0) then
             jmx=4
          else
             jmx=1
          endif

          do j=1,jmx
             plist(i,j,2)=plist(i,j,2)+yshift
             plist(i,j,3)=plist(i,j,3)+zshift
          enddo
       enddo
    else

       do i=0,imax
          do j=0,jmax
             do k=0,kmax
                pp(i,j,k,2)=pp(i,j,k,2)+yshift
                pp(i,j,k,3)=pp(i,j,k,3)+zshift
             enddo
          enddo
       enddo

       do i=1,nextra          ! Note: rotates even if extra=.false.
          if (extratype(i).eq.1) then
             jmx=1
          else
             jmx=4
          endif
          
          do j=1,jmx
             rextra(i,j,2)=rextra(i,j,2)+yshift
             rextra(i,j,3)=rextra(i,j,3)+zshift
          enddo
       enddo
       
    endif

    do ic=1,3
       do i=1,4
          do j=1,2
             bb(ic,i,j,2)=bb(ic,i,j,2)+yshift
             bb(ic,i,j,3)=bb(ic,i,j,3)+zshift
          enddo
       enddo
    enddo

    do ic=1,3
       do j=1,2
          ff(ic,j,2)=ff(ic,j,2)+yshift
          ff(ic,j,3)=ff(ic,j,3)+zshift
       enddo
    enddo


  end subroutine shift

! ----------------------------------------------------------------------

  subroutine rotx(ang)
    use global
!    use const, only : WM_WDEBUG
!    use user32

    real rp(3),ang,theta,rmag
    integer i,j,k,ret,ic,jmx

!    write (message,*) 'rotx called'
!    ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

    if (fastdraw) then
       do i=1,npoly
          if (preg(i).gt.0) then
             jmx=4
          else
             jmx=1
          endif

          do j=1,jmx
             rp=plist(i,j,:)
             theta=atan2(rp(2),rp(3))
             rmag=sqrt(rp(2)**2+rp(3)**2)
             
             plist(i,j,3)=rmag*cos(theta+ang)
             plist(i,j,2)=rmag*sin(theta+ang)
          enddo
       enddo
    else

       do i=0,imax
          do j=0,jmax
             do k=0,kmax
                
                rp=pp(i,j,k,:)
                theta=atan2(rp(2),rp(3))
                rmag=sqrt(rp(2)**2+rp(3)**2)
                
                pp(i,j,k,3)=rmag*cos(theta+ang)
                pp(i,j,k,2)=rmag*sin(theta+ang)
             enddo
          enddo
       enddo

       do i=1,nextra
          if (extratype(i).eq.1) then
             jmx=1
          else
             jmx=4
          endif
          
          do j=1,jmx
             rp=rextra(i,j,:)
             theta=atan2(rp(2),rp(3))
             rmag=sqrt(rp(2)**2+rp(3)**2)
             
             rextra(i,j,3)=rmag*cos(theta+ang)
             rextra(i,j,2)=rmag*sin(theta+ang)
          enddo
       enddo
       
    endif

    do ic=1,3
       do i=1,4
          do j=1,2
             rp=bb(ic,i,j,:)
             theta=atan2(rp(2),rp(3))
             rmag=sqrt(rp(2)**2+rp(3)**2)
  
             bb(ic,i,j,3)=rmag*cos(theta+ang)
             bb(ic,i,j,2)=rmag*sin(theta+ang)
          enddo
       enddo
    enddo

    do ic=1,3
       do j=1,2
          rp=ff(ic,j,:)
          theta=atan2(rp(2),rp(3))
          rmag=sqrt(rp(2)**2+rp(3)**2)
          
          ff(ic,j,3)=rmag*cos(theta+ang)
          ff(ic,j,2)=rmag*sin(theta+ang)
       enddo
    enddo

  end subroutine rotx

! ----------------------------------------------------------------------

  subroutine roty(ang)
    use global
!    use const, only : WM_WDEBUG
!    use user32

    real rp(3),ang,theta,rmag
    integer i,j,k,ret,ic,jmx

!    write (message,*) 'roty called'
!    ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

    if (fastdraw) then
       do i=1,npoly
          if (preg(i).gt.0) then
             jmx=4
          else
             jmx=1
          endif

          do j=1,jmx
             rp=plist(i,j,:)
             theta=atan2(rp(3),rp(1))
             rmag=sqrt(rp(3)**2+rp(1)**2)
             
             plist(i,j,1)=rmag*cos(theta+ang)
             plist(i,j,3)=rmag*sin(theta+ang)
          enddo
       enddo
    else

       do i=0,imax
          do j=0,jmax
             do k=0,kmax
                rp=pp(i,j,k,:)
                theta=atan2(rp(3),rp(1))
                rmag=sqrt(rp(3)**2+rp(1)**2)
                
                pp(i,j,k,1)=rmag*cos(theta+ang)
                pp(i,j,k,3)=rmag*sin(theta+ang)
             enddo
          enddo
       enddo

       do i=1,nextra
          if (extratype(i).eq.1) then
             jmx=1
          else
             jmx=4
          endif
          
          do j=1,jmx
             rp=rextra(i,j,:)
             theta=atan2(rp(3),rp(1))
             rmag=sqrt(rp(3)**2+rp(1)**2)
             
             rextra(i,j,1)=rmag*cos(theta+ang)
             rextra(i,j,3)=rmag*sin(theta+ang)
          enddo
       enddo

    endif

    do ic=1,3
       do i=1,4
          do j=1,2
             rp=bb(ic,i,j,:)
             theta=atan2(rp(3),rp(1))
             rmag=sqrt(rp(3)**2+rp(1)**2)
             
             bb(ic,i,j,1)=rmag*cos(theta+ang)
             bb(ic,i,j,3)=rmag*sin(theta+ang)
          enddo
       enddo
    enddo

    do ic=1,3
       do j=1,2
          rp=ff(ic,j,:)
          theta=atan2(rp(3),rp(1))
          rmag=sqrt(rp(3)**2+rp(1)**2)
          
          ff(ic,j,1)=rmag*cos(theta+ang)
          ff(ic,j,3)=rmag*sin(theta+ang)
       enddo
    enddo

  end subroutine roty

! ----------------------------------------------------------------------

  subroutine rotz(ang)
    use global
!    use const, only : WM_WDEBUG
!    use user32

    real rp(3),ang,theta,rmag
    integer i,j,k,ret,ic,jmx

!    write (message,*) 'rotz called'
!    ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

    if (fastdraw) then
       do i=1,npoly
          if (preg(i).gt.0) then
             jmx=4
          else
             jmx=1
          endif

          do j=1,jmx
             rp=plist(i,j,:)
             theta=atan2(rp(2),rp(1))
             rmag=sqrt(rp(2)**2+rp(1)**2)
             
             plist(i,j,1)=rmag*cos(theta+ang)
             plist(i,j,2)=rmag*sin(theta+ang)
          enddo
       enddo
    else
    
       do i=0,imax
          do j=0,jmax
             do k=0,kmax
                rp=pp(i,j,k,:)
                theta=atan2(rp(2),rp(1))
                rmag=sqrt(rp(2)**2+rp(1)**2)
                
                pp(i,j,k,1)=rmag*cos(theta+ang)
                pp(i,j,k,2)=rmag*sin(theta+ang)
             enddo
          enddo
       enddo

       do i=1,nextra
          if (extratype(i).eq.1) then
             jmx=1
          else
             jmx=4
          endif
          
          do j=1,jmx
             rp=rextra(i,j,:)
             theta=atan2(rp(2),rp(1))
             rmag=sqrt(rp(2)**2+rp(1)**2)
             
             rextra(i,j,1)=rmag*cos(theta+ang)
             rextra(i,j,2)=rmag*sin(theta+ang)
          enddo
       enddo

    endif

    do ic=1,3
       do i=1,4
          do j=1,2
             rp=bb(ic,i,j,:)
             theta=atan2(rp(2),rp(1))
             rmag=sqrt(rp(2)**2+rp(1)**2)
             
             bb(ic,i,j,1)=rmag*cos(theta+ang)
             bb(ic,i,j,2)=rmag*sin(theta+ang)
          enddo
       enddo
    enddo

    do ic=1,3
       do j=1,2
          rp=ff(ic,j,:)
          theta=atan2(rp(2),rp(1))
          rmag=sqrt(rp(2)**2+rp(1)**2)
             
          ff(ic,j,1)=rmag*cos(theta+ang)
          ff(ic,j,2)=rmag*sin(theta+ang)
       enddo
    enddo

  end subroutine rotz
end module geo

