! Part of Zinc FE package. Author: John Blackburn

module interpolate

contains

subroutine interp(ii)

! ----------------------------------------------------------------------
!     Set initial values in u array
!     rnode gives node positions
!     ireg, iregfix to tell if node is fixed (for a particular variable)
!     and regfix its fixed value
!
!     The function reads specifications from a private file, NOT unit 1
! ----------------------------------------------------------------------

  use common
  use indexQ
  use matrices, only : iregfix
  implicit none

  integer ii,i,j,k,ireg1,i2,j2,k2,ireg2
  double precision x,y,z,top,bot,xr,yr,zr,fr,d2

! ----------------------------------------------------------------------
!     set variable values
! ----------------------------------------------------------------------

  do i=0,imax
     do j=0,jmax
        do k=0,kmax
           x=rnode(i,j,k,1)
           y=rnode(i,j,k,2)
           z=rnode(i,j,k,3)
           ireg1=ireg(i,j,k)
        
           if (iregfix(ireg1,ii).eq.0) then
                 
! ----------------------------------------------------------------------
!     simple interpolation, but could be slow
!     Interpolate variable on the list and set others to zero
! ----------------------------------------------------------------------

              top=0
              bot=0
              do i2=0,imax
                 do j2=0,jmax
                    do k2=0,kmax
                       ireg2=ireg(i2,j2,k2)
                       if (iregfix(ireg2,ii).eq.1) then
                          xr=rnode(i2,j2,k2,1)
                          yr=rnode(i2,j2,k2,2)
                          zr=rnode(i2,j2,k2,3)
                          fr=vec(indQ(i2,j2,k2,ii))
                          
                          d2=(x-xr)**2+(y-yr)**2+(z-zr)**2
                          if (d2.eq.0) then
                             write (*,*) 'interp: coincident node'
                             call exit(1)
                          else
                             top=top+fr/d2
                             bot=bot+1/d2
                          endif
                       endif
                    enddo
                 enddo
              enddo
              
              vec(indQ(i,j,k,ii))=top/bot
           endif
        enddo
     enddo
  enddo                     ! over i,j,k
  
end subroutine interp

end module interpolate
