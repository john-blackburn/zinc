module integrate

  implicit none

  private
  public lockface

contains

! ######################################################################

  subroutine lockface(itabsrt,itabsrtface,icoord,coord,facelst)

! ----------------------------------------------------------------------
! Given sorted notes itabsrt (element) and itabsrtface
! find icoord, coord indicating which face of itabsrt is the same
! as itabsrtface.
! Note that the elface below are in order so that parameters eg (xi,eta)
! on face increase from 1st to 3rd point. Thus compatible
! with face integration
! ----------------------------------------------------------------------

    use iofile, only : prterr

    integer itabsrt(8,3),itabsrtface(4,3),facelst(4),icoord
    double precision coord
    
    integer elface(6,4),icoords(6),itriplet(3),jtriplet(3)
    integer i,j,ielface,starti,startj,tot,ind
    double precision coords(6)
    
    elface(1,:)=[2,3,7,6]; icoords(1)=1; coords(1)=+1.0
    elface(2,:)=[1,4,8,5]; icoords(2)=1; coords(2)=-1.0
    
    elface(3,:)=[4,3,7,8]; icoords(3)=2; coords(3)=+1.0
    elface(4,:)=[1,2,6,5]; icoords(4)=2; coords(4)=-1.0
    
    elface(5,:)=[5,6,7,8]; icoords(5)=3; coords(5)=+1.0
    elface(6,:)=[1,2,3,4]; icoords(6)=3; coords(6)=-1.0
    
    do ielface=1,6
       
! ----------------------------------------------------------------------
! for ielface, find first equivalent point
! ----------------------------------------------------------------------
     
       starti=-1
       
       outer: do i=1,4
          itriplet=itabsrtface(i,:)
          
          do j=1,4
             jtriplet=itabsrt(elface(ielface,j),:)
             if (all(itriplet.eq.jtriplet)) then
                starti=i
                startj=j
                exit outer
             endif
          enddo
       enddo outer
       
! ----------------------------------------------------------------------
! loop over a face and check if all triplets are the same
! between this face of the element and the face
! ----------------------------------------------------------------------
     
       if (starti.ne.-1) then
          
          tot=0
          do ind=0,3
             i=starti+ind
             j=startj+ind
             
             if (i.gt.4) i=i-4
             if (i.lt.1) i=i+4
             
             if (j.gt.4) j=j-4
             if (j.lt.1) j=j+4
             
             itriplet=itabsrtface(i,:)
             jtriplet=itabsrt(elface(ielface,j),:)
             if (all(itriplet.eq.jtriplet)) tot=tot+1
          enddo
          
          if (tot.eq.4) then
             icoord=icoords(ielface)
             coord = coords(ielface)
             facelst=elface(ielface,:)
             return
          endif
          
! ----------------------------------------------------------------------
! Try it again in the opposite direction
! ----------------------------------------------------------------------

          tot=0
          do ind=0,3
             i=starti+ind
             j=startj-ind
             
             if (i.gt.4) i=i-4
             if (i.lt.1) i=i+4
             
             if (j.gt.4) j=j-4
             if (j.lt.1) j=j+4
             
             itriplet=itabsrtface(i,:)
             jtriplet=itabsrt(elface(ielface,j),:)
             if (all(itriplet.eq.jtriplet)) tot=tot+1
          enddo
          
          if (tot.eq.4) then
             icoord=icoords(ielface)
             coord = coords(ielface)
             facelst=elface(ielface,:)
             return
          endif
       endif
    enddo
    
    call prterr('Did not find equivalent face: surface integrals must "point" into the simulation volume')
  
  end subroutine lockface

! ######################################################################

end module integrate
