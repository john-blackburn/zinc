! Part of Zinc FE package. Author: John Blackburn

module heap

! heapsort
! heapify
! siftdown
! swap

contains

subroutine heapsort(a,count)

! ----------------------------------------------------------------------
!     subroutine to sort integer array a of length count into increasing
!     order. User must supply 
!     subroutine swap(a,i,j) 
!     integer a(0:*),i,j
!     to swap elements i and j over (this sub can also swap auxiliary arrays
!     passed to it by common block, if needed)
! ----------------------------------------------------------------------

  integer a(0:*),count,end

  print *,'Sorting Q matrix'

  call heapify(a,count)
  
  end=count-1

  do while (end.gt.0)

     call swap(a,end,0)
     end=end-1
     call siftdown(a,0,end)
         
  enddo

end subroutine heapsort

! ######################################################################

subroutine heapify(a,count)

  integer a(0:*),count,start
  
  start=count/2-1

  do while (start.ge.0)
     call siftdown(a,start,count-1)
     start=start-1
  enddo

end subroutine heapify

! ######################################################################

subroutine siftdown(a,start,end)

  integer a(0:*),start,end,root,child
  
  root=start
  
  do while (root*2+1.le.end)
     
     child=root*2+1
     
     if (child.lt.end.and.a(child).lt.a(child+1)) child=child+1
         
     if (a(root).lt.a(child)) then
        call swap(a,root,child)
        root=child
     else
        return
     endif

  enddo

end subroutine siftdown

! ######################################################################

subroutine swap(ia,i,j)

! ----------------------------------------------------------------------
! This routine is called by subroutine heapsort
! swap elements i and j of array ia
! can also swap other arrays in same way, note i,j are from 0:count
! Note the definition of ia must be ia(0:*)
! Uses jQ, Qval from common blocks. These are indexed from 1:*
! ----------------------------------------------------------------------

  use common
  implicit none
  
  integer ia(0:*),itemp,i,j
  integer i1,j1
  double precision temp
  
  itemp=ia(i)
  ia(i)=ia(j)
  ia(j)=itemp
  
  i1=i+1
  j1=j+1
  
  itemp=jQ(i1)
  jQ(i1)=jQ(j1)
  jQ(j1)=itemp
  
  temp=Qval(i1)
  Qval(i1)=Qval(j1)
  Qval(j1)=temp
  
end subroutine swap

! ######################################################################

function sorted(ia,count)

  logical sorted
  integer ia(:),count,i

  do i=1,count-1
     if (ia(i)>ia(i+1)) then
        sorted=.false.
        return
     endif
  enddo

  sorted=.true.
end function sorted

end module heap
