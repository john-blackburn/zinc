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

  real a(0:*)
  integer count,end

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

  real a(0:*)
  integer count,start
  
  start=count/2-1

  do while (start.ge.0)
     call siftdown(a,start,count-1)
     start=start-1
  enddo

end subroutine heapify

! ######################################################################

subroutine siftdown(a,start,end)

  real a(0:*)
  integer start,end,root,child
  
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

end module heap
