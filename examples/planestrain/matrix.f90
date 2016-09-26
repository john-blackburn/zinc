module matrix

! function inverse(a)        :  invert matrix (square)
! function det(a)  recursive :  determinant of matrix (square)
! subroutine wmat(a)         :  write matrix (need not be square)

contains

  function inverse(a)

! ======================================================================
! function to calculate inverse of matrix a
! ======================================================================

    implicit none

    integer,parameter :: long=selected_real_kind(12,300)

    real(long) a(:,:),inverse(size(a,1),size(a,2))

    real(long), dimension(size(a,1),size(a,2)) :: c
    real(long), dimension(size(a,1)-1,size(a,2)-1) :: minor
    
    integer, dimension(size(a,1)-1) :: rows, cols
    integer i,j,k,n,ind
    
! ======================================================================
! check matrix is square
! ======================================================================

    if (size(a,1).ne.size(a,2)) then
       write (*,*) 'inverse: a is not square'
       stop
    endif

! ======================================================================
! Calculate cofactor matrix. Take det of minors * sign factor
! ======================================================================

    n=size(a,1)
    
    do i=1,n
       do j=1,n
          
          ind=0
          do k=1,n
             if (k.ne.i) then
                ind=ind+1
                rows(ind)=k
             endif
          enddo
          
          ind=0
          do k=1,n
             if (k.ne.j) then
                ind=ind+1
                cols(ind)=k
             endif
          enddo
          
          minor=a(rows,cols)
          c(i,j)=det(minor)*(-1)**(i+j)
       enddo
    enddo

! ======================================================================
! inv(a)=c^T/det(a)
! ======================================================================

    inverse=transpose(c)/det(a)

  end function inverse

! ######################################################################

  recursive function det(a) result(d)

! ======================================================================
! function to calculate det of a
! note result clause is needed for direct function recursion
! ======================================================================

    implicit none

    integer,parameter :: long=selected_real_kind(12,300)
    
    real(long) d,a(:,:)
    
    real(long) minor(size(a,1)-1,size(a,2)-1)
    integer n,j,k,ind,cols(size(a,2)-1)

! ======================================================================
! Check a is square. Expand along top row
! Special case for n=1, n=2 (which can terminate the recursion)
! ======================================================================

    if (size(a,1).ne.size(a,2)) then
       write (*,*) 'inverse: a is not square'
       stop
    endif
    
    n=size(a,1)

    if (n.eq.1) then
       d=a(1,1)
    else if (n.eq.2) then
       d=a(1,1)*a(2,2)-a(2,1)*a(1,2)
    else
       d=0
       do j=1,n
          
          ind=0
          do k=1,n
             if (k.ne.j) then
                ind=ind+1
                cols(ind)=k
             endif
          enddo
          
          minor=a(2:,cols)
          
          d=d+a(1,j)*det(minor)*(-1)**(1+j)
       enddo
    endif
    
  end function det

! ######################################################################

  subroutine wmat(iunit,name,a)

    implicit none

    integer,parameter :: long=selected_real_kind(12,300)

    real(long) a(:,:)
    character(len=*) name
    integer iunit
    
    integer n,m,i,j
    
    n=size(a,1)
    m=size(a,2)
    
    write (iunit,*)
    write (iunit,*) name
    
    do i=1,n
       write (iunit,'(100e13.5)') (a(i,j),j=1,m)
    enddo
  end subroutine wmat

end module matrix

