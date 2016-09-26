! Part of Zinc FE package. Author: John Blackburn

module indexQ

! function indQ
! subroutine archive
! subroutine symmetrise
! subroutine closeup
! subroutine expand

contains

! ######################################################################

function indQ(i,j,k,ii)
  use common
  implicit none
  integer i,j,k,ii,indQ

  if (NDIM==3) then
     indQ=(k*(imax+1)*(jmax+1)+j*(imax+1)+i)*nvar+ii
  else if (NDIM==2) then
     indQ=(j*(imax+1)+i)*nvar+ii
  else
     indQ=i*nvar+ii
  endif

end function indQ

! ######################################################################

subroutine archive

! ----------------------------------------------------------------------
! find start and end of Q rows (now sorted in row order)
! Must be able to cope when there are missing rows due to some rows
! having only zero values. This could happen if C, a are zero for
! all elements connected to a node
! Rows with no non-zeros have irowst=1, irowed=0
! ----------------------------------------------------------------------

  use common
  implicit none

  integer i,ip,irow
  
  write (*,*) 'Archiving Q matrix'
  
  nnod=(imax+1)*(jmax+1)*(kmax+1)
  ndof=nnod*nvar
  
  do i=1,ndof
     irowst(i)=1
     irowed(i)=0
  enddo
  
  irow=0
  
  do ip=1,lenQ
     if (iQ(ip).ne.irow) then
        if (irow.gt.0) irowed(irow)=ip-1
        irow=iQ(ip)
        irowst(irow)=ip
     endif
  enddo
  
  irowed(irow)=lenQ
  
end subroutine archive

! ######################################################################

subroutine symmetrise
      
! ----------------------------------------------------------------------
! Symmetrize Qij -> (Qij+Qji)/2, write matrix
! This checks pattern of non-zeros is symmetrical and numerically
! sets values, but cannot alter pattern if it is wrong. In that case
! it exits with failure message. Loop must be j=1..N, not j=i+1..N
! for complete check on pattern
! ----------------------------------------------------------------------

  use common
  implicit none

  integer iQ1,jQ1,ip,ip2
  double precision av
  
  write (*,*) 'Symmetrising Q matrix...'
  
  nnod=(imax+1)*(jmax+1)*(kmax+1)
  ndof=nnod*nvar
  
  do iQ1=1,ndof
     do ip=irowst(iQ1),irowed(iQ1)
        
        jQ1=jQ(ip)
        do ip2=irowst(jQ1),irowed(jQ1)
           if (jQ(ip2).eq.iQ1) then
              av=0.5*(Qval(ip)+Qval(ip2))
              Qval(ip)=av
              Qval(ip2)=av
              goto 3
           endif
        enddo
        
        write (*,*) 'Error, could not find Q',jQ1,',',iQ1
        call exit(1)
        
3       continue
     enddo
  enddo
      
  write (*,*) 'Done.'

end subroutine symmetrise

! ######################################################################

subroutine closeup

! ----------------------------------------------------------------------
! Prepare vecred (ndofred) from vec
! Also close up RR
! Also close up irowst/ed
! NB does not destroy vec
! ----------------------------------------------------------------------

  use common
  implicit none

  integer i,j,ip,k,ind,ialpha,ii

  vecred=vec
  ndofred=ndof
  iunkred=iunk

! ----------------------------------------------------------------------
! Prepare lookup table so we can recover the unknowns from vec
! In same order as indQ
! ----------------------------------------------------------------------

  ind=0
  do k=0,kmax
     do j=0,jmax
        do i=0,imax
           do ii=1,nvar
              ialpha=indQ(i,j,k,ii)
              if (iunkred(ialpha)==1) then
                 ind=ind+1
                 veclookup(ialpha)=ind
              else
                 veclookup(ialpha)=-1    ! specify fixed
              endif
           enddo
        enddo
     enddo
  enddo
     
  do ip=1,lenQ
     iQ(ip)=veclookup(iQ(ip))
     jQ(ip)=veclookup(jQ(ip))
  enddo

! ----------------------------------------------------------------------
! Now, close up vec, RR, also iunk, irowst, irowed.
! ----------------------------------------------------------------------

  i=1
  do while (i<=ndofred)
!     print *,'i,ndofred=',i,ndofred
     if (iunkred(i)==0) then
        
        do j=i,ndofred-1
           iunkred(j)=iunkred(j+1)
           vecred(j)=vecred(j+1)
           RR(j)=RR(j+1)
           irowst(j)=irowst(j+1)
           irowed(j)=irowed(j+1)
        enddo
        ndofred=ndofred-1
     else
        i=i+1
     endif
  enddo

  do i=1,ndofred
     if (iunkred(i)==0) then
        print *,"Error during closeup"
        stop
     endif
  enddo

  do i=1,ndofred
     do ip=irowst(i),irowed(i)
        if (iQ(ip)/=i) then
           print *,'Error during closeup bad iQ:',i,ip,iQ(ip)
           stop
        endif
     enddo
  enddo

end subroutine closeup

! ######################################################################

   subroutine expand

! ----------------------------------------------------------------------
! Prepare vec from vecred. Does not destroy vecred
! ----------------------------------------------------------------------

     use common
     implicit none

     integer i
     do i=1,ndof
        if (veclookup(i) /= -1) vec(i)=vecred(veclookup(i))
     enddo
   end subroutine expand

end module indexQ
