! Part of Zinc FE package. Author: John Blackburn

module update

contains

! ######################################################################

  subroutine update_transient(wrt_conv)

    use common
    implicit none
    
    logical wrt_conv
    integer i,ip,jQ1
    double precision tot,diff,lhs,lhsI,resid,rhs,rhsI,vec2

    diff=0
    do i=1,ndof
       if (iunk(i).eq.1) then
          tot=0
          do ip=irowst(i),irowed(i)
             tot=tot-Qval(ip)*vec(jQ(ip))
          enddo

          vec2=vec(i)+tstep*(tot+RR(i))/kfac(i)
          diff=diff+abs(vec2-vec(i))
          vec(i)=vec2
       endif
    enddo
    diff=diff/ndof

    if (wrt_conv) then

       resid=0
       lhs=0
       rhs=0
        
       do i=1,ndof
          if (iunk(i) == 1) then

             lhsI=0
             rhsI=RR(i)
             
             ! Q(i,j)=Q(iQ(ip),jQ(ip))=Qval(ip)
             do ip=irowst(i),irowed(i)
                jQ1=jQ(ip)
                if (iunk(jQ1) == 1) then
                   lhsI=lhsI+Qval(ip)*vec(jQ1)
                else
                   rhsI=rhsI-Qval(ip)*vec(jQ1)
                endif
             enddo
             
             resid=resid+(lhsI-rhsI)**2
             lhs=lhs+lhsI**2
             rhs=rhs+rhsI**2
          endif
       enddo

       resid=sqrt(resid)
       lhs=sqrt(lhs)
       rhs=sqrt(rhs)
       
       write (*,'(i7,4e14.6)') istep,resid/lhs,lhs,rhs,diff
       write (1,'(i7,4e14.6)') istep,resid/lhs,lhs,rhs,diff
       
!!$       if ((key_sim.eq.0.or.key_sim.eq.1).and.resid/lhs.lt.tol) then
!!$          print '(a,e13.5,a,e13.5)','Convergence condition reached: ',resid/lhs,' < ',tol
!!$          write (1,'(a,e13.5,a,e13.5)') 'Convergence condition reached: ',resid/lhs,' < ',tol
!!$          exit
!!$       endif

    endif

  end subroutine update_transient

end module update
