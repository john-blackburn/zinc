! Part of Zinc FE package. Author: John Blackburn

module solvers

  implicit none
  private
  public solve,freesolvers

  integer leniw,lenw,itol
  integer, parameter :: nsave=10,iunit=1,isym=0
  double precision, allocatable :: rwork(:),Ax(:)
  integer, allocatable :: iwork(:),Ap(:),Ai(:)
  logical :: init=.true.

contains
  
  subroutine sor(ndofred,RR,vecred,lenQ,iQ,jQ,Qval,tol,itmax,iter,residDiag,ierr,iunit, &
       irowst,irowed,omega,nstride)  ! <-- these are SOR specific
    
    integer ndofred,lenQ,iQ(:),jQ(:),itmax,iter,irowst(:),irowed(:),iunit,ierr,nstride
    double precision RR(:),vecred(:),Qval(:),tol,residDiag,omega

    integer i,ip,jQ1
    double precision Qval1,fac,tot,lhs,rhs,lhsI,rhsI
    double precision lhsDiag,rhsDiag,lhsIDiag,rhsIDiag,resid,pivot

    if (iunit>0) write (iunit,'(a5,3a13)') 'iter','resid','lhs','rhs'
    write (*,'(a5,3a13)') 'iter','resid','lhs','rhs'

!    print *,'SOR:'
!    print *,ndofred
    
    do iter=0,itmax

       if (mod(iter,nstride)==0) then
          resid=0; residDiag=0
          lhs=0; lhsDiag=0
          rhs=0; rhsDiag=0

!          print '(a5,9a13)','i','lhsI','rhsI','pivot','lhsIDiag','rhsIDiag','resid','residDiag','lhs','rhs'
          do i=1,ndofred
             lhsI=0
             rhsI=RR(i)

             pivot=0
             do ip=irowst(i),irowed(i)
                jQ1=jQ(ip)
                if (jQ1==-1) then
                   print *,'Bad vecred:',ip
                   stop
                endif

                lhsI=lhsI+Qval(ip)*vecred(jQ1)
                if (jQ1==i) pivot=Qval(ip)
             enddo

             lhsIDiag=lhsI/pivot
             rhsIDiag=rhsI/pivot
             
             resid=resid+(lhsI-rhsI)**2
             lhs=lhs+lhsI**2
             rhs=rhs+rhsI**2

             residDiag=residDiag+(lhsIDiag-rhsIDiag)**2
             lhsDiag=lhsDiag+lhsIDiag**2
             rhsDiag=rhsDiag+rhsIDiag**2

!             print '(i5,9e13.5)',i,lhsI,rhsI,pivot,lhsIDiag,rhsIDiag,resid,residDiag,lhs,rhs
          enddo
!          stop

          resid=sqrt(resid)/sqrt(lhs); residDiag=sqrt(residDiag)/sqrt(lhsDiag)
          lhs=sqrt(lhs); lhsDiag=sqrt(lhsDiag)
          rhs=sqrt(rhs); rhsDiag=sqrt(rhsDiag)

          if (iunit>0) write (iunit,'(i5,6e13.5)') iter,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag
          write (*,'(i5,6e13.5)') iter,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag

          if (residDiag<tol) then
             ierr=0   ! converged
             return
          endif
       endif
       
       do i=1,ndofred
          tot=0
          do ip=irowst(i),irowed(i)
             jQ1=jQ(ip)
             Qval1=Qval(ip)
             
             if (jQ1==i) then
                fac=Qval1
             else
                tot=tot+Qval1*vecred(jQ1)
             endif
          enddo
          
          vecred(i)=(1-omega)*vecred(i)+omega*(RR(i)-tot)/fac
       enddo

    enddo

    ierr=1    ! reached itmax
    
  end subroutine sor

! ######################################################################

  subroutine ssor(ndofred,RR,vecred,lenQ,iQ,jQ,Qval,tol,itmax,iter,residDiag,ierr,iunit, &
       irowst,irowed,omega,nstride)
    
    integer ndofred,lenQ,iQ(:),jQ(:),itmax,iter,irowst(:),irowed(:),iunit,ierr,nstride
    double precision RR(:),vecred(:),Qval(:),tol,residDiag,omega

    integer i,ip,jQ1
    double precision Qval1,fac,tot,lhs,rhs,lhsI,rhsI
    double precision lhsDiag,rhsDiag,lhsIDiag,rhsIDiag,resid,pivot


    if (iunit>0) write (iunit,'(a5,3a13)') 'iter','resid','lhs','rhs'
    write (*,'(a5,3a13)') 'iter','resid','lhs','rhs'

    do iter=0,itmax

       if (mod(iter,nstride)==0) then
          resid=0; residDiag=0
          lhs=0; lhsDiag=0
          rhs=0; rhsDiag=0
          do i=1,ndofred
             lhsI=0
             rhsI=RR(i)
             
             pivot=0
             do ip=irowst(i),irowed(i)
                jQ1=jQ(ip)
                lhsI=lhsI+Qval(ip)*vecred(jQ1)
                if (jQ1==i) pivot=Qval(ip)
             enddo

             lhsIDiag=lhsI/pivot
             rhsIDiag=rhsI/pivot
             
             resid=resid+(lhsI-rhsI)**2
             lhs=lhs+lhsI**2
             rhs=rhs+rhsI**2

             residDiag=residDiag+(lhsIDiag-rhsIDiag)**2
             lhsDiag=lhsDiag+lhsIDiag**2
             rhsDiag=rhsDiag+rhsIDiag**2
          enddo

          resid=sqrt(resid)/sqrt(lhs); residDiag=sqrt(residDiag)/sqrt(lhsDiag)
          lhs=sqrt(lhs); lhsDiag=sqrt(lhsDiag)
          rhs=sqrt(rhs); rhsDiag=sqrt(rhsDiag)

          if (iunit>0) write (iunit,'(i5,6e13.5)') iter,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag
          write (*,'(i5,6e13.5)') iter,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag

          if (residDiag<tol) then
             ierr=0    ! converged
             return
          endif
       endif

       do i=1,ndofred
          tot=0
          do ip=irowst(i),irowed(i)
             jQ1=jQ(ip)
             Qval1=Qval(ip)
             
             if (jQ1==i) then
                fac=Qval1
             else
                tot=tot+Qval1*vecred(jQ1)
             endif
          enddo
          
          vecred(i)=(1-omega)*vecred(i)+omega*(RR(i)-tot)/fac
       enddo
       
       do i=ndofred,1,-1
          tot=0
          do ip=irowst(i),irowed(i)
             jQ1=jQ(ip)
             Qval1=Qval(ip)
             
             if (jQ1==i) then
                fac=Qval1
             else
                tot=tot+Qval1*vecred(jQ1)
             endif
          enddo
          
          vecred(i)=(1-omega)*vecred(i)+omega*(RR(i)-tot)/fac
       enddo

    enddo

    ierr=1        ! reached itmax

  end subroutine ssor

! ######################################################################

   subroutine solve(ndofx,vecx,rhsx,iter,err,ierrslap,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag)

! ----------------------------------------------------------------------
! Solve Q * vecx = rhsx. On input vecx is initial guess, on output holds solution
! The matrix is always Q [stored as iQ,jQ,Qval] but can actually hold Jacobian
! rhsx might be RR or negative residual vector
! ndofx (ndof or ndofred) and lenQ will always be the same if this is called repeatedly for NL
! NOTE: calling the slap routines converts iQ, jQ, Qval to column format
! ie it DESTROYS iQ,jQ,Qval
! ----------------------------------------------------------------------

     use common
     use umfpack_mod
     use iofile

     integer ndofx,iter,ierrslap
     double precision vecx(:),rhsx(:),err

     double precision lhs,rhs,lhsI,rhsI,resid
     double precision lhsDiag,rhsDiag,lhsIDiag,rhsIDiag,residDiag,pivot
     double precision lhsvec(ndofx)   ! auto arrays
     integer i,j,ip,ist,ied,iQx,ret,jQ1
     integer Numeric,Symbolic                ! C pointers

! ----------------------------------------------------------------------
! First call, allocate workspace arrays for SLAP and UMFPACK
! ----------------------------------------------------------------------

!     print *,'Initsolvers:'
!     print *,ndof,ndofred,ndofx

     if (init) then
        print *,'call initsolvers'
        init=.false.
        call initsolvers(ndofx)
     endif

! ----------------------------------------------------------------------
! Run sor, ssor (this file) or SLAP routines
! ----------------------------------------------------------------------

     if (key_sim == 1) then
        print *,'Starting SOR solver'
        call sor(ndofx,rhsx,vecx,lenQ,iQ,jQ,Qval,tol,itmax,iter,err,ierrslap,iunit, &
             irowst,irowed,omega,nstride)
     else if (key_sim == 0) then
        print *,'Starting Symmetric SOR (SSOR) solver'
        call ssor(ndofx,rhsx,vecx,lenQ,iQ,jQ,Qval,tol,itmax,iter,err,ierrslap,iunit, &
             irowst,irowed,omega,nstride)

     else if (key_sim == 3) then
        print *,'Starting GMRES (Diag scaled)'      ! d s d gmr
        call dsdgmr(ndofx,rhsx,vecx,lenQ,iQ,jQ,Qval,isym,nsave,itol,tol,&
             itmax,iter,err,ierrslap,iunit,rwork,lenw,iwork,leniw)
        
     else if (key_sim == 4) then
        print *,'Starting GMRES (Incomplete LU)'    ! d s lu gm
        call dslugm(ndofx,rhsx,vecx,lenQ,iQ,jQ,Qval,isym,nsave,itol,tol,&
             itmax,iter,err,ierrslap,iunit,rwork,lenw,iwork,leniw)
        
     else if (key_sim == 5) then                    ! d s d bcg
        print *,'Starting Biconjugate Gradient (Diag scaled)'
        call dsdbcg(ndofx,rhsx,vecx,lenQ,iQ,jQ,Qval,isym,itol,tol, &
             itmax,iter,err,ierrslap,iunit,rwork,lenw,iwork,leniw)
        
     else if (key_sim == 6) then                    ! d s lu bc
        print *,'Starting Biconjugate Gradient (Incomplete LU)'
        call dslubc(ndofx,rhsx,vecx,lenQ,iQ,jQ,Qval,isym,itol,tol, &
             itmax,iter,err,ierrslap,iunit,rwork,lenw,iwork,leniw)

     else if (key_sim == 7) then                    ! d s d cg
        print *,'Starting Conjugate Gradient method (Diag scaled)'
        print *,'WARNING: this should only be used for symmetric, +ve definite systems'
        
        call dsdcg (ndofx,rhsx,vecx,lenQ,iQ,jQ,Qval,isym,itol,tol, &
             itmax,iter,err,ierrslap,iunit,rwork,lenw,iwork,leniw)
        
     else if (key_sim == 8) then                    ! d s ic cg
        print *,'Starting Conjugate Gradient method (Incomplete Cholesky factorisation)'
        print *,'WARNING: this should only be used for symmetric, +ve definite systems'
        
        call dsiccg(ndofx,rhsx,vecx,lenQ,iQ,jQ,Qval,isym,itol,tol, &
             itmax,iter,err,ierrslap,iunit,rwork,lenw,iwork,leniw)

! ----------------------------------------------------------------------
! Run UMFPACK
! ----------------------------------------------------------------------

     else if (key_sim == 9) then
        print *,'Starting UMFPACK (direct solver)'
!        print *,ndofx,lenQ,init

        do ip=1,lenQ         ! need base 0
           iQ(ip)=iQ(ip)-1
           jQ(ip)=jQ(ip)-1
        enddo

!        allocate (Ap(ndofx+1),Ai(lenQ),Ax(lenQ),stat=allocerr)      ! In column format

        print *,'Set to column format'
        ret=umfpack_di_triplet_to_col(ndofx,ndofx,lenQ,iQ,jQ,Qval,Ap,Ai,Ax,NULL)
        if (ret /= UMFPACK_OK) call prterr("umfpack_di_triplet_to_col failed")

        print *,'Symbolic calculation'
        ret=umfpack_di_symbolic(ndofx, ndofx, Ap, Ai, Ax, Symbolic, NULL, NULL) 
        if (ret /= UMFPACK_OK) call prterr("umfpack_di_symbolic failed")

        print *,'Numeric calculation'
        ret=umfpack_di_numeric(Ap, Ai, Ax, Symbolic, Numeric, NULL, NULL) 
        if (ret /= UMFPACK_OK) call prterr("umfpack_di_numeric failed")

        print *,'Symbolic=',Symbolic
        print *,'Numeric=',Numeric
        
        call umfpack_di_free_symbolic(Symbolic) 

        print *,'Solve'
        ret=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, vecx, rhsx, Numeric, NULL, NULL) 
        if (ret /= UMFPACK_OK) call prterr("umfpack_di_solve failed")

        call umfpack_di_free_numeric(Numeric)

!        deallocate (Ap,Ai,Ax)

        do ip=1,lenQ           ! back to base 1
           iQ(ip)=iQ(ip)+1
           jQ(ip)=jQ(ip)+1
        enddo

        iter=-1        ! irrelevant for Direct solvers
        err=-1
        ierrslap=-1

     else
        call prterr("solve: unknown key_sim option")
     endif

! ----------------------------------------------------------------------
! Calculate BOTH "regular" and reduced/diag residual. In case of SLAP routine, 
! iQ, jQ, Qval have been changed to SLAP column format
! ----------------------------------------------------------------------

     if (key_sim==1.or.key_sim==0.or.key_sim==9) then
        resid=0; residDiag=0
        lhs=0; lhsDiag=0
        rhs=0; rhsDiag=0
        do i=1,ndofx
           lhsI=0
           rhsI=rhsx(i)

           pivot=0
           do ip=irowst(i),irowed(i)
              jQ1=jQ(ip)
              if (jQ1>ndofx) then
                 print *,'error',ip,iQ(ip),i,jQ1
                 stop
              endif

              if (jQ1==i) pivot=Qval(ip)
              lhsI=lhsI+Qval(ip)*vecx(jQ1)
           enddo

           lhsIDiag=lhsI/pivot
           rhsIDiag=rhsI/pivot
           
           resid=resid+(lhsI-rhsI)**2
           lhs=lhs+lhsI**2
           rhs=rhs+rhsI**2

           residDiag=residDiag+(lhsIDiag-rhsIDiag)**2
           lhsDiag=lhsDiag+lhsIDiag**2
           rhsDiag=rhsDiag+rhsIDiag**2
        enddo
        
        resid=sqrt(resid)/sqrt(lhs); residDiag=sqrt(residDiag)/sqrt(lhsDiag)
        lhs=sqrt(lhs); lhsDiag=sqrt(lhsDiag)
        rhs=sqrt(rhs); rhsDiag=sqrt(rhsDiag)
     else
        lhsvec=0
        do j=1,ndofx   ! columns
           ist=jQ(j)
           ied=jQ(j+1)-1
           do ip=ist,ied
              iQx=iQ(ip)
              lhsvec(iQx)=lhsvec(iQx)+Qval(ip)*vecx(j)
           enddo
        enddo

        resid=0
        lhs=0
        rhs=0

        residDiag=0
        lhsDiag=0
        rhsDiag=0

        do i=1,ndofx
           lhsI=lhsvec(i)
           rhsI=rhsx(i)

           lhsIDiag=lhsvec(i)/Qval(jQ(i))
           rhsIDiag=rhsx(i)/Qval(jQ(i))

           resid=resid+(lhsI-rhsI)**2
           lhs=lhs+lhsI**2
           rhs=rhs+rhsI**2

           residDiag=residDiag+(lhsIDiag-rhsIDiag)**2
           lhsDiag=lhsDiag+lhsIDiag**2
           rhsDiag=rhsDiag+rhsIDiag**2
        enddo
        resid=sqrt(resid)/sqrt(lhs)
        lhs=sqrt(lhs)
        rhs=sqrt(rhs)

        residDiag=sqrt(residDiag)/sqrt(lhsDiag)
        lhsDiag=sqrt(lhsDiag)
        rhsDiag=sqrt(rhsDiag)
     endif

     print *,'Direct residual:'
     write (*,'(3e13.5)') resid,lhs,rhs
     print *,'Normalised residual:'
     write (*,'(3e13.5)') residDiag,lhsDiag,rhsDiag
   end subroutine solve

! ######################################################################

   subroutine initsolvers(ndofx)

! ----------------------------------------------------------------------
! Set parameters for SLAP and allocate workspace arrays
! Note that lenQ is guaranteed not to vary simply as a result
! of non-linearity so we can allocate once only
! (alternative would be to make leniw etc "big enough")
! itol=2 means use |D^-1 (Qu-R)|_2/|D^-1 R|_2 for tol
! This is good for removeFixed=NO
! ----------------------------------------------------------------------

     use common, only : lenQ,key_sim
     use iofile

     integer ndofx,allocerr
     
     if (key_sim == 3) then     ! dsdgmr
        lenw=2+ndofx*(nsave+7)+nsave*(nsave+3)
        leniw=31
        itol=0

     else if (key_sim == 4) then   ! dslugm
        lenw=1+ndofx*(nsave+7)+nsave*(nsave+3)+lenQ*2
        leniw=lenQ*2+4*ndofx+32
        itol=0

     else if (key_sim == 5) then   ! dsdbcg
        lenw=8*ndofx+1
        leniw=11
        itol=2

     else if (key_sim == 6) then   ! dslubc
        lenw=lenQ*2+8*ndofx
        leniw=lenQ*2+4*ndofx+12
        itol=2

     else if (key_sim == 7) then   ! dsdcg
        lenw=5*ndofx+1
        leniw=11
        itol=2

     else if (key_sim == 8) then   ! dsiccg
        lenw=lenQ+5*ndofx
        leniw=lenQ+ndofx+11
        itol=2
     endif

     if (key_sim>=3.and.key_sim<=8) then
        allocate (rwork(lenw),iwork(leniw),stat=allocerr)
        if (allocerr.ne.0) call prterr('Failed to allocate SLAP work variables')
     endif

     if (key_sim==9) then
!        print *,'allocating Ap,Ai,Ax'
        allocate (Ap(ndofx+1),Ai(lenQ),Ax(lenQ),stat=allocerr)      ! In column format
        if (allocerr /= 0) call prterr("callslap: Failed to allocate Ap, Ai, Ax")
     endif

   end subroutine initsolvers

! ######################################################################

   subroutine freesolvers

     if (allocated(rwork)) then
        deallocate (rwork,iwork)
     endif

     if (allocated(Ax)) then
        deallocate (Ax,Ap,Ai)
     endif
   end subroutine freesolvers

end module solvers
