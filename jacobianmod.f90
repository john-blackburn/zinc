! Part of Zinc FE package. Author: John Blackburn

module jacobianmod

implicit none

! New functions for non-linear calculation using Newton-Raphson method
! Supports only C, F matrices and only DLL non-linear functions
! (not single line expressions)
! In this case we need to also link dCfun_du, dCfun_ddu
! In this case we need to also link dffun_du, dffun_ddu
!
! Simplified getQR here does NOT handle surface integrals but DOES
! (optionally) integrate material properties over element.
! Does not use helper functions like getelinfo, getelint as these
! are not suitable to new functionality and also obscures what we
! are doing.
! (can add ALL sparse matrix entries including zeros if this causes problems)
!
! Should be able to test proof of principle for cathodetoken
! and then fullfuel which both only use C, f, fixed nodes, no SURFACE

! getJacobianC
! getJacobianF
! getJacobian2C (by elements)
! getJacobian2F (by elements)
! getsol
! addtoQ

! note flags jac_by_element and centre_eval

contains

! ######################################################################

subroutine getJacobian

  use common, only : jac_by_element

  if (jac_by_element) then
     call getJacobian2C
     call getJacobian2F
  else
     call getJacobianC
     call getJacobianF
  endif

end subroutine getJacobian

! ######################################################################

subroutine getJacobianC

! ----------------------------------------------------------------------
! Update Q matrix to equal the Jacobian. iQ, jQ. Q must have already been
! set and heapsorted and archived (irowst,irowed). We assume J has the
! same sparsity structure as Q, if not we are in trouble!
! This routine loops over node pairs (and triplets) and then loops
! over elements connected to these to discover the the Jacobian
! This routine add C-related terms
! ----------------------------------------------------------------------

  use common
  use geom
  use shape
  use matrices
  use indexQ

  integer ia,ja,ka,ig,jg,kg,ib,jb,kb,ie,je,ke,iep,ind,la,lg,lb
  integer inode,jnode,knode,icomb(2,3)
  logical igin,jgin,kgin,ibin,jbin,kbin,iein,jein,kein,foundg,foundb

  integer itab(8,4,3),itabsrt(8,3)
  double precision rtabsrt(8,3),rtab(8,3),ul(8,nvar),dCdu(nvar),dCddu(nvar,3)
  double precision elint(nvar,nvar,nvar),elint1(3,3),elint2(3,3,3)
  double precision jac(3,3),invjac(3,3),jacdet,r(3)
  double precision term2(nvar,nvar),dNa(3),dNg(3),dNb(3)

  integer ixi,ieta,imu,i1,j1,k1,irege
  double precision xi,eta,mu,dur(nvar,3),ur(nvar),tot,trailer
  integer ii,jj,kk,ll,mm,nn,icol,irow,l,ip
  character(EXPRLEN) token

  integer iusel(0:imax-1,0:jmax-1,0:kmax-1)

  integer, parameter :: ies(8)=[1, 1, 1, 1,-1,-1,-1,-1]
  integer, parameter :: jes(8)=[1, 1,-1,-1, 1, 1,-1,-1]
  integer, parameter :: kes(8)=[1,-1, 1,-1, 1,-1, 1,-1]

  print *,'Preparing Jacobian (getJacobianC)'

  iusel=0

  do ia=0,imax

     print *,'i-plane: ',ia,' / ',imax

  do ja=0,jmax
  do ka=0,kmax          ! alpha node loop

     do ig=ia-1,ia+1
     do jg=ja-1,ja+1
     do kg=ka-1,ka+1    ! gamma node loop

        igin = ig>=0.and.ig<=imax
        jgin = jg>=0.and.jg<=jmax
        kgin = kg>=0.and.kg<=kmax

        if (.not.(igin.and.jgin.and.kgin)) cycle ! gamma node must exist

        term2=0           ! will be term2 (for each ii,nn) by end of beta loop

        do ib=ia-1,ia+1
        do jb=ja-1,ja+1
        do kb=ka-1,ka+1

           ibin = ib>=0.and.ib<=imax
           jbin = jb>=0.and.jb<=jmax
           kbin = kb>=0.and.kb<=kmax
           
           if (.not.(ibin.and.jbin.and.kbin)) cycle  ! node must exist

! ----------------------------------------------------------------------
! Loop over elements connected to alpha, gamma and beta (if any exist)
! Note alpha=gamma etc are allowed
! ----------------------------------------------------------------------
        
           elint=0         ! will be sum of integrals when iep loop completes
              
           do iep=1,8
              ie=ies(iep)
              je=jes(iep)
              ke=kes(iep)

              iein = ia+ie>=0.and.ia+ie<=imax
              jein = ja+je>=0.and.ja+je<=jmax
              kein = ka+ke>=0.and.ka+ke<=kmax

              if (.not.(iein.and.jein.and.kein)) cycle    ! element must exist

! ----------------------------------------------------------------------
! Find itab, rtab for element. Use itab to check both gamma and beta
! nodes are attached to candidate element else cycle element loop
! If all is well use sort to determine itabsrt, rtabsrt
! ----------------------------------------------------------------------

              icomb(1,1)=ia
              icomb(1,2)=ja
              icomb(1,3)=ka

              icomb(2,1)=ia+ie
              icomb(2,2)=ja+je
              icomb(2,3)=ka+ke

              call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)

              foundg=.false.
              foundb=.false.

              do ind=1,8
                 if (itab(ind,1,1)==ig.and.itab(ind,1,2)==jg.and.itab(ind,1,3)==kg) foundg=.true.
                 if (itab(ind,1,1)==ib.and.itab(ind,1,2)==jb.and.itab(ind,1,3)==kb) foundb=.true.
              enddo

              if (.not.(foundg.and.foundb)) cycle

              call sort(itab,rtab,itabsrt,rtabsrt)

! ----------------------------------------------------------------------
! Find local numbers for alpha, gamma, beta: la, lg, lb
! ----------------------------------------------------------------------

              la=0; lg=0; lb=0
              do ip=1,8
                 inode=itabsrt(ip,1)
                 jnode=itabsrt(ip,2)
                 knode=itabsrt(ip,3)

                 if (inode==ia.and.jnode==ja.and.knode==ka) la=ip
                 if (inode==ig.and.jnode==jg.and.knode==kg) lg=ip
                 if (inode==ib.and.jnode==jb.and.knode==kb) lb=ip
              enddo

              if (la==0.or.lg==0.or.lb==0) then
                 print *,'Error: getJacobian: failed to find la,lg or lb'
                 stop
              endif

! ----------------------------------------------------------------------
! Find the element's region number and local solution at nodes
! surrounding element
! Update iusel which states how many times element used
! ----------------------------------------------------------------------

              i1=ia+(ie-1)/2
              j1=ja+(je-1)/2
              k1=ka+(ke-1)/2
              irege=iregup(i1,j1,k1)

              do l=1,8
                 inode=itabsrt(l,1)
                 jnode=itabsrt(l,2)
                 knode=itabsrt(l,3)
                 
                 do ii=1,nvar
                    ul(l,ii)=vec(indQ(inode,jnode,knode,ii))
                 enddo
              enddo
              
              iusel(i1,j1,k1)=iusel(i1,j1,k1)+1

! ----------------------------------------------------------------------
! Integrate
! ----------------------------------------------------------------------

              elint1=0            ! For centre_eval
              elint2=0            ! For centre_eval

              do ixi=1,ng_xi
              do ieta=1,ng_eta
              do imu=1,ng_mu

                 xi=gauss(ng_xi,ixi)
                 eta=gauss(ng_eta,ieta)
                 mu=gauss(ng_mu,imu)

                 call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet)

                 call getdN(la,xi,eta,mu,invjac,dNa)   ! dN_alpha/dr
                 call getdN(lg,xi,eta,mu,invjac,dNg)
                 call getdN(lb,xi,eta,mu,invjac,dNb)

                 trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)*wt(ng_mu,imu)

! ----------------------------------------------------------------------
! For centre eval, calculate shape function integrals
! ----------------------------------------------------------------------

                 if (centre_eval) then
                    do jj=1,3
                       do ll=1,3
                          elint1(jj,ll)=elint1(jj,ll)+dNa(jj)*dNb(ll)*Nl(lg,xi,eta,mu)*trailer
                          do mm=1,3
                             elint2(jj,ll,mm)=elint2(jj,ll,mm)+dNa(jj)*dNb(ll)*dNg(mm)*trailer
                          enddo
                       enddo
                    enddo

                 else

! ----------------------------------------------------------------------
! else Calculate local solution, accumulate elint directly
! ----------------------------------------------------------------------
                    
                    call getsol(ul,invjac,rtabsrt, xi,eta,mu, ur,dur,r)
                    
                    do ii=1,nvar
                    do kk=1,nvar
                    do nn=1,nvar

                       do jj=1,3
                       do ll=1,3

                          if (iCC(irege,ii,jj,kk,ll)==3) then
                             
                             token=sCC(irege,ii,jj,kk,ll)

                             call dCfun_du (token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                                  ireg,iregup,rnode,vec,imax,jmax,kmax,dCdu)

                             call dCfun_ddu(token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                                  ireg,iregup,rnode,vec,imax,jmax,kmax,dCddu)
                          
!                             dCdu=  dCfun_du(token,ur,dur)   ! variation with u
!                             dCddu=dCfun_ddu(token,ur,dur)   ! variation with du/dr
                             
                             tot=dCdu(nn)*Nl(lg,xi,eta,mu)
                             do mm=1,3
                                tot=tot+dCddu(nn,mm)*dNg(mm)
                             enddo
                             
                             elint(ii,kk,nn)=elint(ii,kk,nn)+dNa(jj)*dNb(ll)*tot*trailer
                          endif
                          
                       enddo
                       enddo

                    enddo
                    enddo
                    enddo
                 endif
              
              enddo    ! integration loops
              enddo
              enddo

! ----------------------------------------------------------------------
! For centre eval, form elint, which is a sum over integrals for
! each element. Also sum over jj,ll for each ii,kk,nn
! ----------------------------------------------------------------------

              if (centre_eval) then

                 call jacobian(zero,zero,zero, rtabsrt,jac,invjac,jacdet)
                 call getsol(ul,invjac, rtabsrt, zero,zero,zero, ur,dur,r)

                 do ii=1,nvar
                 do kk=1,nvar
                 do nn=1,nvar

                    do jj=1,3
                    do ll=1,3

                       if (iCC(irege,ii,jj,kk,ll)==3) then
                          token=sCC(irege,ii,jj,kk,ll)

!                          dCdu =dCfun_du (token,ur,dur)
!                          dCddu=dCfun_ddu(token,ur,dur)

                          call dCfun_du (token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                               ireg,iregup,rnode,vec,imax,jmax,kmax,dCdu)
                          
                          call dCfun_ddu(token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                               ireg,iregup,rnode,vec,imax,jmax,kmax,dCddu)
                    
                          elint(ii,kk,nn)=elint(ii,kk,nn)+dCdu(nn)*elint1(jj,ll) &
                               +dCddu(nn,1)*elint2(jj,ll,1)+dCddu(nn,2)*elint2(jj,ll,2)+dCddu(nn,3)*elint2(jj,ll,3)

                       endif
                    enddo
                    enddo
                 enddo
                 enddo
                 enddo
              endif

           enddo   ! loop over elements

! ----------------------------------------------------------------------
! We have now summed integrals over all elements connected to alpha, beta, gamma
! nodes. Can now prepare contribution to J(i,alpha,n,gamma)
! ----------------------------------------------------------------------
           
           do ii=1,nvar
           do nn=1,nvar
              do kk=1,nvar
                 term2(ii,nn)=term2(ii,nn)+vec(indQ(ib,jb,kb,kk))*elint(ii,kk,nn)
              enddo
           enddo
           enddo

        enddo    ! loop over beta
        enddo
        enddo

! ----------------------------------------------------------------------
! Now add term2 to Q
! ----------------------------------------------------------------------

        do ii=1,nvar
           do nn=1,nvar
              irow=indQ(ia,ja,ka,ii)
              icol=indQ(ig,jg,kg,nn)

!              print *,'irow=',ia,ja,ka,ii
!              print *,'icol=',ig,jg,kg,nn
              call addtoQ(irow,icol,term2(ii,nn))
           enddo
        enddo

     enddo     ! loop over gamma
     enddo
     enddo

  enddo        ! loop over alpha
  enddo
  enddo

  print *,'El usage: ',minval(iusel),' to ',maxval(iusel)

end subroutine getJacobianC

! ######################################################################

subroutine getJacobianF

! ----------------------------------------------------------------------
! As getJacobianC but now add f term
! ----------------------------------------------------------------------

  use common
  use shape
  use geom
  use matrices
  use indexQ

  integer ia,ja,ka,ig,jg,kg,ie,je,ke,iep,ind,la,lg
  integer inode,jnode,knode,irege
  logical igin,jgin,kgin,iein,jein,kein,foundg

  integer icomb(2,3),itab(8,4,3),itabsrt(8,3)
  double precision rtabsrt(8,3),rtab(8,3),ul(8,nvar),dfdu(nvar),dfddu(nvar,3)
  double precision elint(nvar,nvar),elint1,elint2(3),dNa(3),dNg(3)
  double precision jac(3,3),invjac(3,3),jacdet,tot,r(3)

  integer ixi,ieta,imu,ip
  double precision xi,eta,mu,dur(nvar,3),ur(nvar),trailer
  integer ii,mm,nn,icol,irow,l,i1,j1,k1
  character(EXPRLEN) token

  integer iusel(0:imax-1,0:jmax-1,0:kmax-1)

  integer, parameter :: ies(8)=[1, 1, 1, 1,-1,-1,-1,-1]
  integer, parameter :: jes(8)=[1, 1,-1,-1, 1, 1,-1,-1]
  integer, parameter :: kes(8)=[1,-1, 1,-1, 1,-1, 1,-1]

  iusel=0

  print *,'Preparing Jacobian (getJacobianF)'

  do ia=0,imax

     print *,'i-plane: ',ia,' / ',imax

  do ja=0,jmax
  do ka=0,kmax

     do ig=ia-1,ia+1
     do jg=ja-1,ja+1
     do kg=ka-1,ka+1

        igin = ig>=0.and.ig<=imax
        jgin = jg>=0.and.jg<=jmax
        kgin = kg>=0.and.kg<=kmax

        if (.not.(igin.and.jgin.and.kgin)) cycle ! node must exist

        elint=0
        do iep=1,8
           ie=ies(iep)
           je=jes(iep)
           ke=kes(iep)

           iein = ia+ie>=0.and.ia+ie<=imax
           jein = ja+je>=0.and.ja+je<=jmax
           kein = ka+ke>=0.and.ka+ke<=kmax

           if (.not.(iein.and.jein.and.kein)) cycle    ! element must exist

           icomb(1,1)=ia
           icomb(1,2)=ja
           icomb(1,3)=ka
           
           icomb(2,1)=ia+ie
           icomb(2,2)=ja+je
           icomb(2,3)=ka+ke
           
           call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
           
           foundg=.false.
           do ind=1,8
              if (itab(ind,1,1)==ig.and.itab(ind,1,2)==jg.and.itab(ind,1,3)==kg) foundg=.true.
           enddo
           
           if (.not.foundg) cycle

           call sort(itab,rtab,itabsrt,rtabsrt)

! ----------------------------------------------------------------------
! Discover if alpha, beta and gamma nodes are all within the element
! (alpha is bound to be but check anyway)
! ----------------------------------------------------------------------

           la=0; lg=0
           do ip=1,8
              inode=itabsrt(ip,1)
              jnode=itabsrt(ip,2)
              knode=itabsrt(ip,3)
              
              if (inode==ia.and.jnode==ja.and.knode==ka) la=ip
              if (inode==ig.and.jnode==jg.and.knode==kg) lg=ip
           enddo

           if (la==0.or.lg==0) then
              print *,'Error: getJacobian: failed to find la,lg or lb'
              stop
           endif
           
! ----------------------------------------------------------------------
! Find the element's region number and local solution at nodes
! surrounding element
! Update iusel which states how many times element used
! ----------------------------------------------------------------------

           i1=ia+(ie-1)/2
           j1=ja+(je-1)/2
           k1=ka+(ke-1)/2
           irege=iregup(i1,j1,k1)

           do l=1,8
              inode=itabsrt(l,1)
              jnode=itabsrt(l,2)
              knode=itabsrt(l,3)
              
              do ii=1,nvar
                 ul(l,ii)=vec(indQ(inode,jnode,knode,ii))
              enddo
           enddo
           
           iusel(i1,j1,k1)=iusel(i1,j1,k1)+1

! ----------------------------------------------------------------------
! Integrate
! ----------------------------------------------------------------------

           elint1=0
           elint2=0

           do ixi=1,ng_xi
           do ieta=1,ng_eta
           do imu=1,ng_mu

              xi=gauss(ng_xi,ixi)
              eta=gauss(ng_eta,ieta)
              mu=gauss(ng_mu,imu)

              call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet)
              
              call getdN(la,xi,eta,mu,invjac,dNa)   ! dN_alpha/dr
              call getdN(lg,xi,eta,mu,invjac,dNg)
              
              trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)*wt(ng_mu,imu)

              if (centre_eval) then
                 elint1=elint1+Nl(la,xi,eta,mu)*Nl(lg,xi,eta,mu)*trailer
                 do mm=1,3
                    elint2(mm)=elint2(mm)+Nl(la,xi,eta,mu)*dNg(mm)*trailer
                 enddo
              else

                 call getsol(ul,invjac,rtabsrt, xi,eta,mu, ur,dur,r)
                               
                 do ii=1,nvar
                 do nn=1,nvar
                       
                    if (iff(irege,ii)==3) then
                       token=sff(irege,ii)
                       
!                       dfdu=dffun_du(token,ur,dur)
!                       dfddu=dffun_ddu(token,ur,dur)

                       call dffun_du (token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                            ireg,iregup,rnode,vec,imax,jmax,kmax,dfdu)

                       call dffun_ddu(token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                            ireg,iregup,rnode,vec,imax,jmax,kmax,dfddu)
                       
                       tot=dfdu(nn)*Nl(lg,xi,eta,mu)
                       do mm=1,3
                          tot=tot+dfddu(nn,mm)*dNg(mm)
                       enddo
                       
                       elint(ii,nn)=elint(ii,nn)+Nl(la,xi,eta,mu)*tot*trailer
                    endif
                 enddo
                 enddo

              endif

           enddo   ! end integral
           enddo
           enddo

! ----------------------------------------------------------------------
! In case of centre_eval
! ----------------------------------------------------------------------

           if (centre_eval) then

              call jacobian(zero,zero,zero,rtabsrt,jac,invjac,jacdet)
              call getsol(ul,invjac,rtabsrt, zero,zero,zero, ur,dur,r)

              do ii=1,nvar
              do nn=1,nvar

                 if (iff(irege,ii)==3) then
                    token=sff(irege,ii)

!                    dfdu=dffun_du(token,ur,dur)
!                    dfddu=dffun_ddu(token,ur,dur)

                    call dffun_du (token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                         ireg,iregup,rnode,vec,imax,jmax,kmax,dfdu)

                    call dffun_ddu(token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                         ireg,iregup,rnode,vec,imax,jmax,kmax,dfddu)
                    
                    elint(ii,nn)=elint(ii,nn)+dfdu(nn)*elint1 &
                         +dfddu(nn,1)*elint2(1) &
                         +dfddu(nn,2)*elint2(2) &
                         +dfddu(nn,3)*elint2(3)
                 endif

              enddo
              enddo
           endif

        enddo      ! end element loop

        do ii=1,nvar
           do nn=1,nvar
              irow=indQ(ia,ja,ka,ii)
              icol=indQ(ig,jg,kg,nn)
              call addtoQ(irow,icol,-elint(ii,nn))
           enddo
        enddo

     enddo           ! gamma
     enddo
     enddo

  enddo              ! alpha
  enddo
  enddo

  print *,'El usage: ',minval(iusel),' to ',maxval(iusel)

end subroutine getJacobianF

! ######################################################################

subroutine getJacobian2C

! ----------------------------------------------------------------------
! Update Q matrix to equal the Jacobian. iQ, jQ. Q must have already been
! set and heapsorted and archived (irowst,irowed). We assume J has the
! same sparsity structure as Q, if not we are in trouble!
! This routine loops over elements then considers pairs or triplets
! at edges of each element to fill Jacobian
! ----------------------------------------------------------------------

  use common
  use shape
  use geom
  use matrices
  use indexQ

  integer i,j,k,in,jn,kn,irege,inode,jnode,knode
  integer icomb(2,3),itab(8,4,3),itabsrt(8,3)
  double precision rtabsrt(8,3),rtab(8,3),ul(8,nvar),dCdu(nvar),dCddu(nvar,3)
  double precision elint(nvar,nvar,nvar),elint1(3,3),elint2(3,3,3),jac(3,3),invjac(3,3),jacdet
  integer la,lg,lb,ia,ja,ka,ig,jg,kg,ib,jb,kb
  integer ixi,ieta,imu
  double precision xi,eta,mu,dur(nvar,3),ur(nvar),tot,trailer,dNa(3),dNb(3),dNg(3),r(3)
  integer ii,jj,kk,ll,mm,nn,irow,icol,l
  character(EXPRLEN) token
  
! ----------------------------------------------------------------------
! Loop over elements
! ----------------------------------------------------------------------

  print *,'Preparing Jacobian (getJacobian2C)'

  do i=0,imax-1

     print *,'i-plane: ',i,' / ',imax-1

  do j=0,jmax-1
  do k=0,kmax-1

     in=i+1
     jn=j+1
     kn=k+1

     irege=iregup(i,j,k)

     icomb(1,1)=i
     icomb(1,2)=j
     icomb(1,3)=k

     icomb(2,1)=in
     icomb(2,2)=jn
     icomb(2,3)=kn

     call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
     call sort(itab,rtab,itabsrt,rtabsrt)

! ----------------------------------------------------------------------
! Get local solution at nodes surrounding element
! ----------------------------------------------------------------------

     do l=1,8
        inode=itabsrt(l,1)
        jnode=itabsrt(l,2)
        knode=itabsrt(l,3)
        
        do ii=1,nvar
           ul(l,ii)=vec(indQ(inode,jnode,knode,ii))
        enddo
     enddo

! ----------------------------------------------------------------------
! Consider every possible triplet of nodes surrounding elements (C matrix)
! ----------------------------------------------------------------------

     do la=1,8
     do lg=1,8
     do lb=1,8
        ia=itabsrt(la,1)
        ja=itabsrt(la,2)
        ka=itabsrt(la,3)
        
        ig=itabsrt(lg,1)
        jg=itabsrt(lg,2)
        kg=itabsrt(lg,3)
        
        ib=itabsrt(lb,1)
        jb=itabsrt(lb,2)
        kb=itabsrt(lb,3)

! ----------------------------------------------------------------------
! Calculate integrals over the element
! ----------------------------------------------------------------------

        elint=0
        elint1=0
        elint2=0
        
        do ixi=1,ng_xi
        do ieta=1,ng_eta
        do imu=1,ng_mu
           xi=gauss(ng_xi,ixi)
           eta=gauss(ng_eta,ieta)
           mu=gauss(ng_mu,imu)

           call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet)

           call getdN(la,xi,eta,mu,invjac,dNa)   ! dN_alpha/dr
           call getdN(lg,xi,eta,mu,invjac,dNg)
           call getdN(lb,xi,eta,mu,invjac,dNb)

           trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)*wt(ng_mu,imu)

! ----------------------------------------------------------------------
! For each ii,nn,kk combination, accumulate integral for C
! (sum over jj,ll spatial directions)
! All of this is for the alpha, gamma node (beta summed automatically)
! ----------------------------------------------------------------------

           if (centre_eval) then
              do jj=1,3
              do ll=1,3
                 elint1(jj,ll)=elint1(jj,ll)+dNa(jj)*dNb(ll)*Nl(lg,xi,eta,mu)*trailer
                 do mm=1,3
                    elint2(jj,ll,mm)=elint2(jj,ll,mm)+dNa(jj)*dNb(ll)*dNg(mm)*trailer
                 enddo
              enddo
              enddo
           else

              call getsol(ul,invjac,rtabsrt, xi,eta,mu, ur,dur,r)

              do ii=1,nvar
              do nn=1,nvar
              do kk=1,nvar

                 do jj=1,3
                 do ll=1,3

                    ! Allow only token/DLL non-linear for now

                    if (iCC(irege,ii,jj,kk,ll)==3) then

                       token=sCC(irege,ii,jj,kk,ll)

!                    dCdu=  dCfun_du(token,ur,dur)   ! variation with u
!                    dCddu=dCfun_ddu(token,ur,dur)   ! variation with du/dr

                       call dCfun_du (token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                            ireg,iregup,rnode,vec,imax,jmax,kmax,dCdu)

                       call dCfun_ddu(token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                            ireg,iregup,rnode,vec,imax,jmax,kmax,dCddu)
                    
                       tot=dCdu(nn)*Nl(lg,xi,eta,mu)
                       do mm=1,3
                          tot=tot+dCddu(nn,mm)*dNg(mm)
                       enddo
                       
                       elint(ii,kk,nn)=elint(ii,kk,nn)+dNa(jj)*dNb(ll)*tot*trailer
                    endif

                 enddo
                 enddo
              enddo
              enddo
              enddo
           
           endif

        enddo            ! element integral (Gaussian quad)
        enddo
        enddo

! ----------------------------------------------------------------------
! For centre_eval, form elint
! ----------------------------------------------------------------------

        if (centre_eval) then

           call jacobian(zero,zero,zero, rtabsrt,jac,invjac,jacdet)
           call getsol(ul,invjac, rtabsrt, zero,zero,zero, ur,dur,r)

!           print *,'Centre eval!'

           do ii=1,nvar
           do kk=1,nvar
           do nn=1,nvar

              do jj=1,3
              do ll=1,3

              if (iCC(irege,ii,jj,kk,ll)==3) then

                 token=sCC(irege,ii,jj,kk,ll)

!                 dCdu=dCfun_du(token,ur,dur)
!                 dCddu=dCfun_ddu(token,ur,dur)

                 call dCfun_du(token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                      ireg,iregup,rnode,vec,imax,jmax,kmax,dCdu)

                 call dCfun_ddu(token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                      ireg,iregup,rnode,vec,imax,jmax,kmax,dCddu)

!                 print *,jj,ll,elint1(jj,ll),elint2(jj,ll,:)

                 elint(ii,kk,nn)=elint(ii,kk,nn)+dCdu(nn)*elint1(jj,ll) &
                      +dCddu(nn,1)*elint2(jj,ll,1)+dCddu(nn,2)*elint2(jj,ll,2)+dCddu(nn,3)*elint2(jj,ll,3)
              endif

              enddo
              enddo

!              print *,ii,kk,nn,elint(ii,kk,nn)

           enddo
           enddo
           enddo

!           stop

        endif

! ----------------------------------------------------------------------
! Add contributions to the [alpha,i] and [gamma,n] dofs
! There is an implied sum over [beta,k]
! ----------------------------------------------------------------------

        do ii=1,nvar
        do nn=1,nvar

           tot=0
           do kk=1,nvar
              tot=tot+vec(indQ(ib,jb,kb,kk))*elint(ii,kk,nn)
           enddo

           irow=indQ(ia,ja,ka,ii)
           icol=indQ(ig,jg,kg,nn)
           
           call addtoQ(irow,icol,tot)

        enddo
        enddo

     enddo     ! loop over node triplets
     enddo
     enddo

  enddo        ! Element loop
  enddo
  enddo

end subroutine getJacobian2C

! ######################################################################

subroutine getJacobian2F

  use common
  use shape
  use geom
  use matrices
  use indexQ

  integer i,j,k,in,jn,kn,irege,inode,jnode,knode
  integer icomb(2,3),itab(8,4,3),itabsrt(8,3)
  double precision rtabsrt(8,3),ul(8,nvar),dfdu(nvar),dfddu(nvar,3),rtab(8,3)
  double precision elint(nvar,nvar),elint1,elint2(3),jac(3,3),invjac(3,3),jacdet
  integer la,lg,ia,ja,ka,ig,jg,kg,icol,irow,l
  integer ixi,ieta,imu
  double precision xi,eta,mu,dur(nvar,3),ur(nvar),tot,trailer,dNa(3),dNg(3),r(3)
  integer ii,mm,nn
  character(EXPRLEN) token

! ----------------------------------------------------------------------
! Loop over elements
! ----------------------------------------------------------------------

  print *,'Preparing Jacobian (getJacobian2F)'

  do i=0,imax-1

     print *,'i-plane: ',i,' / ',imax-1

  do j=0,jmax-1
  do k=0,kmax-1

     in=i+1
     jn=j+1
     kn=k+1

     irege=iregup(i,j,k)

     icomb(1,1)=i
     icomb(1,2)=j
     icomb(1,3)=k

     icomb(2,1)=in
     icomb(2,2)=jn
     icomb(2,3)=kn

     call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
     call sort(itab,rtab,itabsrt,rtabsrt)

! ----------------------------------------------------------------------
! Get local solution at nodes surrounding element
! ----------------------------------------------------------------------

     do l=1,8
        inode=itabsrt(l,1)
        jnode=itabsrt(l,2)
        knode=itabsrt(l,3)
        
        do ii=1,nvar
           ul(l,ii)=vec(indQ(inode,jnode,knode,ii))
        enddo
     enddo

! ----------------------------------------------------------------------
! Consider every possible pair of nodes surrounding elements (C matrix)
! ----------------------------------------------------------------------

     do la=1,8
     do lg=1,8
        ia=itabsrt(la,1)
        ja=itabsrt(la,2)
        ka=itabsrt(la,3)
        
        ig=itabsrt(lg,1)
        jg=itabsrt(lg,2)
        kg=itabsrt(lg,3)
        
! ----------------------------------------------------------------------
! Calculate integrals over the element
! ----------------------------------------------------------------------

        elint=0
        elint1=0
        elint2=0

        do ixi=1,ng_xi
        do ieta=1,ng_eta
        do imu=1,ng_mu
           xi=gauss(ng_xi,ixi)
           eta=gauss(ng_eta,ieta)
           mu=gauss(ng_mu,imu)

           call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet)

           call getdN(la,xi,eta,mu,invjac,dNa)   ! dN_alpha/dr
           call getdN(lg,xi,eta,mu,invjac,dNg)

           trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)*wt(ng_mu,imu)

           if (centre_eval) then
              elint1=elint1+Nl(la,xi,eta,mu)*Nl(lg,xi,eta,mu)*trailer
              do mm=1,3
                 elint2(mm)=elint2(mm)+Nl(la,xi,eta,mu)*dNg(mm)*trailer
              enddo
           else
  
! ----------------------------------------------------------------------
! For each ii, nn combination, accumulate integral for f
! ----------------------------------------------------------------------

           call getsol(ul,invjac,rtabsrt, xi,eta,mu, ur,dur,r)

           do ii=1,nvar
           do nn=1,nvar
                    
              if (iff(irege,ii)==3) then
                 token=sff(irege,ii)
                 
!                 dfdu=dffun_du(token,ur,dur)
!                 dfddu=dffun_ddu(token,ur,dur)

                 call dffun_du (token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                      ireg,iregup,rnode,vec,imax,jmax,kmax,dfdu)
                 
                 call dffun_ddu(token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                      ireg,iregup,rnode,vec,imax,jmax,kmax,dfddu)
                 
                 tot=dfdu(nn)*Nl(lg,xi,eta,mu)
                 do mm=1,3
                    tot=tot+dfddu(nn,mm)*dNg(mm)
                 enddo
                    
                 elint(ii,nn)=elint(ii,nn)+Nl(la,xi,eta,mu)*tot*trailer
              endif
           enddo
           enddo

           endif

        enddo   ! element integral
        enddo
        enddo

! ----------------------------------------------------------------------
! In case of centre_eval, form integral here
! ----------------------------------------------------------------------

        if (centre_eval) then

           call jacobian(zero,zero,zero, rtabsrt,jac,invjac,jacdet)
           call getsol(ul,invjac,rtabsrt, zero,zero,zero, ur,dur,r)

           do ii=1,nvar
              do nn=1,nvar

                 if (iff(irege,ii)==3) then
                    token=sff(irege,ii)
!                    dfdu=dffun_du(token,ur,dur)
!                    dfddu=dffun_ddu(token,ur,dur)

                    call dffun_du (token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                         ireg,iregup,rnode,vec,imax,jmax,kmax,dfdu)
                    
                    call dffun_ddu(token,r(1),r(2),r(3),ur,dur,nvar,istep, &
                         ireg,iregup,rnode,vec,imax,jmax,kmax,dfddu)

                    elint(ii,nn)=elint(ii,nn)+dfdu(nn)*elint1 &
                         +dfddu(nn,1)*elint2(1)+dfddu(nn,2)*elint2(2)+dfddu(nn,3)*elint2(3)
                 endif
              enddo
           enddo

        endif

! ----------------------------------------------------------------------
! Accumulate onto Q
! ----------------------------------------------------------------------

        do ii=1,nvar
        do nn=1,nvar

           irow=indQ(ia,ja,ka,ii)
           icol=indQ(ig,jg,kg,nn)

           call addtoQ(irow,icol,-elint(ii,nn))
        enddo
        enddo

     enddo         ! gamma
     enddo         ! alpha

  enddo            ! element loop
  enddo
  enddo

end subroutine getJacobian2F

! ######################################################################

subroutine getsol(ul,invjac,rtabsrt, xi,eta,mu,ur,dur,r)

! ----------------------------------------------------------------------
! Given corner solution ul and xi,eta,mu and invjac(xi,eta,mu) calculate
! solution at (xi,eta,mu) for element. Also return r=(x,y,z) coord of point
! ----------------------------------------------------------------------

  use common, only : nvar
  use shape
  use geom

  double precision ul(:,:),invjac(3,3),xi,eta,mu,ur(:),dur(:,:),r(:)
  double precision Nl1,dN(3),rtabsrt(:,:)
  integer l,ii,jj

  ur=0
  do l=1,8
     Nl1=Nl(l,xi,eta,mu)
     do ii=1,nvar
        ur(ii)=ur(ii)+ul(l,ii)*Nl1
     enddo
  enddo
  
  dur=0
  do l=1,8
     call getdN(l,xi,eta,mu,invjac,dN)
     do ii=1,nvar
        do jj=1,3
           dur(ii,jj)=dur(ii,jj)+ul(l,ii)*dN(jj)
        enddo
     enddo
  enddo

  r=0
  do l=1,8
     r=r+rtabsrt(l,:)*Nl(l,xi,eta,mu)
  enddo

end subroutine getsol

! ######################################################################

subroutine addtoQ(irow,icol,val)

  use common, only : irowst,irowed,Qval,iQ,jQ,iunk

  integer irow,icol,ip
  double precision val
  logical found

  if (iunk(irow) == 0) return   ! no contribution to jacobian in fixed node

  if (val /= 0) then

     if (irowst(irow) == 1.and.irowed(irow) == 0) then
        print *,'getJacobian: Found empty row, fail!'
        stop
     endif

     found=.false.
     do ip=irowst(irow),irowed(irow)
!        print *,ip,irow,icol,iQ(ip),jQ(ip)
        if (jQ(ip)==icol) then
           found=.true.
           Qval(ip)=Qval(ip)+val
           exit
        endif
     enddo

     if (.not.found) then
        print *,'addtoQ, element not found',irow,icol,val

        do ip=irowst(irow),irowed(irow)
           print *,ip,iQ(ip),jQ(ip)
        enddo
        stop
     endif
  endif

end subroutine addtoQ

end module jacobianmod
