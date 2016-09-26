! Part of Zinc FE package. Author: John Blackburn

! Note, I do not have USE statements in the global scope of modules
! as this leads to a cascade. (not if we use global PRIVATE)
! The only thing in the modules global scope should be private internal
! storage for the contained subroutines (and IMPLICIT)
! uses evalexpr (or evaltoken) in CCval, aaval, ffval

module qpr

! subroutine getQ
! subroutine getQ2
! subroutine getQ3
! subroutine getR
! subroutine getk
! function energy
! subroutine getelint
! function elint2d

  implicit none

contains

subroutine getQ(rmFixed)

! ----------------------------------------------------------------------
! Form Q matrix (iQ,jQ,Qval) and its length lenQ (on common)
! from rnode,iregup,C,aa,ff (on common)
! iunk and vec must be prepared before calling this subroutine
! ----------------------------------------------------------------------

  use common
  use geom
  use indexQ
  use matrices
  use util
  use iofile, only : prterr

  logical rmFixed

  integer icomb(2,NDIM),ies(8),jes(8),kes(8),off1s(4),off2s(4)
  integer itab(2**NDIM,NDIM+1,NDIM),itabsrt(2**NDIM,NDIM),iface(6,4)
  double precision rtab(2**NDIM,NDIM),rtabsrt(2**NDIM,NDIM),elint(NDIM,NDIM)
  logical inin,jnin,knin,iein,jein,kein,wtout,elin
  logical found,literal

  integer iusel(0:imax-1,0:jmax-1,0:kmax-1) ! auto array

  integer indf,i,j,k,in,jn,kn,lenQst,iep,i1,j1,k1
  integer n,m,l,ifnd,ifc,ialpha,ibeta
  integer ie,je,ke,ii,jj,kk,ll,ip
  integer inode,jnode,knode,irege

  double precision tot,elintNN,ur(nvar),dur(nvar,NDIM),rc(NDIM),norm(3)
  double precision centre1(3),centre2(3),facecentre(3)
  integer ifirst,iusel1,iusmn,iusmx,l1,l2,idir,irege1,irege2
  integer off1,off2,iaax,iCCx,iqqx
  save ies,jes,kes,iface
  
! ----------------------------------------------------------------------
! Set iface for debugging and ies etc for matrix assembly
! ----------------------------------------------------------------------

  data ((iface(ifc,i),i=1,4),ifc=1,6)/&  ! debug only
       1,2,4,3,&
       5,6,8,7,&
       2,3,6,7,&
       1,4,5,8,&
       1,2,5,6,&
       4,3,8,7/

  data ies/ 1,-1, 1,-1, 1,-1, 1,-1/     ! For 1D only use first 2 ies
  data jes/ 1, 1,-1,-1, 1, 1,-1,-1/     ! For 2D only use first 4 ies, jes
  data kes/ 1, 1, 1, 1,-1,-1,-1,-1/

  data off1s/1, 1,-1,-1/                ! For 3D surfaces
  data off2s/1,-1, 1,-1/

! ----------------------------------------------------------------------
! zero iusel which records number of times each element was used
! ----------------------------------------------------------------------

  print *,'Calculating Q matrix'

  iusel(:,:,:)=0
  if (rmFixed) RR(:)=0  ! zero here since we are adding to RR in this routine

! ----------------------------------------------------------------------
! Prepare Qmatrix
! loop over nodes and neighbour nodes
! loop over elements ie=(-1,+1) etc
! reserved: ijk,ijkn,ijke,ifnd,indf,lenQ
! ----------------------------------------------------------------------

  indf=0
  lenQ=0

  do k=0,kmax    ! for 2D kmax=0
  write (*,*) 'k-plane:',k,'/',kmax
  do j=0,jmax
  do i=0,imax

  wtout=key_db.eq.2.and.&
       i.ge.idbmin.and.i.le.idbmax.and.&
       j.ge.jdbmin.and.j.le.jdbmax.and.&
       k.ge.kdbmin.and.k.le.kdbmax
  
  do in=i-1,i+1
  do jn=j-1,j+1
  do kn=k-1,k+1

     inin=(in.ge.0.and.in.le.imax)
     jnin=(jn.ge.0.and.jn.le.jmax) ! for 1D j=0, jn=[-1,0,1] so only allow jn=0
     knin=(kn.ge.0.and.kn.le.kmax) ! for 2D k=0, kn=[-1,0,1] so only allow kn=0

     if (inin.and.jnin.and.knin) then
        
        lenQst=lenQ+1
        
        do iep=1,2**NDIM
           ie=ies(iep)
           je=jes(iep)
           ke=kes(iep)
           
           iein=(i+ie.ge.0.and.i+ie.le.imax)
           jein=(j+je.ge.0.and.j+je.le.jmax).or.NDIM==1
           kein=(k+ke.ge.0.and.k+ke.le.kmax).or.NDIM<=2
              
           if (iein.and.jein.and.kein) then

! ----------------------------------------------------------------------
! set "combination" array and generate itab first col
! check if (in,jn,kn) is in itab
! ----------------------------------------------------------------------

              if (wtout) then
                 write (1,*) '====================================='
                 write (1,*) 'i,j,k=',i,j,k
                 write (1,*) 'in,jn,kn=',in,jn,kn
                 write (1,*) 'ie,je,ke=',ie,je,ke
              endif
              
              icomb(1,1)=i
              icomb(2,1)=i+ie

              if (NDIM>=2) then
                 icomb(1,2)=j
                 icomb(2,2)=j+je
              endif

              if (NDIM==3) then
                 icomb(1,3)=k
                 icomb(2,3)=k+ke
              endif

              call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
              
              ifnd=0
              if (NDIM==3) then
                 do ip=1,8
                    if (itab(ip,1,1).eq.in.and.itab(ip,1,2).eq.jn.and.itab(ip,1,3).eq.kn) ifnd=1
                 enddo
              else if (NDIM==2) then
                 do ip=1,4
                    if (itab(ip,1,1)==in.or.itab(ip,1,2)==jn) ifnd=1
                 enddo
              else
                 do ip=1,2
                    if (itab(ip,1,1)==in) ifnd=1
                 enddo
              endif

              if (wtout) then
                 do n=1,8
                    write (1,*) (itab(n,1,l),l=1,3)
                 enddo
              endif
              
! ----------------------------------------------------------------------
! If (in,jn,kn) is in element, scan itab to find neighbour nodes
! Thus, fill itab col 2,3,4
! ----------------------------------------------------------------------

              if (ifnd.eq.1) then

                 if (wtout) then
                    write (1,*) 'itab filled'
                    do n=1,8
                       write (1,*) ((itab(n,m,l),l=1,3),'|',m=1,4)
                    enddo
                 endif

                 call sort(itab,rtab,itabsrt,rtabsrt)

! ----------------------------------------------------------------------
! Write debug info
! ----------------------------------------------------------------------
                 
                 if (wtout) then
                    do n=1,8
                       write (1,*) (itabsrt(n,m),m=1,3)
                    enddo
                    
                    indf=indf+1
                    
                    open (2,file='label'//i2c(indf)//'.gnu',&
                         status='unknown')

                    write (2,'(a)') '# prepared by '//name
                    write (2,*) '# i,j,k=',i,j,k
                    write (2,*) '# in,jn,kn=',in,jn,kn
                    write (2,*) '# ie,je,ke=',ie,je,ke
                    
                    do n=1,8
                       write (2,*) 'set label "',n,'" at ',&
                            rtabsrt(n,1),',',rtabsrt(n,2),',',rtabsrt(n,3)
                    enddo

                    close (2)

                    open (2,file='line'//i2c(indf)//'.out',status='unknown')
                    write (2,'(a)') '# prepared by '//name

                    do ifc=1,6
                       write (2,*) (rtabsrt(iface(ifc,1),l),l=1,3)
                       write (2,*) (rtabsrt(iface(ifc,2),l),l=1,3)
                       write (2,*)
                       write (2,*) (rtabsrt(iface(ifc,3),l),l=1,3)
                       write (2,*) (rtabsrt(iface(ifc,4),l),l=1,3)
                       write (2,*)
                       write (2,*)
                    enddo

                    close (2)
                 endif

! ----------------------------------------------------------------------
! find (i,j,k) and (in,jn,kn) in itabsrt, hence get local numbers
! corresponding to these, l1, l2
! ----------------------------------------------------------------------

                 l1=0; l2=0
                 if (NDIM==3) then
                    do ip=1,8
                       inode=itabsrt(ip,1)
                       jnode=itabsrt(ip,2)
                       knode=itabsrt(ip,3)
                       if (inode.eq.i.and.jnode.eq.j.and.knode.eq.k) l1=ip
                       if (inode.eq.in.and.jnode.eq.jn.and.knode.eq.kn) l2=ip
                    enddo
                 else if (NDIM==2) then
                    do ip=1,4
                       inode=itabsrt(ip,1)
                       jnode=itabsrt(ip,2)
                       if (inode.eq.i.and.jnode.eq.j) l1=ip
                       if (inode.eq.in.and.jnode.eq.jn) l2=ip
                    enddo
                 else
                    do ip=1,2
                       inode=itabsrt(ip,1)
                       if (inode.eq.i) l1=ip
                       if (inode.eq.in) l2=ip
                    enddo
                 endif

                 if (l1==0.or.l2==0) call prterr('Error in getQ, bad l1,l2')

! ----------------------------------------------------------------------
! get region number of element
! if ie = 1, i1=i; if ie=-1, i1=i-1
! ----------------------------------------------------------------------

                 i1=i+(ie-1)/2
                 j1=j+(je-1)/2
                 k1=k+(ke-1)/2            ! for 2D, k=0, ke=1 => k1=0
                 irege=iregup(i1,j1,k1)

                 iusel(i1,j1,k1)=iusel(i1,j1,k1)+1
                 
                 if (wtout) then
                    write (1,*) 'l1,l2,i1,j1,k1,irege=',&
                         l1,l2,i1,j1,k1,irege
                 endif

! ----------------------------------------------------------------------
! Integrate over element using Gauss Quad
! sub jacobian returns jac, invjac and jacdet
! elint(jj,ll)=int_e dN/dxj dN/dxl, elintNN = int_e N N
! where the N's are associated with ijk and in,jn,kn resp.
! ----------------------------------------------------------------------

                 call getelint(l1,l2,rtabsrt,elint,elintNN)

! ----------------------------------------------------------------------
! In non-linear case, we need to get the ur, dur for centre
! of this element and its centroid (xc,yc,zc)
! ----------------------------------------------------------------------

                              rc(1)=sum(rtabsrt(:,1))/2**NDIM
                 if (NDIM>=2) rc(2)=sum(rtabsrt(:,2))/2**NDIM
                 if (NDIM==3) rc(3)=sum(rtabsrt(:,3))/2**NDIM

                 if (regCNL(irege).or.regANL(irege)) then
                    call getelinfo(itabsrt,rtabsrt,ur,dur)
                 endif

! ----------------------------------------------------------------------
! add contributions to matrix elements
! lenQst was set at lenQ+1 inside the ijk, ijkn loops
! We need only search from lenQst to lenQ since we will not access
! previous elements of Q again. Because Q depends on ijk, ijkn.
! Matrix element has two contributions for this finite element,
! One due to C which is called Q in report and one due to a which is called P
! Here Q+P -> Q
! We do not add zero entries. However entries coincidentally zero due to
! Non-linearity are not considered as zero and will be added to Q
! We do not add fixed dofs but instead accumulate onto RR.
! Note that RR will have more content added in getR (due to f)
! ----------------------------------------------------------------------

                 do ii=1,nvar
                    do kk=1,nvar
                       ialpha=indQ(i,j,k,ii)
                       if (rmFixed.and.iunk(ialpha)==0) cycle
                       ibeta=indQ(in,jn,kn,kk)

                       literal=.true.
                       tot=0
                       do jj=1,NDIM
                          do ll=1,NDIM
                             tot=tot+CCval(irege,ii,jj,kk,ll,ur,dur,rc,iCCx)*elint(jj,ll)
                             if (iCCx /= 0) literal=.false.
                          enddo
                       enddo
                       
                       tot=tot+aaval(irege,ii,kk,ur,dur,rc,iaax)*elintNN
                       if (iaax /=0) literal=.false.

                       if (.not.(tot==0.0.and.literal).or.newton=="YES") then

                          if (rmFixed.and.iunk(ibeta)==0) then
                             RR(ialpha)=RR(ialpha)-tot*vec(ibeta)
                          else

                             found=.false.
                             do ip=lenQst,lenQ
                                if (iQ(ip).eq.ialpha.and.jQ(ip).eq.ibeta) then
                                   Qval(ip)=Qval(ip)+tot
                                   found=.true.
                                   exit
                                endif
                             enddo

                             if (.not.found) then
                                lenQ=lenQ+1
                                iQ(lenQ)=ialpha
                                jQ(lenQ)=ibeta
                                Qval(lenQ)=tot
                             endif

                          endif

                       endif
                    enddo
                 enddo
                 
! ----------------------------------------------------------------------
! end of loops
! ----------------------------------------------------------------------

              endif         ! element has both nodes, ifnd=1
           endif            ! element exists
        enddo               ! element search around node ijk (iep)

! ----------------------------------------------------------------------
! Now do surface integration Q term. Get 12 possible faces
! ----------------------------------------------------------------------

        if (ntable.gt.0) then

        if (NDIM/=3) call prterr('surfaces not supported for 1D or 2D')

        do idir=1,3
        do iep=1,4
   
           off1=off1s(iep)
           off2=off2s(iep)
           
           if (idir.eq.1) then
              elin=i+off1.le.imax.and.i+off1.ge.0.and.&
                   j+off2.le.jmax.and.j+off2.ge.0
           else if (idir.eq.2) then
              elin=i+off1.le.imax.and.i+off1.ge.0.and.&
                   k+off2.le.kmax.and.k+off2.ge.0
           else
              elin=j+off1.le.jmax.and.j+off1.ge.0.and.&
                   k+off2.le.kmax.and.k+off2.ge.0
           endif
           
           if (elin) then            ! face must exist
              
              if (idir.eq.1) then
                 rtabsrt(1,:)=rnode(i,j,k,:)
                 rtabsrt(2,:)=rnode(i+off1,j,k,:)
                 rtabsrt(3,:)=rnode(i+off1,j+off2,k,:)
                 rtabsrt(4,:)=rnode(i,j+off2,k,:)
                 
                 itabsrt(1,:)=[i,j,k]
                 itabsrt(2,:)=[i+off1,j,k]
                 itabsrt(3,:)=[i+off1,j+off2,k]
                 itabsrt(4,:)=[i,j+off2,k]
                 
                 facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4

                 i1=i+(off1-1)/2
                 j1=j+(off2-1)/2
                 
                 if (k.eq.kmax) then
                    irege1=lookind('ZMAX')
                    centre1=facecentre
                 else
                    irege1=iregup(i1,j1,k)
                    centre1=centre(i1,j1,k)
                 endif
                 
                 if (k.eq.0) then
                    irege2=lookind('ZMIN')
                    centre2=facecentre
                 else
                    irege2=iregup(i1,j1,k-1)
                    centre2=centre(i1,j1,k-1)
                 endif

              else if (idir.eq.2) then
                 rtabsrt(1,:)=rnode(i,j,k,:)
                 rtabsrt(2,:)=rnode(i+off1,j,k,:)
                 rtabsrt(3,:)=rnode(i+off1,j,k+off2,:)
                 rtabsrt(4,:)=rnode(i,j,k+off2,:)
                 
                 itabsrt(1,:)=[i,j,k]
                 itabsrt(2,:)=[i+off1,j,k]
                 itabsrt(3,:)=[i+off1,j,k+off2]
                 itabsrt(4,:)=[i,j,k+off2]

                 facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4
                 
                 i1=i+(off1-1)/2
                 k1=k+(off2-1)/2
                 
                 if (j.eq.jmax) then
                    irege1=lookind('YMAX')
                    centre1=facecentre
                 else
                    irege1=iregup(i1,j,k1)
                    centre1=centre(i1,j,k1)
                 endif
                 
                 if (j.eq.0) then
                    irege2=lookind('YMIN')
                    centre2=facecentre
                 else
                    irege2=iregup(i1,j-1,k1)
                    centre2=centre(i1,j-1,k1)
                 endif
              else
                 rtabsrt(1,:)=rnode(i,j,k,:)
                 rtabsrt(2,:)=rnode(i,j+off1,k,:)
                 rtabsrt(3,:)=rnode(i,j+off1,k+off2,:)
                 rtabsrt(4,:)=rnode(i,j,k+off2,:)
                 
                 itabsrt(1,:)=[i,j,k]
                 itabsrt(2,:)=[i,j+off1,k]
                 itabsrt(3,:)=[i,j+off1,k+off2]
                 itabsrt(4,:)=[i,j,k+off2]

                 facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4
                 
                 j1=j+(off1-1)/2
                 k1=k+(off2-1)/2
                 
                 if (i.eq.imax) then
                    irege1=lookind('XMAX')
                    centre1=facecentre
                 else
                    irege1=iregup(i,j1,k1)
                    centre1=centre(i,j1,k1)
                 endif
                 
                 if (i.eq.0) then
                    irege2=lookind('XMIN')
                    centre2=facecentre
                 else
                    irege2=iregup(i-1,j1,k1)
                    centre2=centre(i-1,j1,k1)
                 endif
              endif

! ----------------------------------------------------------------------
! if face type is specified in input
! ----------------------------------------------------------------------

              if (lookup(irege1,irege2).ne.0) then    ! face type must be specified in input file
      
                 l1=1
                 l2=0

                 do ip=1,4
                    if (itabsrt(ip,1).eq.in.and.itabsrt(ip,2).eq.jn.and.itabsrt(ip,3).eq.kn) l2=ip
                 enddo

                 if (l2.ne.0) then       ! face must connect to beta=(in,jn,kn)

                    call getfaceinfo(irege1,irege2,centre1,centre2,itabsrt,rtabsrt,ur,rc,norm)

                    elintNN=elint2d(l1,l2,rtabsrt)

                    do ii=1,nvar
                       do kk=1,nvar
                          
                          ialpha=indQ(i,j,k,ii)
                          if (rmFixed.and.iunk(ialpha)==0) cycle
                          ibeta=indQ(in,jn,kn,kk)
                          
                          tot=qqval(irege1,irege2,ii,kk,ur,rc,norm,iqqx)*elintNN   ! iqqx returned

                          if (.not.(tot==0.0.and.iqqx<1).or.newton=="YES") then   ! non-coincidental zero

                             if (rmFixed.and.iunk(ibeta)==0) then
                                RR(ialpha)=RR(ialpha)-tot*vec(ibeta)
                             else
                                found=.false.
                                do ip=lenQst,lenQ
                                   if (iQ(ip).eq.ialpha.and.jQ(ip).eq.ibeta) then
                                      Qval(ip)=Qval(ip)+tot
                                      found=.true.
                                      exit
                                   endif
                                enddo

                                if (.not.found) then
                                   lenQ=lenQ+1
                                   iQ(lenQ)=ialpha
                                   jQ(lenQ)=ibeta
                                   Qval(lenQ)=tot
                                endif
                             endif
                          endif

                       enddo
                    enddo
         
                 endif   ! face connected to beta
              endif      ! face specified
           endif         ! face exists
        enddo            ! iep
        enddo            ! idir

        endif            ! ntable>0

! ----------------------------------------------------------------------
! end loop over neighbours and nodes
! ----------------------------------------------------------------------

     endif                  ! neighbour node exists
  enddo                     ! loop over in,jn,kn
  enddo
  enddo

  enddo                     ! loop over ijk
  enddo
  enddo

! ----------------------------------------------------------------------
! Write out iusel which records how many times each el used
! and other info
! ----------------------------------------------------------------------

  if (key_db.eq.1) then

     ifirst=1
     do i=0,imax-1
        do j=0,max(0,jmax-1)
           do k=0,max(0,kmax-1)    ! just k=0 for 2D
              
              iusel1=iusel(i,j,k)
              
              if (ifirst.eq.1) then
                 iusmn=iusel1
                 iusmx=iusel1
                 ifirst=0
              else
                 iusmn=min(iusmn,iusel1)
                 iusmx=max(iusmx,iusel1)
              endif
              
           enddo
        enddo
     enddo
     

     do i=0,min(1,imax-1)
        do j=0,jmax-1
           do k=0,kmax-1
              write (1,*) 'iusel',i,j,k,iusel(i,j,k)
           enddo
        enddo
     enddo

     write (1,*) 'el usage=',iusmn,' to',iusmx
  endif
  
end subroutine getQ

! ######################################################################

subroutine getQ2(rmFixed)

! ----------------------------------------------------------------------
! Form Q matrix (iQ,jQ,Qval) and its length lenQ (on common)
! from rnode,iregup,C,aa,ff (on common)
! Alternative routine to form matrix Q by looping over elements
! rather than node pairs
! Set lenQst by considering planes of i=const. We need not revisit
! matrix elements from two planes back or further.
! NOW with nonlinear support!
! ----------------------------------------------------------------------

  use common
  use iofile, only : prterr

  integer endprevrow,lenQst,i,j,k
  logical rmFixed

  if (rmFixed) call prterr('getQ2: rmFixed not supported')

! ----------------------------------------------------------------------
! Loop over elements and form matrix
! ----------------------------------------------------------------------

  print *,'Calculating Q matrix'

  lenQ=0
  lenQst=1

  if (NDIM==1) then
     do i=0,imax-1
        call getQ2x(i,0,0,lenQst)
     enddo

  else if (NDIM==2) then

     do i=0,imax-1
        do j=0,jmax-1
           call getQ2x(i,j,0,lenQst)
        enddo

        if (i /= 0) then
           lenQst=endprevrow+1
        endif
        
        endprevrow=lenQ
     enddo
     
  else if (NDIM==3) then

     do i=0,imax-1

        write (*,*) 'i-plane ',i,' /',imax-1
        
        do j=0,jmax-1
           do k=0,kmax-1
!              print *,i,j,k,lenQ,lenQst
              call getQ2x(i,j,k,lenQst)
           enddo
        enddo

        if (i.ne.0) then
           lenQst=endprevrow+1
        endif
        
        endprevrow=lenQ
     enddo
  endif

end subroutine getQ2

! ######################################################################

subroutine getQ2x(i,j,k,lenQst)

  use geom
  use common
  use indexQ
  use matrices
  
  integer icomb(2,NDIM),itab(2**NDIM,NDIM+1,NDIM),itabsrt(2**NDIM,NDIM)
  double precision rtab(2**NDIM,NDIM),rtabsrt(2**NDIM,NDIM),elint(NDIM,NDIM)
  double precision ur(nvar),dur(nvar,NDIM)
  
  integer lenQst,i,j,k,irege,l1,l2,i1,j1,k1,i2,j2,k2
  integer ii,kk,jj,ll,ialpha,ibeta,ip,iCCx,iaax
  double precision tot,elintNN,rc(NDIM)
  logical found,literal

  j1=0; j2=0         ! if 1D stay zero
  k1=0; k2=0         ! if 1D, 2D stay zero

  if (NDIM==3) then
     irege=iregup(i,j,k)
  else if (NDIM==2) then
     irege=iregup(i,j,0)
  else
     irege=iregup(i,0,0)
  endif
  
  icomb(1,1)=i
  icomb(2,1)=i+1

  if (NDIM>=2) then
     icomb(1,2)=j
     icomb(2,2)=j+1
  endif
  
  if (NDIM==3) then
     icomb(1,3)=k
     icomb(2,3)=k+1
  endif
  
  call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
  call sort(itab,rtab,itabsrt,rtabsrt)

! ----------------------------------------------------------------------
! In non-linear case, we need to get the ur, dur for centre
! of this element and its centroid (xc,yc,zc)
! ----------------------------------------------------------------------

               rc(1)=sum(rtabsrt(:,1))/2**NDIM
  if (NDIM>=2) rc(2)=sum(rtabsrt(:,2))/2**NDIM
  if (NDIM==3) rc(3)=sum(rtabsrt(:,3))/2**NDIM

  if (regCNL(irege).or.regANL(irege)) then
     call getelinfo(itabsrt,rtabsrt,ur,dur)
  endif

! ----------------------------------------------------------------------
! Loop (i1,j1,k1) and (i2,j2,k2) over every pairing of points
! on this hexahedron
! ----------------------------------------------------------------------

  do l1=1,2**NDIM
  do l2=1,2**NDIM
     i1=itabsrt(l1,1)
     i2=itabsrt(l2,1)
     
     if (NDIM>=2) then
        j1=itabsrt(l1,2)
        j2=itabsrt(l2,2)
     endif
     
     if (NDIM==3) then
        k1=itabsrt(l1,3)
        k2=itabsrt(l2,3)
     endif

! ----------------------------------------------------------------------
! Integrate over element using Gauss Quad
! sub jacobian returns jac, invjac and jacdet
! elint(jj,ll)=int_e dN/dxj dN/dxl, elintNN = int_e N N
! where the N's are associated with ijk and in,jn,kn resp.
! ----------------------------------------------------------------------

     call getelint(l1,l2,rtabsrt,elint,elintNN)

! ----------------------------------------------------------------------
! add contributions to matrix elements
! lenQst was set at lenQ+1 inside the ijk, ijkn loops
! We need only search from lenQst to lenQ since we will not access
! previous elements of Q again. Because Q depends on ijk, ijkn
! Matrix element has two contributions for this element,
! One due to C which is called Q in report and one due to a which is P
! Here Q+P -> Q
! To check lenQst, change do ip=1,lenQ and uncomment if(ip.lt.lenQst)..
! ----------------------------------------------------------------------

     do ii=1,nvar
        do kk=1,nvar
           ialpha=indQ(i1,j1,k1,ii)
           ibeta= indQ(i2,j2,k2,kk)
              
           literal=.true.
           tot=0
           do jj=1,NDIM
              do ll=1,NDIM
                 tot=tot+CCval(irege,ii,jj,kk,ll,ur,dur,rc,iCCx)*elint(jj,ll)
                 if (iCCx /=0) literal=.false.
              enddo
           enddo

           tot=tot+aaval(irege,ii,kk,ur,dur,rc,iaax)*elintNN
           if (iaax /= 0) literal=.false.

           if (.not.(tot == 0.and.literal).or.newton=="YES") then
                  
              found=.false.
              do ip=lenQst,lenQ
                 if (iQ(ip).eq.ialpha.and.jQ(ip).eq.ibeta) then
!  (debug)                 if (ip.lt.lenQst) write (*,*) ip,lenQst
                    Qval(ip)=Qval(ip)+tot
                    found=.true.
                    exit
                 endif
              enddo
              
              if (.not.found) then
                 
                 lenQ=lenQ+1
                 iQ(lenQ)=ialpha
                 jQ(lenQ)=ibeta
                 Qval(lenQ)=tot
              endif

           endif
        enddo
     enddo

  enddo                  ! end loop over nodes for given element
  enddo
         
end subroutine getQ2x

! ######################################################################

subroutine getQ3

! ----------------------------------------------------------------------
! Form Q matrix (iQ,jQ,Qval) and its length lenQ (on common)
! from rnode,iregup,C,aa,ff (on common)
! Alternative routine to form matrix Q by looping over elements
! rather than node pairs
! Set lenQst by considering planes of i=const. We need not revisit
! matrix elements from two planes back or further.
! NOW with nonlinear support!
! ----------------------------------------------------------------------

  use geom
  use common
  use indexQ
  use tree
  use matrices
  
  integer icomb(2,3),itab(8,4,3),itabsrt(8,3)
  double precision rtab(8,3),rtabsrt(8,3),elint(3,3)
  double precision ur(nvar),dur(nvar,3),rc(3) ! auto arrays
  
  double precision elintstore(8,8,3,3),elintNNstore(8,8)
  integer iQx,jQx(8),nentry
  double precision Qvalx(8),elintNN

  integer i,j,k,in,jn,kn,irege,l1,l2,ii,jj,kk,ll,ialpha,ibeta,i1,j1,k1
  integer i2,j2,k2,iCCx,iaax
  double precision tot

! ----------------------------------------------------------------------
! Loop over elements and form matrix
! ----------------------------------------------------------------------

  print *,'Calculating Q matrix'

  leniQ=0
  lenjQ=0

  left=0
  right=0
  data=0

  do i=0,imax-1

  write (*,*) 'i-plane ',i,' /',imax-1

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
! In non-linear case, we need to get the ur, dur for centre
! of this element and its centroid rc
! ----------------------------------------------------------------------

     rc(1)=sum(rtabsrt(:,1))/8
     rc(2)=sum(rtabsrt(:,2))/8
     rc(3)=sum(rtabsrt(:,3))/8

     if (regCNL(irege).or.regANL(irege)) then
        call getelinfo(itabsrt,rtabsrt,ur,dur)
     endif

! ----------------------------------------------------------------------
! Loop (i1,j1,k1) and (i2,j2,k2) over every pairing of points
! on this hexahedron
! ----------------------------------------------------------------------

     do l1=1,8
     do l2=1,8
        
! ----------------------------------------------------------------------
! Integrate over element using Gauss Quad
! sub jacobian returns jac, invjac and jacdet
! elint(jj,ll)=int_e dN/dxj dN/dxl, elintNN = int_e N N
! where the N's are associated with ijk and in,jn,kn resp.
! ----------------------------------------------------------------------

        call getelint(l1,l2,rtabsrt,elint,elintNN)
        
        elintstore(l1,l2,:,:)=elint(:,:)
        elintNNstore(l1,l2)=elintNN

     enddo
     enddo
        
! ----------------------------------------------------------------------
! add contributions to matrix elements
! lenQst was set at lenQ+1 inside the ijk, ijkn loops
! We need only search from lenQst to lenQ since we will not access
! previous elements of Q again. Because Q depends on ijk, ijkn
! Matrix element has two contributions for this element,
! One due to C which is called Q in report and one due to a which is P
! Here Q+P -> Q
! To check lenQst, change do ip=1,lenQ and uncomment if(ip.lt.lenQst)..
! ----------------------------------------------------------------------

     do ii=1,nvar
        do kk=1,nvar

           do l1=1,8

              i1=itabsrt(l1,1)
              j1=itabsrt(l1,2)
              k1=itabsrt(l1,3)

              ialpha=indQ(i1,j1,k1,ii)
              iQx=ialpha

              nentry=0
              do l2=1,8

                 i2=itabsrt(l2,1)
                 j2=itabsrt(l2,2)
                 k2=itabsrt(l2,3)

                 ibeta= indQ(i2,j2,k2,kk)
              
                 tot=0
                 do jj=1,3
                    do ll=1,3
                       tot=tot+CCval(irege,ii,jj,kk,ll,ur,dur,rc,iCCx)*elint(jj,ll)
                    enddo
                 enddo

                 tot=tot+aaval(irege,ii,kk,ur,dur,rc,iaax)*elintNN

                 if (tot.ne.0) then
                    nentry=nentry+1
                    jQx(nentry)=ibeta
                    Qvalx(nentry)=tot
                 endif

              enddo               ! loop over l2
              if (nentry.gt.0) call insertQ(nentry,iQx,jQx,Qvalx)
           enddo                  ! loop over l1

        enddo
     enddo                        ! loop over ii and kk
         
  enddo                     
  enddo
  enddo                     ! end loop over i,j,k elements

end subroutine getQ3

! ######################################################################

subroutine getR(append)

! ----------------------------------------------------------------------
! Form the R vector. Uses rnode, ff, iff, u (non-linear)
! ----------------------------------------------------------------------

  use common
  use geom
  use shape
  use indexQ
  use matrices
  use util
  use iofile, only : prterr

  logical append
  
  integer icomb(2,NDIM),ies(8),jes(8),kes(8),itab(2**NDIM,4,NDIM)
  integer itabsrt(2**NDIM,NDIM)

  double precision rtab(2**NDIM,NDIM),rtabsrt(2**NDIM,NDIM)
  double precision jac(NDIM,NDIM),invjac(NDIM,NDIM),jacdet
  double precision ur(nvar),dur(nvar,NDIM),rc(NDIM),norm(3)

  logical iein,jein,kein,elin
  double precision elintN,xi,eta,mu,trailer

  integer i,j,k,iep,ie,je,ke,i1,j1,k1,irege,ip,ixi,ieta,imu,ii,l1
  integer ialpha,inode,jnode,knode,l,ic,off1,off2
  integer off1s(4),off2s(4),idir,irege1,irege2
  double precision drdxi(3),drdeta(3)
  double precision centre1(3),centre2(3),facecentre(3)

  save ies,jes,kes,off1s,off2s

  data ies/ 1,-1, 1,-1, 1,-1, 1,-1/  
  data jes/ 1, 1,-1,-1, 1, 1,-1,-1/
  data kes/ 1, 1, 1, 1,-1,-1,-1,-1/

  data off1s/1, 1,-1,-1/
  data off2s/1,-1, 1,-1/

! ----------------------------------------------------------------------
! Loop over nodes. getQ must have been called already so RR has non-zero
! entries by now
! ----------------------------------------------------------------------

  print *,'Calculating R matrix'

  if (.not.append) RR=0            ! start from scratch

  do k=0,kmax
  write (*,*) 'k-plane:',k,'/',kmax
  do j=0,jmax
  do i=0,imax

     do iep=1,2**NDIM
        ie=ies(iep)
        je=jes(iep)
        ke=kes(iep)
            
        iein=(i+ie.ge.0.and.i+ie.le.imax)
        jein=(j+je.ge.0.and.j+je.le.jmax).or.NDIM==1
        kein=(k+ke.ge.0.and.k+ke.le.kmax).or.NDIM<=2
        
! ----------------------------------------------------------------------
! get region number of element
! if ie = 1, i1=i; if ie=-1, i1=i-1
! If connecting element exists, get itabsrt, rtabsrt to define it
! Get l1, local index number of node (i,j,k)
! ----------------------------------------------------------------------

        if (iein.and.jein.and.kein) then
               
           i1=i+(ie-1)/2
           j1=j+(je-1)/2
           k1=k+(ke-1)/2
           irege=iregup(i1,j1,k1)

           icomb(1,1)=i
           icomb(2,1)=i+ie

           if (NDIM>=2) then
              icomb(1,2)=j
              icomb(2,2)=j+je
           endif

           if (NDIM==3) then
              icomb(1,3)=k
              icomb(2,3)=k+ke
           endif

           call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
           call sort(itab,rtab,itabsrt,rtabsrt)
           
           l1=0
           
           if (NDIM==3) then
              do ip=1,8
                 inode=itabsrt(ip,1)
                 jnode=itabsrt(ip,2)
                 knode=itabsrt(ip,3)
                 if (inode.eq.i.and.jnode.eq.j.and.knode.eq.k) l1=ip
              enddo
           else if (NDIM==2) then
              do ip=1,4
                 inode=itabsrt(ip,1)
                 jnode=itabsrt(ip,2)
                 if (inode.eq.i.and.jnode.eq.j) l1=ip
              enddo
           else
              do ip=1,2
                 inode=itabsrt(ip,1)
                 if (inode.eq.i) l1=ip
              enddo
           endif

           if (l1==0) call prterr('Error in getR, bad l1')

! ----------------------------------------------------------------------
! Calculate elintN = int N_alpha^(e) dV, alpha -> l1
! for given node and element
! ----------------------------------------------------------------------

           elintN=0
           
           do ixi=1,ng_xi
              do ieta=1,ng_eta
                 do imu=1,ng_mu          ! ng_mu=1 for 2D
                    xi=gauss(ng_xi,ixi)
                    eta=gauss(ng_eta,ieta)
                    mu=gauss(ng_mu,imu)
                    
                    call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet)
                    
                    if (NDIM==3) then
                       trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)*wt(ng_mu,imu)
                    else if (NDIM==2) then
                       trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)
                    else
                       trailer=jacdet*wt(ng_xi,ixi)
                    endif

                    elintN=elintN+Nl(l1,xi,eta,mu)*trailer
                    
                 enddo
              enddo
           enddo
           
! ----------------------------------------------------------------------
! In non-linear case, we need to get the ur, dur for centre
! of this element
! ----------------------------------------------------------------------

                        rc(1)=sum(rtabsrt(:,1))/2**NDIM
           if (NDIM>=2) rc(2)=sum(rtabsrt(:,2))/2**NDIM
           if (NDIM==3) rc(3)=sum(rtabsrt(:,3))/2**NDIM

           if (regFNL(irege)) then
              call getelinfo(itabsrt,rtabsrt,ur,dur)
           endif

! ----------------------------------------------------------------------
! Add this elements contribution to RR(ii,i,j,k)->RR(ialpha)
! ----------------------------------------------------------------------
               
           do ii=1,nvar
              ialpha=indQ(i,j,k,ii)
              RR(ialpha)=RR(ialpha)+ffval(irege,ii,ur,dur,rc)*elintN
           enddo
        endif               ! if element exists
     enddo                  ! loop over connecting elements

! ----------------------------------------------------------------------
! Add surface G term. Get 12 faces
! idir=1: i,j. idir=2 i,k. idir=3 j,k
! ----------------------------------------------------------------------

     if (ntable > 0) then

     do idir=1,3
     do iep=1,4
   
        off1=off1s(iep)
        off2=off2s(iep)

        if (idir.eq.1) then
           elin=i+off1.le.imax.and.i+off1.ge.0.and.&
                j+off2.le.jmax.and.j+off2.ge.0
        else if (idir.eq.2) then
           elin=i+off1.le.imax.and.i+off1.ge.0.and.&
                k+off2.le.kmax.and.k+off2.ge.0
        else
           elin=j+off1.le.jmax.and.j+off1.ge.0.and.&
                k+off2.le.kmax.and.k+off2.ge.0
        endif

        if (elin) then            ! face must exist

           if (idir.eq.1) then
              rtabsrt(1,:)=rnode(i,j,k,:)
              rtabsrt(2,:)=rnode(i+off1,j,k,:)
              rtabsrt(3,:)=rnode(i+off1,j+off2,k,:)
              rtabsrt(4,:)=rnode(i,j+off2,k,:)

              itabsrt(1,:)=[i,j,k]
              itabsrt(2,:)=[i+off1,j,k]
              itabsrt(3,:)=[i+off1,j+off2,k]
              itabsrt(4,:)=[i,j+off2,k]
              
              facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4

              i1=i+(off1-1)/2
              j1=j+(off2-1)/2
              
              if (k.eq.kmax) then
                 irege1=lookind('ZMAX')
                 centre1=facecentre
              else
                 irege1=iregup(i1,j1,k)
                 centre1=centre(i1,j1,k)
              endif
              
              if (k.eq.0) then
                 irege2=lookind('ZMIN')
                 centre2=facecentre
              else
                 irege2=iregup(i1,j1,k-1)
                 centre2=centre(i1,j1,k-1)
              endif

           else if (idir.eq.2) then
              rtabsrt(1,:)=rnode(i,j,k,:)
              rtabsrt(2,:)=rnode(i+off1,j,k,:)
              rtabsrt(3,:)=rnode(i+off1,j,k+off2,:)
              rtabsrt(4,:)=rnode(i,j,k+off2,:)

              itabsrt(1,:)=[i,j,k]
              itabsrt(2,:)=[i+off1,j,k]
              itabsrt(3,:)=[i+off1,j,k+off2]
              itabsrt(4,:)=[i,j,k+off2]

              facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4

              i1=i+(off1-1)/2
              k1=k+(off2-1)/2
              
              if (j.eq.jmax) then
                 irege1=lookind('YMAX')
                 centre1=facecentre
              else
                 irege1=iregup(i1,j,k1)
                 centre1=centre(i1,j,k1)
              endif
              
              if (j.eq.0) then
                 irege2=lookind('YMIN')
                 centre2=facecentre
              else
                 irege2=iregup(i1,j-1,k1)
                 centre2=centre(i1,j-1,k1)
              endif
           else
              rtabsrt(1,:)=rnode(i,j,k,:)
              rtabsrt(2,:)=rnode(i,j+off1,k,:)
              rtabsrt(3,:)=rnode(i,j+off1,k+off2,:)
              rtabsrt(4,:)=rnode(i,j,k+off2,:)

              itabsrt(1,:)=[i,j,k]
              itabsrt(2,:)=[i,j+off1,k]
              itabsrt(3,:)=[i,j+off1,k+off2]
              itabsrt(4,:)=[i,j,k+off2]
              
              facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4

              j1=j+(off1-1)/2
              k1=k+(off2-1)/2

              if (i.eq.imax) then
                 irege1=lookind('XMAX')
                 centre1=facecentre
              else
                 irege1=iregup(i,j1,k1)
                 centre1=centre(i,j1,k1)
              endif
              
              if (i.eq.0) then
                 irege2=lookind('XMIN')
                 centre2=facecentre
              else
                 irege2=iregup(i-1,j1,k1)
                 centre2=centre(i-1,j1,k1)
              endif
           endif

! ----------------------------------------------------------------------
! if face specified, add R contribution
! ----------------------------------------------------------------------

           if (lookup(irege1,irege2).ne.0) then    ! face type must be specified in input file
              
              elintN=0
              do ixi=1,ng_xi
                 do ieta=1,ng_eta
        
                    xi=gauss(ng_xi,ixi)
                    eta=gauss(ng_eta,ieta)
        
                    drdxi=0
                    drdeta=0

                    do l=1,4
                       do ic=1,3
                          drdxi(ic) =drdxi(ic) +rtabsrt(l,ic)* dN2d_dxi(l,xi,eta)
                          drdeta(ic)=drdeta(ic)+rtabsrt(l,ic)*dN2d_deta(l,xi,eta)
                       enddo
                    enddo
        
                    trailer=mag(cross(drdxi,drdeta))*wt(ng_xi,ixi)*wt(ng_eta,ieta)
        
                    elintN=elintN+N2d(1,xi,eta)*trailer
                 enddo
              enddo
              
              call getfaceinfo(irege1,irege2,centre1,centre2,itabsrt,rtabsrt,ur,rc,norm)
              
              do ii=1,nvar
                 ialpha=indQ(i,j,k,ii)
                 RR(ialpha)=RR(ialpha)+ggval(irege1,irege2,ii,ur,rc,norm)*elintN
              enddo
              
           endif      ! face specified
        endif         ! face exists
     enddo
     enddo

     endif            ! ntable > 0
      
  enddo                     ! loop over nodes
  enddo
  enddo

end subroutine getR

! ######################################################################

subroutine getk

! ----------------------------------------------------------------------
! Form the kfac vector for transient. Uses rnode, kfac
! Note that this is very similar to getR, kfac entries for same
! node but different variable number are the same.
! ----------------------------------------------------------------------

  use common
  use geom
  use shape
  use indexQ

  integer icomb(2,3),ies(8),jes(8),kes(8),itab(8,4,3),itabsrt(8,3)

  double precision rtab(8,3),rtabsrt(8,3)
  double precision xi,eta,mu,jac(3,3),invjac(3,3),jacdet

  logical iein,jein,kein

  integer i,j,k,iep,ie,je,ke,ip,ixi,ieta,imu,ii,ialpha
  integer inode,jnode,knode,l1
  double precision elintN,trailer
  
  save ies,jes,kes
  
  data ies/ 1, 1, 1, 1,-1,-1,-1,-1/
  data jes/ 1, 1,-1,-1, 1, 1,-1,-1/
  data kes/ 1,-1, 1,-1, 1,-1, 1,-1/

! ----------------------------------------------------------------------
! Zero kk vector and loop over nodes
! ----------------------------------------------------------------------

  print *,'Calculating transient coefficients'

  do i=1,ndof
     kfac(i)=0
  enddo

  do k=0,kmax
  write (*,*) 'k-plane:',k,'/',kmax
  do j=0,jmax
  do i=0,imax

     do iep=1,8
        ie=ies(iep)
        je=jes(iep)
        ke=kes(iep)
            
        iein=(i+ie.ge.0.and.i+ie.le.imax)
        jein=(j+je.ge.0.and.j+je.le.jmax)
        kein=(k+ke.ge.0.and.k+ke.le.kmax)
        
! ----------------------------------------------------------------------
! get region number of element
! if ie = 1, i1=i; if ie=-1, i1=i-1
! If connecting element exists, get itabsrt, rtabsrt to define it
! Get l1, local index number of node (i,j,k)
! ----------------------------------------------------------------------

        if (iein.and.jein.and.kein) then
           icomb(1,1)=i
           icomb(1,2)=j
           icomb(1,3)=k
           
           icomb(2,1)=i+ie
           icomb(2,2)=j+je
           icomb(2,3)=k+ke
           
           call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
           call sort(itab,rtab,itabsrt,rtabsrt)
           
           do ip=1,8
              inode=itabsrt(ip,1)
              jnode=itabsrt(ip,2)
              knode=itabsrt(ip,3)
              if (inode.eq.i.and.jnode.eq.j.and.knode.eq.k) l1=ip
           enddo

! ----------------------------------------------------------------------
! Calculate elintN = int N_alpha^(e) dV, alpha -> l1
! for given node and element
! ----------------------------------------------------------------------

           elintN=0
           
           do ixi=1,ng_xi
              do ieta=1,ng_eta
                 do imu=1,ng_mu
                    xi=gauss(ng_xi,ixi)
                    eta=gauss(ng_eta,ieta)
                    mu=gauss(ng_mu,imu)
                    
                    call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet) 
                    trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)*wt(ng_mu,imu)
                    
                    elintN=elintN+Nl(l1,xi,eta,mu)*trailer
                    
                 enddo
              enddo
           enddo
           
! ----------------------------------------------------------------------
! Add this elements contribution to kk(i,j,k)->RR(ialpha)
! ----------------------------------------------------------------------
               
           do ii=1,nvar
              ialpha=indQ(i,j,k,ii)
              kfac(ialpha)=kfac(ialpha)+elintN
           enddo
        endif               ! if element exists
     enddo                  ! loop over connecting elements
     
  enddo                     ! loop over nodes
  enddo
  enddo

end subroutine getk

! ######################################################################
      
function energy(enreg)

! ----------------------------------------------------------------------
! Return energy for given solution in global "vec"
! also return energy of element regions, enreg without using Q. 
! Can be compared to 0.5 vec^T Q vec
! Uses u,C,iregup,ijkmax,gauss,wt from common module
! calls jacobian
! Need du(nvar,3). du(1,2)=du1/dx2 etc
!
! u(l,ii) is u_ii at local node l
! U=sum(el) int sum 0.5 ui,j Cijkl uk,l dV
!  +sum(el) int sum 0.5 ui aik uk
!  -sum(el) int sum fi ui
!
! ui=sum_local ui(local)*N(local)
! ui,j=sum_local ui(local)*dN(local)/dxj
! Note: several local arrays are automatic
! ----------------------------------------------------------------------

  use common
  use geom
  use shape
  use indexQ
  use matrices
  use util
  use iofile, only : prterr

  double precision energy,enreg(:) ! nreg on module common

  integer itab(2**NDIM,NDIM+1,NDIM),itabsrt(2**NDIM,NDIM),icomb(2,NDIM),ic,idir,ind,irege1,irege2
  double precision rtab(2**NDIM,NDIM),rtabsrt(2**NDIM,NDIM)
  double precision ul(2**NDIM,nvar),dN(NDIM),du(nvar,NDIM),ur(nvar)
  double precision r0(3),u0(nvar),du0(nvar,NDIM),facecentre(3),centre1(3),centre2(3),norm0(3)

  double precision Nl1,drdxi(3),drdeta(3)
  double precision jac(NDIM,NDIM),invjac(NDIM,NDIM),jacdet

  double precision Ce(nvar,NDIM,nvar,NDIM),ae(nvar,nvar),fe(nvar)  ! auto arrays

  double precision en1,xi,eta,mu,trailer
  integer i,j,k,irege,ixi,ieta,imu,l,ii,jj,kk,ll,inode,jnode,knode
  integer itop,jtop,ktop,iCCx,iaax,iqqx

  energy=0
  enreg=0

  if (NDIM==3) then
     itop=imax-1
     jtop=jmax-1
     ktop=kmax-1
  else if (NDIM==2) then
     itop=imax-1
     jtop=jmax-1
     ktop=0
  else
     itop=imax-1
     jtop=0
     ktop=0
  endif

  do i=0,itop
  do j=0,jtop
  do k=0,ktop

! ----------------------------------------------------------------------
! Get region number and local C and a mat props
! ----------------------------------------------------------------------

     irege=iregup(i,j,k)
     
     Ce=CC(irege,:,:,:,:)
     ae=aa(irege,:,:)
     fe=ff(irege,:)

! ----------------------------------------------------------------------
! Get the sorted element nodes. Get local u values in ul(l,ii)
! ----------------------------------------------------------------------

     icomb(1,1)=i
     icomb(2,1)=i+1

     if (NDIM>=2) then
        icomb(1,2)=j
        icomb(2,2)=j+1
     endif

     if (NDIM==3) then
        icomb(1,3)=k
        icomb(2,3)=k+1
     endif

     call getel(icomb,rnode,imax,jmax,kmax,itab,rtab)
     call sort(itab,rtab,itabsrt,rtabsrt)
         
     do l=1,2**NDIM
        inode=itabsrt(l,1)
        if (NDIM>=2) jnode=itabsrt(l,2)
        if (NDIM==3) knode=itabsrt(l,3)
        
        do ii=1,nvar
           ul(l,ii)=vec(indQ(inode,jnode,knode,ii))
        enddo
     enddo

! ----------------------------------------------------------------------
! for non-linear elements, get element info
! ----------------------------------------------------------------------

                  r0(1)=sum(rtabsrt(:,1))/2**NDIM
     if (NDIM>=2) r0(2)=sum(rtabsrt(:,2))/2**NDIM
     if (NDIM==3) r0(3)=sum(rtabsrt(:,3))/2**NDIM

     if (regCNL(irege).or.regANL(irege)) then
        call getelinfo(itabsrt,rtabsrt,u0,du0)
     endif

! ----------------------------------------------------------------------
! Element integral loops
! Integrate over element, get du(ii,jj), ur(ii)
! ----------------------------------------------------------------------

     do ixi=1,ng_xi
     do ieta=1,ng_eta
     do imu=1,ng_mu
        xi=gauss(ng_xi,ixi)
        eta=gauss(ng_eta,ieta)
        mu=gauss(ng_mu,imu)

        call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet)

        if (NDIM==3) then
           trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)*wt(ng_mu,imu)
        else if (NDIM==2) then
           trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)
        else
           trailer=jacdet*wt(ng_xi,ixi)
        endif

        du(:,:)=0
            
        do l=1,2**NDIM
           call getdN(l,xi,eta,mu,invjac,dN)
           do ii=1,nvar
              do jj=1,NDIM
                 du(ii,jj)=du(ii,jj)+ul(l,ii)*dN(jj)
              enddo
           enddo
        enddo

        ur(:)=0
            
        do l=1,2**NDIM
           Nl1=Nl(l,xi,eta,mu)
           do ii=1,nvar
              ur(ii)=ur(ii)+ul(l,ii)*Nl1
           enddo
        enddo

! ----------------------------------------------------------------------
! Calculate contributions (times 2) due to C and a terms for element
! ----------------------------------------------------------------------

        en1=0
        do ii=1,nvar
           do jj=1,NDIM
              do kk=1,nvar
                 do ll=1,NDIM
                    en1=en1+0.5*du(ii,jj)*CCval(irege,ii,jj,kk,ll,u0,du0,r0,iCCx)*du(kk,ll)*trailer
                 enddo
              enddo
           enddo
        enddo

        do ii=1,nvar
           do kk=1,nvar              
              en1=en1+0.5*ur(ii)*aaval(irege,ii,kk,u0,du0,r0,iaax)*ur(kk)*trailer
           enddo
        enddo
        
        do ii=1,nvar
           en1=en1-ffval(irege,ii,u0,du0,r0)*ur(ii)*trailer
        enddo

        enreg(irege)=enreg(irege)+en1
        energy=energy+en1

     enddo
     enddo
     enddo                     ! element integral

  enddo
  enddo
  enddo                     ! sum of all elements

! ----------------------------------------------------------------------
! energy in surfaces, directions
! idir=1 i,j. idir=2 i,k. idir=3 j,k
! ----------------------------------------------------------------------

  if (ntable > 0) then

  if (NDIM /= 3) call prterr('SURFACE only supported for 3D')

  do idir=1,3

  if (idir.eq.1) then
     itop=imax-1
     jtop=jmax-1
     ktop=kmax
  else if (idir.eq.2) then
     itop=imax-1
     ktop=kmax-1
     jtop=jmax
  else
     jtop=jmax-1
     ktop=kmax-1
     itop=imax
  endif

! ----------------------------------------------------------------------
! Loop over all faces (kface, jface, iface). get face geometry
! ----------------------------------------------------------------------

  do i=0,itop
  do j=0,jtop
  do k=0,ktop

     if (idir.eq.1) then
        rtabsrt(1,:)=rnode(i,j,k,:)
        rtabsrt(2,:)=rnode(i+1,j,k,:)
        rtabsrt(3,:)=rnode(i+1,j+1,k,:)
        rtabsrt(4,:)=rnode(i,j+1,k,:)
        
        itabsrt(1,:)=[i,j,k]
        itabsrt(2,:)=[i+1,j,k]
        itabsrt(3,:)=[i+1,j+1,k]
        itabsrt(4,:)=[i,j+1,k]

        facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4
        
        if (k.eq.kmax) then
           irege1=lookind('ZMAX')
           centre1=facecentre
        else
           irege1=iregup(i,j,k)
           centre1=centre(i,j,k)
        endif

        if (k.eq.0) then
           irege2=lookind('ZMIN')
           centre2=facecentre
        else
           irege2=iregup(i,j,k-1)
           centre2=centre(i,j,k-1)
        endif

     else if (idir.eq.2) then
        rtabsrt(1,:)=rnode(i,j,k,:)
        rtabsrt(2,:)=rnode(i+1,j,k,:)
        rtabsrt(3,:)=rnode(i+1,j,k+1,:)
        rtabsrt(4,:)=rnode(i,j,k+1,:)
        
        itabsrt(1,:)=[i,j,k]
        itabsrt(2,:)=[i+1,j,k]
        itabsrt(3,:)=[i+1,j,k+1]
        itabsrt(4,:)=[i,j,k+1]

        facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4

        if (j.eq.jmax) then
           irege1=lookind('YMAX')
           centre1=facecentre
        else
           irege1=iregup(i,j,k)
           centre1=centre(i,j,k)
        endif

        if (j.eq.0) then
           irege2=lookind('YMIN')
           centre2=facecentre
        else
           irege2=iregup(i,j-1,k)
           centre2=centre(i,j-1,k)
        endif
     else
        rtabsrt(1,:)=rnode(i,j,k,:)
        rtabsrt(2,:)=rnode(i,j+1,k,:)
        rtabsrt(3,:)=rnode(i,j+1,k+1,:)
        rtabsrt(4,:)=rnode(i,j,k+1,:)
        
        itabsrt(1,:)=[i,j,k]
        itabsrt(2,:)=[i,j+1,k]
        itabsrt(3,:)=[i,j+1,k+1]
        itabsrt(4,:)=[i,j,k+1]

        facecentre=(rtabsrt(1,:)+rtabsrt(2,:)+rtabsrt(3,:)+rtabsrt(4,:))/4

        if (i.eq.imax) then
           irege1=lookind('XMAX')
           centre1=facecentre
        else
           irege1=iregup(i,j,k)
           centre1=centre(i,j,k)
        endif

        if (i.eq.0) then
           irege2=lookind('XMIN')
           centre2=facecentre
        else
           irege2=iregup(i-1,j,k)
           centre2=centre(i-1,j,k)
        endif
     endif

! ----------------------------------------------------------------------
! if face has energy, get local node solutions ul
! get trailer and integrate energy
! ----------------------------------------------------------------------

     ind=lookup(irege1,irege2)
     
     if (ind.ne.0) then

        do l=1,4
           inode=itabsrt(l,1)
           jnode=itabsrt(l,2)
           knode=itabsrt(l,3)
           
           do ii=1,nvar
              ul(l,ii)=vec(indQ(inode,jnode,knode,ii))
           enddo
        enddo

        call getfaceinfo(irege1,irege2,centre1,centre2,itabsrt,rtabsrt,u0,r0,norm0)

        do ixi=1,ng_xi
        do ieta=1,ng_eta
           xi=gauss(ng_xi,ixi)
           eta=gauss(ng_eta,ieta)

           drdxi=0
           drdeta=0

           do ic=1,3
              do l=1,4
                 drdxi(ic)= drdxi(ic)+ rtabsrt(l,ic)*dN2d_dxi(l,xi,eta)
                 drdeta(ic)=drdeta(ic)+rtabsrt(l,ic)*dN2d_deta(l,xi,eta)
              enddo
           enddo

           trailer=mag(cross(drdxi,drdeta))*wt(ng_xi,ixi)*wt(ng_eta,ieta)

           ur=0
           do ii=1,nvar
              do l=1,4
                 ur(ii)=ur(ii)+ul(l,ii)*N2d(l,xi,eta)
              enddo
           enddo

           do ii=1,nvar
              do kk=1,nvar
                 energy=energy+0.5*ur(ii)*qqval(irege1,irege2,ii,kk,u0,r0,norm0,iqqx)*ur(kk)*trailer    ! previously qq(ind,ii,kk)
              enddo
           enddo
           
           do ii=1,nvar
              energy=energy-ur(ii)*ggval(irege1,irege2,ii,u0,r0,norm0)*trailer                     ! previously gg(ind,ii)
           enddo
        enddo
        enddo
     endif

  enddo
  enddo
  enddo

enddo

endif  ! ntable > 0

end function energy

! ######################################################################

subroutine getelint(l1,l2,rtabsrt,elint,elintNN)

  use common
  use geom
  use shape

  integer l1,l2
  double precision rtabsrt(2**NDIM,NDIM),elint(NDIM,NDIM),elintNN,dN1(NDIM),dN2(NDIM)
  integer ixi,ieta,imu
  double precision xi,eta,mu,trailer,jac(NDIM,NDIM),invjac(NDIM,NDIM),jacdet

  integer jj,ll

  elintNN=0
  elint(:,:)=0
                 
  do ixi=1,ng_xi
     do ieta=1,ng_eta
        do imu=1,ng_mu
           xi=gauss(ng_xi,ixi)
           eta=gauss(ng_eta,ieta)
           mu=gauss(ng_mu,imu)
           
           call jacobian(xi,eta,mu,rtabsrt,jac,invjac,jacdet) 
           call getdN(l1,xi,eta,mu,invjac,dN1)
           call getdN(l2,xi,eta,mu,invjac,dN2)
                    
           if (NDIM==3) then
              trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)*wt(ng_mu,imu)
           else if (NDIM==2) then
              trailer=jacdet*wt(ng_xi,ixi)*wt(ng_eta,ieta)
           else
              trailer=jacdet*wt(ng_xi,ixi)
           endif

           elintNN=elintNN &
                +Nl(l1,xi,eta,mu)*Nl(l2,xi,eta,mu)*trailer
           
           do jj=1,NDIM
              do ll=1,NDIM
                 elint(jj,ll)=elint(jj,ll)+dN1(jj)*dN2(ll)*trailer
              enddo
           enddo
           
        enddo
     enddo
  enddo
  
end subroutine getelint

! ######################################################################
     
function elint2d(l1,l2,rtabsrt)

  use common
  use shape
  use util

  double precision rtabsrt(8,3)
  integer l1,l2
  
  double precision elint2d,xi,eta,drdxi(3),drdeta(3),trailer
  integer ixi,ieta,l,ic

  elint2d=0
  do ixi=1,ng_xi
     do ieta=1,ng_eta
        
        xi=gauss(ng_xi,ixi)
        eta=gauss(ng_eta,ieta)
        
        drdxi=0
        drdeta=0
        do l=1,4
           do ic=1,3
              drdxi(ic) =drdxi(ic) +rtabsrt(l,ic)* dN2d_dxi(l,xi,eta)
              drdeta(ic)=drdeta(ic)+rtabsrt(l,ic)*dN2d_deta(l,xi,eta)
           enddo
        enddo
        
        trailer=mag(cross(drdxi,drdeta))*wt(ng_xi,ixi)*wt(ng_eta,ieta)
        
        elint2d=elint2d+N2d(l1,xi,eta)*N2d(l2,xi,eta)*trailer
     enddo
  enddo
end function elint2d

end module qpr
