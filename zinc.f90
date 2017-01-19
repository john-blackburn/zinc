! Part of Zinc FE package. Author: John Blackburn

program zinc

!----------------------------------------------------------------------
! Zinc Is Not Comsol! (much better)
!
! Program to simulate problem:
! -(Cijkl uk,l),j+aik uk=fi; nj(Cijkl uk,l)+qik uk=gi boundary
!
! Reads mesh info from Metamesh output (MTF file)
!
! main, getQ, getQ2 and getR are the only subprogs
! to write to module common in 'common.f90'
!
! Q matrix is ndof*ndof, ndof=nvar*nnod
! (stored sparse: about nvar**2*nnod*27 entries in Qval,iQ,jQ)
! vec, iunk, irowst, irowed, RR are ndof long
! ndof=nvar*nnod, nnod=(imax+1)*(jmax+1)*(kmax+1)
!
! Program input: file.zin, file.mtf (from zmesh via file.min), 
! file.rst, file.con, file.dll (opened using LoadLibrary)
! Output: file.zou, file.zls
! use defparam in .con file and for initialization
! uses evalexpr for initialization (could replace with evaltoken)
! The only modules that use evaluate module are qpr.f90 and iofile.f90

! Terminology: In theory manual we have p(alpha,i) with alpha global
! number of node and i var number.
! Here we call p, alpha. I.e. alpha goes over all degress of freedom.
! Same for beta.
! var numbers i,j,k,l are referred as ii,jj,kk,ll
!
! Now extended to 1D and 2D. The rule is that each subprogram is smart
! enough to do the right thing based on NDIM. Caller never needs to
! use getelint2D, getelint3D etc. (and no function overloading) 
! Typically subprograms will simply
! have different loop limits based on ndim but some if statements
! will be needed to avoid unnecessary calculations. 
! Some subprograms will need to be substantially rewritten 3 ways
! Where functions have the call sequence (i,j,k,...) etc this will stay but 
! eg k will be ignored. Actual array sizes are exact not "big enough"
! eg invjac(ndim,ndim) and dur(nvar,ndim). Use array bounds check for bugs
! ----------------------------------------------------------------------

  use common
  use strings, only : lowercase
  use evaluate
  use qpr
  use indexQ
  use update
  use iofile
  use interpolate
  use heap
  use matrices
  use solvers
  use jacobianmod
  use cuboid

  implicit none

  double precision redge(3,2),r(3)
  double precision, allocatable :: enreg(:)
  
  integer key_sym,istatus,length,ios,term
  integer i,j,k,ic,jc,kc,mreg,mregup,ireg1,nunk,ind,ialpha
  integer ifirst,iQ1,jQ1,ip,ii,itot,mtf_nregel,mtf_nregnd
  integer i1,j1,k1,ii1,lfname,allocerr(30),tind,kk

  character(1000) currfile,line
  character(20) temp,cdum

  double precision val,x,y,z,enQ,en,u1,err,resid0,s
  double precision resid,lhs,rhs,residDiag,lhsDiag,rhsDiag
  double precision lhsI,rhsI,lhsIDiag,rhsIDiag,pivot

  integer iter,ierrslap,ndofx,iback
  real t1,t2
  logical jactest,found,backtrack

! ----------------------------------------------------------------------
! Set Gaussian integration points and weights
! ----------------------------------------------------------------------

  gauss(1,1)= 0d0;                 wt(1,1)=2d0

  gauss(2,1)=-0.577350269189626d0; wt(2,1)=1d0                   ! -sqrt(1/3)
  gauss(2,2)= 0.577350269189626d0; wt(2,2)=1d0

  gauss(3,1)= 0d0;                 wt(3,1)=0.888888888888889d0   ! 8/9
  gauss(3,2)=-0.774596669241483d0; wt(3,2)=0.555555555555556d0   ! -sqrt(3/5), 5/9
  gauss(3,3)= 0.774596669241483d0; wt(3,3)=0.555555555555556d0

  gauss(4,1)=-0.339981043584856d0; wt(4,1)=0.652145154862546d0
  gauss(4,2)= 0.339981043584856d0; wt(4,2)=0.652145154862546d0
  gauss(4,3)=-0.861136311594053d0; wt(4,3)=0.347854845137454d0
  gauss(4,4)= 0.861136311594053d0; wt(4,4)=0.347854845137454d0

  gauss(5,1)= 0d0;                 wt(5,1)=0.568888888888889d0
  gauss(5,2)=-0.538469310105683d0; wt(5,2)=0.478628670499366d0
  gauss(5,3)= 0.538469310105683d0; wt(5,3)=0.478628670499366d0
  gauss(5,4)=-0.906179845938664d0; wt(5,4)=0.236926885056189d0
  gauss(5,5)= 0.906179845938664d0; wt(5,5)=0.236926885056189d0

! ----------------------------------------------------------------------
! These parameters are set here and cannot be set by the user
! Later we might make these available to the user
! ----------------------------------------------------------------------

  centre_eval=.true.            ! true = evaluate mat props at centre of element
  jac_by_element=.true.         ! true = loop over elements not node triplets
  jactest=.false.               ! true = also calculate jacobian by perturbation
  key_sym=0                     ! set key_sym=0 for now. Ie, do not symmetrize
  backtrack=.false.

! ----------------------------------------------------------------------
! Get single command prompt. Strip .zin off end if found
! ----------------------------------------------------------------------

  call get_command_argument(1,fname,length,istatus)

  if (istatus.gt.0) then
     print *,'Unable to retrieve command line argument!'
     call exit(1)
  else if (istatus.eq.-1) then
     print *,'Filename too long, max 1000 characters'
     call exit(1)
  endif

  fname=lowercase(fname)
  lfname=len_trim(fname)

  if (lfname.eq.0) then
     print *,name,' Error: no input file specified!'
     call exit(1)
  endif

  do i=lfname,1,-1
     if (fname(i:i+3).eq.'.zin') then
        fname=fname(:i-1)
        exit
     endif
  enddo

  write (*,*) 'This is '//trim(name)
  print *,'Processing files beginning:'
  print *,trim(fname)

  date=' '
  time=' '
  call date_and_time(date,time)
  call cpu_time(t1)

! ----------------------------------------------------------------------
! Read the ZIN and CON files. Allocate arrays related to CAFQG, init, BC
! read files and fill these arrays. Set other vars on common.
! Set up procedure pointers if DLL to be linked (cfun...BCfun)
! also sets nvar_rst, nregel, nregnd
! ----------------------------------------------------------------------

  call readmatrix
  call usersetup(.false.)    ! no scanfun needed

! ----------------------------------------------------------------------
! We use 3 Gauss quadrature points in each active direction (might increase)
! Set to 1 if inactive
! ----------------------------------------------------------------------

  if (NDIM==3) then
     ng_xi=3
     ng_eta=3
     ng_mu=3
  else if (NDIM==2) then
     ng_xi=3
     ng_eta=3
     ng_mu=1
  else
     ng_xi=3
     ng_eta=1
     ng_mu=1
  endif

  print *,'ng=',ng_xi,ng_eta,ng_mu

  if (key_sim /=2 .and.(.not.isNL)) nstep=1   ! sanity check

! ----------------------------------------------------------------------
! Read IJKMAX from MTF file
! for 2D KMAX=0
! for 1D JMAX=KMAX=0
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.mtf',term,'old','ASCII')

  itot=0
  do
     call getline(1,ios,line)
     if (ios.ne.0) call prterr('Unexpected end of file')

     if (index(line,'IMAX').ne.0) then
        read (line,*,iostat=ios) cdum,imax
        itot=itot+1
     else if (index(line,'JMAX').ne.0) then
        read (line,*,iostat=ios) cdum,jmax
        itot=itot+1
     else if (index(line,'KMAX').ne.0) then
        read (line,*,iostat=ios) cdum,kmax
        itot=itot+1
     endif
     
     if (ios.ne.0) call prterr('Incomprehensible imax, jmax or kmax')
     if (itot.eq.3) exit
  enddo

! ----------------------------------------------------------------------
! Check NDIM versus MTF file
! ----------------------------------------------------------------------

  if (NDIM==2.and.kmax /= 0) then
     call prterr('Bad 2D MTF file, expect KMAX=0.')
  endif

  if (NDIM==1.and.(kmax /= 0.or.jmax /= 0)) then
     call prterr('Bad 1D MTF file, expect KMAX=0, JMAX=0.')
  endif

  if (NDIM==3.and.(kmax==0.or.jmax==0.or.imax==0)) then
     call prterr('Bad 3D MTF file, IMAX, JMAX, KMAX must be non-zero')
  endif

  if (NDIM==2.and.(jmax==0.or.imax==0)) then
     call prterr('Bad 2D MTF file, IMAX, JMAX must be non-zero')
  endif

! ----------------------------------------------------------------------
! Set nnod, ndof. 
! Allocate arrays vec,rnode,ireg,iregup,Qval,iQ,jQ,
! enreg,iusel,iregfix,regfix,
! irowst,irowed,RR,iunk
! for 2D k=0, for 1D j=k=0
! ----------------------------------------------------------------------

  nel=imax*jmax*kmax
  nnod=(imax+1)*(jmax+1)*(kmax+1)
  ndof=nnod*nvar

  allocerr=0

  allocate (irowst(ndof),irowed(ndof),stat=allocerr(1))
  allocate (rnode(0:imax,0:jmax,0:kmax,NDIM),stat=allocerr(2))
  allocate (ireg(0:imax,0:jmax,0:kmax),stat=allocerr(3))
  allocate (iregup(0:imax,0:jmax,0:kmax),stat=allocerr(4))

  if (key_Q.eq.1.or.key_Q.eq.2.or.key_Q.eq.4.or.key_Q.eq.5) then
     allocate (Qval(nvar**2*nnod*3**NDIM),iQ(nvar**2*nnod*3**NDIM),&
          jQ(nvar**2*nnod*3**NDIM),stat=allocerr(5))
  else if (key_Q.eq.3) then
     allocate(jQ(nvar**2*nnod*27+nvar**2*nel*8),&
          Qval(nvar**2*nnod*27+nvar**2*nel*8),&
          iQ(ndof),left(ndof),right(ndof),data(ndof),stat=allocerr(5))
  endif

  allocate (veclookup(ndof))
  allocate (enreg(nregel),stat=allocerr(6))
  allocate (vec(ndof),vecred(ndof),stat=allocerr(7))
  allocate (RR(ndof),iunk(ndof),iunkred(ndof),kfac(ndof),stat=allocerr(8))

  if (newton=='YES') then
     allocate (Qval0(size(Qval)),vecp(size(vec)),stat=allocerr(9))
     allocate (dvec(ndof),Nres(ndof),stat=allocerr(10))
  endif

  do i=1,10
     if (allocerr(i).ne.0) then
        write (temp,*) i
        call prterr('Failed to allocate array:'//trim(temp))
     endif
  enddo
  
! ----------------------------------------------------------------------
! Read geometry data. scale up
! ----------------------------------------------------------------------

  mtf_nregnd=0
  mtf_nregel=0

  do k=0,kmax
     do j=0,jmax
        do i=0,imax

           do
              call getline(1,ios,line)
              if (ios.ne.0) call prterr('Unexpected end of file')
              if (index(line,'REGNO').eq.0.and.index(line,'===').eq.0) exit
           enddo

           if (NDIM==3) then
              read (line,*,iostat=ios) ic,jc,kc,x,y,z,mreg,mregup
           else if (NDIM==2) then
              read (line,*,iostat=ios) ic,jc,x,y,mreg,mregup
           else
              read (line,*,iostat=ios) ic,x,mreg,mregup
           endif

           if (ios.ne.0) call prterr('Incomprehensible mesh file line')

           if (NDIM==3.and..not.(i.eq.ic.and.j.eq.jc.and.k.eq.kc)) then
              call prterr('Unexpected node: vary i first, then j, then k')
           endif

           if (NDIM==2.and..not.(i.eq.ic.and.j.eq.jc)) then
              call prterr('Unexpected node: vary i first, then j')
           endif

           if (NDIM==1.and..not.(i.eq.ic)) then
              call prterr('Unexpected node: vary i first')
           endif
           
           if ( .not.(mreg  .ge.1.and.mreg  .le.nregnd &
                .and. mregup.ge.0.and.mregup.le.nregel)) then
              print *,'MTF: ',i,j,k,mreg,mregup
              call prterr('Region number out of range in .mtf file:')
           endif
           
                        rnode(i,j,k,1)=x*scale
           if (NDIM>=2) rnode(i,j,k,2)=y*scale
           if (NDIM==3) rnode(i,j,k,3)=z*scale

           ireg(i,j,k)=mreg
           iregup(i,j,k)=mregup

           mtf_nregnd=max(mtf_nregnd,mreg)
           mtf_nregel=max(mtf_nregel,mregup)
        enddo
     enddo
  enddo

  print *,'Found ',mtf_nregel,' element regions'
  print *,'Found ',mtf_nregnd,' nodal regions'

  if (mtf_nregel /= nregel) call prterr('No. element regions does not match with ZIN file')
  if (mtf_nregnd /= nregnd) call prterr('No. nodal regions does not match with ZIN file')
  
  close (1)

! ----------------------------------------------------------------------
! If restart requested, read solution from restart file
! ----------------------------------------------------------------------

  if (restart /= 'NONE') then

     print *,'Loading restart file: ',trim(restart)

     allocate (uinit(0:imax,0:jmax,0:kmax,nvar_rst),stat=allocerr(1))
     if (allocerr(1)/=0) call prterr('Could not initialise uinit')

     term=1
     call initfile(1,restart,term,'old',rst_format)
  
     if (rst_format=="ASCII") then
     
        do
           call getline(1,ios,line)
           if (ios.ne.0) call prterr('Unexpected end of file, searching for "final"')
           if (index(line,'FINAL').ne.0) exit
        enddo

        ! allocate uinit, rs. defparam rs

        do i=0,imax
           do j=0,jmax
              do k=0,kmax
                 do ii=1,nvar_rst
                    
                    call getline(1,ios,line)
                    if (ios.ne.0) call prterr('Unexpected end of file')
                    
                    if (NDIM==3) then
                       read (line,*,iostat=ios) i1,j1,k1,ii1,u1
                       if (ios.ne.0) call prterr('Problem reading .rst file')
                    
                       if (.not.(i.eq.i1.and.j.eq.j1.and.k.eq.k1.and.ii.eq.ii1)) then
                          call prterr('Restart file out of sequence: vary k first, then j, then i.' &
                               //'Must be congruent with .MTF file')
                       endif
                    else if (NDIM==2) then
                       read (line,*,iostat=ios) i1,j1,ii1,u1
                       if (ios.ne.0) call prterr('Problem reading .rst file')
                    
                       if (.not.(i.eq.i1.and.j.eq.j1.and.ii.eq.ii1)) then
                          call prterr('Restart file out of sequence: vary k first, then j, then i.' &
                               //'Must be congruent with .MTF file')
                       endif
                    else
                       read (line,*,iostat=ios) i1,ii1,u1
                       if (ios.ne.0) call prterr('Problem reading .rst file')
                    
                       if (.not.(i.eq.i1.and.ii.eq.ii1)) then
                          call prterr('Restart file out of sequence: vary k first, then j, then i.' &
                               //'Must be congruent with .MTF file')
                       endif
                    endif

                    uinit(i,j,k,ii)=u1
                 
                 enddo
              enddo
           enddo
        enddo
        
     else               ! Binary
        
        read (1)        ! Ignore strapline
        do
           read (1,iostat=ios) i
           if (ios /= 0) call prterr('Unexpected end of file searching for final snapshot')
           if (i == -2) exit           ! i j k ii val: i here will not be negative
        enddo
        
        do i=0,imax
           do j=0,jmax
              do k=0,kmax
                 do ii=1,nvar_rst
                    if (NDIM==3) then
                       read (1,iostat=ios) i1,j1,k1,ii1,u1
                       if (ios.ne.0) call prterr('Problem reading .rst file')
                       
                       if (.not.(i.eq.i1.and.j.eq.j1.and.k.eq.k1.and.ii.eq.ii1)) then
                          call prterr('Restart file out of sequence: vary k first, then j, then i')
                       endif
                    else if (NDIM==2) then
                       read (1,iostat=ios) i1,j1,ii1,u1
                       if (ios.ne.0) call prterr('Problem reading .rst file')
                       
                       if (.not.(i.eq.i1.and.j.eq.j1.and.ii.eq.ii1)) then
                          call prterr('Restart file out of sequence: vary k first, then j, then i')
                       endif
                    else
                       read (1,iostat=ios) i1,ii1,u1
                       if (ios.ne.0) call prterr('Problem reading .rst file')
                       
                       if (.not.(i.eq.i1.and.ii.eq.ii1)) then
                          call prterr('Restart file out of sequence: vary k first, then j, then i')
                       endif
                    endif

                    uinit(i,j,k,ii)=u1
                 enddo
              enddo
           enddo
        enddo
        
     endif

     close (1)

  endif

! ----------------------------------------------------------------------
! Set initial values either by interpolation or as specified
! by user (evaluate expression if needed)
! EVERY var on EVERY node is set a value here
! ----------------------------------------------------------------------

  do ii=1,nvar
     if (inits(ii).eq.'INTERP') then
        do i=0,imax
           do j=0,jmax
              do k=0,kmax
                 ireg1=ireg(i,j,k)
                 if (iregfix(ireg1,ii).eq.1) then
                    vec(indQ(i,j,k,ii))=regfix(ireg1,ii)
                 else if (iregfix(ireg1,ii).gt.1) then
                    call prterr('Fixed nodes must be set constant if INTERP is to be used')
                 endif
              enddo
           enddo
        enddo
        
        call interp(ii)          ! sets ii var only for non-fixed nodes
                                 ! in terms of ii fixed values
     else
        
        do i=0,imax
           do j=0,jmax
              do k=0,kmax
                              x=rnode(i,j,k,1)
                 if (NDIM>=2) y=rnode(i,j,k,2)
                 if (NDIM==3) z=rnode(i,j,k,3)
                 
                 ireg1=ireg(i,j,k)
                 
                 if (iregfix(ireg1,ii) == -1) then           ! init per region
                    vec(indQ(i,j,k,ii))=regfix(ireg1,ii)
                 else if (iregfix(ireg1,ii) == -2) then      ! init per region (x,y,z)
                                 call defparam('X',x)
                    if (NDIM>=2) call defparam('Y',y)
                    if (NDIM==3) call defparam('Z',z)
                    
                    call evaltoken(tregfix(ireg1,ii),val)
                    
                    if (ierr.ne.0) then
                       print *,'ERROR: Could not resolve initial expression for variable ',ii
                       call exit(1)
                    endif
                    
                    vec(indQ(i,j,k,ii))=val
                    
                 else if (initi(ii)==0) then               ! global init
                    vec(indQ(i,j,k,ii))=init(ii)
                 else if (initi(ii)==1) then               ! global init (x,y,z)
                                 call defparam('X',x)
                    if (NDIM>=2) call defparam('Y',y)
                    if (NDIM==3) call defparam('Z',z)
                    
                    call evaltoken(initt(ii),val)
                    
                    if (ierr.ne.0) then
                       print *,'ERROR: Could not resolve initial expression for variable ',ii
                       call exit(1)
                    endif
                    
                    vec(indQ(i,j,k,ii))=val
                 else if (initi(ii)==2) then              ! global init eg R4*2+x
                                 call defparam('X',x)
                    if (NDIM>=2) call defparam('Y',y)
                    if (NDIM==3) call defparam('Z',z)
                    
                    do kk=1,nvar_rst
                       call defparam(rs(kk),uinit(i,j,k,kk))
                    enddo

                    call evaltoken(initt(ii),val)

                    if (ierr.ne.0) then
                       print *,'ERROR: Could not resolve initial expression for variable ',ii
                       call exit(1)
                    endif
                    
                    vec(indQ(i,j,k,ii))=val
                 else
                    call prterr('Illegal option during initialisation')
                 endif
                 
              enddo
           enddo
        enddo
     endif
  enddo

  if (allocated(uinit)) deallocate(uinit)

! ----------------------------------------------------------------------
! set fixed nodes either linear or NL
! ----------------------------------------------------------------------

  call fixnodes

! ----------------------------------------------------------------------
! prepare iunk cf vec showing which components of nodes are unknown
! iregfix(1..nreg,ii)=1 if u_ii fixed. regfix(1..nreg,ii) value of fix
! Also prepare nunk, total number of unknowns
! ----------------------------------------------------------------------

  nunk=0
  ind=0
  do k=0,kmax
     do j=0,jmax
        do i=0,imax
           do ii=1,nvar
              ind=ind+1
              ialpha=indQ(i,j,k,ii)
              
              if (ialpha.ne.ind) then
                 write (*,*) 'ERROR: wrong index'
                 call exit(1)
              endif
              
              ireg1=ireg(i,j,k)
              if (iregfix(ireg1,ii).le.0) then   ! <0 means init per region
                 nunk=nunk+1
                 iunk(ialpha)=1
              else
                 iunk(ialpha)=0
              endif
           enddo
        enddo
     enddo
  enddo
  
! ----------------------------------------------------------------------
! Find cuboid which contains entire simulation. store in redge(3,2)
! ----------------------------------------------------------------------

  ifirst=1
  do i=0,imax
     do j=0,jmax
        do k=0,kmax
           
           do ic=1,NDIM
              r(ic)=rnode(i,j,k,ic)
           enddo
           
           if (ifirst.eq.1) then
              do ic=1,NDIM
                 redge(ic,1)=r(ic)
                 redge(ic,2)=r(ic)
              enddo
              ifirst=0
           else
              do ic=1,NDIM
                 redge(ic,1)=min(r(ic),redge(ic,1))
                 redge(ic,2)=max(r(ic),redge(ic,2))
              enddo
           endif
        enddo
     enddo
  enddo
  
! ----------------------------------------------------------------------
! Open output file .zls and write basic geometry details
! Write out unknowns if key_db=1
! ----------------------------------------------------------------------

  currfile=trim(fname)//'.zls'

  term=1
  call initfile(1,currfile,term,'unknown','ASCII')

  if (NDIM==3) then
     write (1,*) 'imax,jmax,kmax=',imax,jmax,kmax
     write (1,*) '(Note that i=0,1,...,imax etc)'
     
     write (1,*) 'xrange:',redge(1,1),redge(1,2)
     write (1,*) 'yrange:',redge(2,1),redge(2,2)
     write (1,*) 'zrange:',redge(3,1),redge(3,2)
  else if (NDIM==2) then
     write (1,*) 'imax,jmax=',imax,jmax
     write (1,*) '(Note that i=0,1,...,imax etc)'
     
     write (1,*) 'xrange:',redge(1,1),redge(1,2)
     write (1,*) 'yrange:',redge(2,1),redge(2,2)
  else
     write (1,*) 'imax=',imax
     write (1,*) '(Note that i=0,1,...,imax etc)'
     
     write (1,*) 'xrange:',redge(1,1),redge(1,2)
  endif

  if (key_db.eq.1) then
     write (1,*) 'Sample of iunk vector'
     write (1,'(6a6)') 'alpha','i','j','k','ii','iunk'

     do k=0,min(kmax,1)
        do j=0,jmax
           do i=0,imax
              do ii=1,nvar
                 ialpha=indQ(i,j,k,ii)
                 write (1,'(6i6)') ialpha,i,j,k,ii,iunk(ialpha)
              enddo
           enddo
        enddo
     enddo
  endif

! ----------------------------------------------------------------------
! Write out matrix details to unit 1
! ----------------------------------------------------------------------

  call wrtmatrix(1)

! ----------------------------------------------------------------------
! *** This section is for transient problems, key_sim=2 ***
! ----------------------------------------------------------------------

  if (key_sim == 2) then

! ----------------------------------------------------------------------
! Form the Q matrix and R and kfac vectors
! In non-linear case, these will depend on initial state from vec
! symmetrise and archive the Q matrix and write debug info
! ----------------------------------------------------------------------
     
     call getk                 ! only need to do this once
     
     if (key_Q == 1) then
        call getQ(.false.)          ! do not remove fixed
     else if (key_Q == 2) then
        call getQ2(.false.)
     else if (key_Q == 3) then
        call getQ3                  ! deprecated
     else if (key_Q == 4) then
        call getQcuboid(.false.)
     else if (key_Q == 5) then
        call getQcuboid2(.false.)
     else
        call prterr('Illegal key_Q')
     endif
  
     if (key_Q /=3) then
        if (.not.sorted(iQ,lenQ)) call heapsort(iQ,lenQ)    ! heapsort calls swap
        call archive
     endif

     if (key_sym == 1) call symmetrise
  
     if (key_Q==4.or.key_Q==5) then
        call getRcuboid(.false.)
     else
        call getR(.false.)   ! starts from scratch
     endif

! ----------------------------------------------------------------------
! Write info about Q to *.ZLS. Write out some of the Q matrix
! in debugging case
! ----------------------------------------------------------------------
  
     call Qinfo(.false.)   ! not closed up

! ----------------------------------------------------------------------
! Check that empty rows correspond to fixed nodes
! ----------------------------------------------------------------------

     if (key_Q /= 3) then
        do i=1,ndof
           if (irowst(i) == 1.and.irowed(i) == 0.and.iunk(i) == 1) then
              write (*,*) 'Variable ',i,' has empty row but not fixed'
              call exit(1)
           endif
        enddo
     endif
     
! ----------------------------------------------------------------------
! Write out initial energy calculated without Q, R from func energy
! ----------------------------------------------------------------------
     
     en=energy(enreg)
     write (*,*) 'Initl energy       =',en
     write (1,*) 'Initl energy       =',en

! ----------------------------------------------------------------------
! Calculate energy from Q, R and vec
! ----------------------------------------------------------------------

     enQ=0
     do ip=1,lenQ
        iQ1=iQ(ip)
        jQ1=jQ(ip)
        
        enQ=enQ+0.5*vec(iQ1)*Qval(ip)*vec(jQ1)
     enddo
  
     do i=1,ndof
        enQ=enQ-vec(i)*RR(i)
     enddo
  
     write (1,*) 'Initl energy from Q=',enQ
     write (*,*) 'Initl energy from Q=',enQ

! ----------------------------------------------------------------------
! open .zou and write header. May be in Binary of ASCII format
! ----------------------------------------------------------------------

     term=1
     currfile=trim(fname)//'.zou'
     call initfile(2,currfile,term,'unknown',zou_format)

! **********************************************************************
! ----------------------------------------------------------------------
! Main loop for transient, iterate vec using Q, R, kfac matrices
! Q is matrix, iQ(ip),jQ(ip),Qval(ip), ip=1..lenQ
! Q_ij * vec_j = R_i for i is in unknown set
! ----------------------------------------------------------------------
! **********************************************************************

     write (*,'(a7,4a14)')   'istep','|resid|/|lhs|','|lhs|','|rhs|','|diff|'
     write (1,'(//a7,4a14)') 'istep','|resid|/|lhs|','|lhs|','|rhs|','|diff|'
     
     tind=0
     do istep=1,nstep

! ----------------------------------------------------------------------
! Perform one update step. Write convergence every nstride steps
! ----------------------------------------------------------------------

        call update_transient(mod(istep,nstride)==0)  ! true means write convergence

! ----------------------------------------------------------------------
! In transient case, write snapshots
! write in paraview format if requested
! ----------------------------------------------------------------------

        if (mod(istep,nstride)==0) then

           if (zou_format=="ASCII") then
              write (2,*)
              write (2,'(a,i10,a,e13.5)') 'istep=',istep,' t=',istep*tstep
           else
              write (2) -1,istep,istep*tstep
           endif

           call snapshot(2)

           if (key_export == 1) then
              tind=tind+1
              call pvsnapshot(10,tind)
           endif
        endif

! ----------------------------------------------------------------------
! update Q, R matrices if necessary (non-linear case)
! Could write more sophisticated routines in future, same prog structure
! ----------------------------------------------------------------------

        if (mod(istep,nstrideNL) == 0.and.istep.ne.nstep) then
           if (isQNL) then
              
              if (key_Q == 1) then
                 call getQ(.false.)
              else if (key_Q == 2) then
                 call getQ2(.false.)
              else if (key_Q == 3) then
                 call getQ3
              else if (key_Q == 4) then
                 call getQcuboid(.false.)
              else if (key_Q == 5) then
                 call getQcuboid2(.false.)
              endif

              if (key_Q /= 3) then
                 if (.not.sorted(iQ,lenQ)) call heapsort(iQ,lenQ)
                 call archive
              endif

              if (key_sym == 1) call symmetrise
           endif
        
           if (isRNL) then
              if (key_Q==4.or.key_Q==5) then
                 call getRcuboid(.false.)
              else
                 call getR(.false.)
              endif
           endif

           if (isfixNL) call fixnodes

        endif

     enddo
     
! ----------------------------------------------------------------------
! calculate u from final vec and hence get final energy from 
! energy func (which does not use Q)
! ----------------------------------------------------------------------

     en=energy(enreg)
        
     write (1,*) 'Final Energy       =',en
     write (*,*) 'Final Energy       =',en
     
     write (1,'(a8,a13)') 'Region','Energy'
     do i=1,nregel
        write (1,'(i8,e13.5)') i,enreg(i)
     enddo
  
! ----------------------------------------------------------------------
! Calculate final energy from Q. Write out and close .zls
! ----------------------------------------------------------------------

     enQ=0
     do ip=1,lenQ
        iQ1=iQ(ip)
        jQ1=jQ(ip)
        
        enQ=enQ+0.5*vec(iQ1)*Qval(ip)*vec(jQ1)
     enddo
  
     do i=1,ndof
        enQ=enQ-vec(i)*RR(i)
     enddo
  
     write (1,*) 'Final energy from Q=',enQ
     write (*,*) 'Final energy from Q=',enQ

     call cpu_time(t2)
     write (1,*) 'Time elaspsed',t2-t1,' seconds'
     close (1)            ! .zls file

! ----------------------------------------------------------------------
! Write final state to .ZOU solution file. Close file (2).
! Also write residual file and paraview
! ----------------------------------------------------------------------

     write (*,*) 'Writing final state to .zou'
  
     if (zou_format=="ASCII") then
        write (2,*)
        write (2,*) '(final) istep=',nstep,' t=',nstep*tstep
     else
        write (2) -2,nstep,nstep*tstep
     endif

     if (nodecheck=="YES") then
        call nodecheck_snapshot(2)
     else
        call snapshot(2)
     endif

     if (residual/="NONE") then
     
        term=1
        call initfile(3,trim(fname)//'.resid',term,'unknown',zou_format)

        if (zou_format=="ASCII") then
           write (3,*)
           write (3,*) '(final) istep=',nstep,' t=',nstep*tstep
        else
           write (3) -2,nstep,nstep*tstep
        endif

        if (residual=='UNKNOWNS') then
           call resid_snapshot(3)
        else
           call resid_snapshot2(3)
        endif

        close (3)
     endif

     if (key_export == 1) then
        if (key_sim == 1.or.key_sim == 0) then
           call pvsnapshot(10,0)
        else
           tind=tind+1
           call pvsnapshot(10,tind)
        endif
     endif

     close (2)

! ----------------------------------------------------------------------
! OR, use SOR, SLAP, UMFPACK to solve the equations and iterate only in NL situations
! ----------------------------------------------------------------------
  
  else

     if (removeFixed=='YES') then
        ndofx=ndofred
     else
        ndofx=ndof
     endif

! ----------------------------------------------------------------------
! Do either Newton or regular non-linear run (includes linear with nstep=1)
! Do not close up for Newton
! ----------------------------------------------------------------------

     if (newton=='YES') then

        do istep=1,nstep

           if (key_Q == 1) then
              call getQ(.false.)
           else if (key_Q == 2) then
              call getQ2(.false.)
           else if (key_Q==4) then
              call getQcuboid(.false.)
           else if (key_Q==5) then
              call getQcuboid2(.false.)              
           else
              call prterr('Illegal key_Q')
           endif
              
           if (.not.sorted(iQ,lenQ)) call heapsort(iQ,lenQ)    ! heapsort calls swap
           call archive
           if (key_sym == 1) call symmetrise

           if (key_Q==4.or.key_Q==5) then
              call getRcuboid(.false.)
           else
              call getR(.false.)
           endif

           if (istep==1) call Qinfo(.false.)
              
           call fixQ
           call fixR

           resid=0; residDiag=0
           lhs=0; lhsDiag=0
           rhs=0; rhsDiag=0
           
           do i=1,ndof
              lhsI=0
              rhsI=RR(i)

              pivot=0
              Nres(i)=RR(i)
              do ip=irowst(i),irowed(i)
                 jQ1=jQ(ip)
                 Nres(i)=Nres(i)-Qval(ip)*vec(jQ(ip))

                 if (jQ1==i) pivot=Qval(ip)
                 lhsI=lhsI+Qval(ip)*vec(jQ1)
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

           resid=sqrt(resid)/sqrt(lhs)
           lhs=sqrt(lhs)
           rhs=sqrt(rhs)

           residDiag=sqrt(residDiag)/sqrt(lhsDiag)
           lhsDiag=sqrt(lhsDiag)
           rhsDiag=sqrt(rhsDiag)

           print *,'RESIDUAL=',istep,sqrt(sum(Nres(:ndof)**2))

           print '(10a13)','NL step','iterations','sol. error','error index','resid','lhs','rhs','residDiag','lhsDiag','rhsDiag'
           print '(a52,6e13.5)','OVERALL RESIDUAL',resid,lhs,rhs,residDiag,lhsDiag,rhsDiag

           write (1,'(10a13)') 'NL step','iterations','sol. error','error index','resid','lhs','rhs','residDiag','lhsDiag','rhsDiag'
           write (1,'(a52,6e13.5)') 'OVERALL RESIDUAL',resid,lhs,rhs,residDiag,lhsDiag,rhsDiag

           resid0=resid
           Qval0=Qval           ! need to store Q matrix to backtracking

           if (key_Q==4.or.key_Q==5) then
              call getJacobianCuboid
           else
              call getJacobian        ! Q -> J
           endif

           if (istep==1) then
              write (1,*) 'The following relates to the Jacobian:'
              call Qinfo(.false.)
           endif

           dvec=0                    ! initial guess
        
           ! solve J*dvec = Nres
!             call sor(ndof,Nres,dvec,lenQ,iQ,jQ,Qval,tol,itmax,iter,err,ierrslap,0, irowst,irowed,omega,nstride)

           print *,'Starting solver'
           call solve(ndof,dvec,Nres,itmax,err,ierrslap,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag)

           print '(10a13)','NL step','iterations','sol. error','error index','resid','lhs','rhs','residDiag','lhsDiag','rhsDiag'
           print '(2i13,e13.5,i13,6e13.5)',istep,iter,err,ierrslap,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag
           write (1,'(2i13,e13.5,i13,6e13.5)') istep,iter,err,ierrslap,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag

! ----------------------------------------------------------------------
! Do backtracking line search. The new solution must actually
! reduce the L2-norm residual of Qu=R. Displacement might be lower than dvec
! ----------------------------------------------------------------------

           if (backtrack) then
              s=1.0
              found=.false.

              do iback=1,32
                 vecp(:ndof)=vec(:ndof)+s*dvec(:ndof)
                 
                 do i=1,ndof
                    Nres(i)=RR(i)
                    do ip=irowst(i),irowed(i)
                       Nres(i)=Nres(i)-Qval0(ip)*vecp(jQ(ip))
                    enddo
                 enddo
                 
                 resid=sqrt(sum(Nres(:ndof)**2))
                 
                 print *,'backtracking: ',iback,resid
                 
                 if (resid < (1-0.01*s)*resid0) then
                    found=.true.
                    exit
                 endif
                 
                 s=s/2
              enddo
              
              if (.not.found) then
                 print *,'Stuck backtracking, will return best solution'
                 exit
              endif
              
              vec(:ndof)=vecp(:ndof)

           else
              vec(:ndof)=vec(:ndof)+dvec(:ndof)
           endif
              
        enddo
     else

! ----------------------------------------------------------------------
! Simple non-linear (with possible removal of fixed dofs)
! ----------------------------------------------------------------------

        do istep=1,nstep

           if (istep > 1.and.isfixNL) call fixnodes

! ----------------------------------------------------------------------
! (re)calculate Q,R. At first step, output energies for checking
! ----------------------------------------------------------------------
                 
           if (istep == 1.or.isQNL.or.isRNL) then
              if (key_Q == 1) then
                 call getQ(removeFixed=='YES')            ! based on new "vec"
              else if (key_Q == 2) then
                 call getQ2(removeFixed=='YES')
              else if (key_Q == 4) then
                 call getQcuboid(removeFixed=='YES')
              else if (key_Q == 5) then
                 call getQcuboid2(removeFixed=='YES')
              else
                 call prterr('Illegal key_Q')
              endif
                       
              if (.not.sorted(iQ,lenQ)) call heapsort(iQ,lenQ)
              call archive
              if (key_sym == 1) call symmetrise
              
              if (key_Q==4.or.key_Q==5) then
                 call getRcuboid(removeFixed=='YES')
              else
                 call getR(removeFixed=='YES')
              endif

! energy checks
              if (istep == 1) then

                 en=energy(enreg)
                 write (*,*) 'Initl energy       =',en  ! correct energy
                 write (1,*) 'Initl energy       =',en

                 if (removeFixed=="NO") then
                    enQ=0
                    do ip=1,lenQ
                       iQ1=iQ(ip)
                       jQ1=jQ(ip)
                       
                       enQ=enQ+0.5*vec(iQ1)*Qval(ip)*vec(jQ1)
                    enddo
                    
                    do i=1,ndof
                       enQ=enQ-vec(i)*RR(i)
                    enddo
                 
                    write (1,*) 'Initl energy from Q=',enQ      ! not correct for removeFixed=='YES'
                    write (*,*) 'Initl energy from Q=',enQ
                 endif

                 if (removeFixed=='YES') then
                    call closeup             ! creates reduced vector vecred
                    call expand              ! round trip. vec should be same (vecred still available)

                    en=energy(enreg)
                    write (*,*) 'Initl energy (round trip)  =',en
                    write (1,*) 'Initl energy (round trip)  =',en
                 endif

              endif
! end of checks

              if (removeFixed=='NO') then
                 call fixQ            ! for fixed nodes, equation reads u=value
                 call fixR
              endif

              if (istep==1) call Qinfo(removeFixed=='YES')

              do i=1,ndofx
                 if (isnan(RR(i))) then
                    print *,'ERROR: NaN in RR(',i,')'
                    stop
                 endif
              enddo
              
              do ip=1,lenQ
                 if (isnan(Qval(ip))) then
                    print *,'ERROR: NaN in Q: ',iQ(ip),jQ(ip)
                    stop
                 endif
              enddo
              
           endif

! ----------------------------------------------------------------------
! solve equations, output, loop
!    subroutine solve(ndofx,vecx,rhsx,iter,err,ierrslap)
! ----------------------------------------------------------------------
                 
           if (removeFixed=='YES') then
              if (istep>1) call closeup
              call solve(ndofred,vecred,RR,iter,err,ierrslap,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag)
              call expand           ! new vec
           else
              call solve(ndof,vec,RR,iter,err,ierrslap,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag)
           endif

           print '(10a13)','NL step','iterations','sol. error','error index','resid','lhs','rhs','residDiag','lhsDiag','rhsDiag'
           print '(2i13,e13.5,i13,6e13.5)',istep,iter,err,ierrslap,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag

           write (1,'(10a13)') 'NL step','iterations','sol. error','error index','resid','lhs','rhs','residDiag','lhsDiag','rhsDiag'
           write (1,'(2i13,e13.5,i13,6e13.5)') istep,iter,err,ierrslap,resid,lhs,rhs,residDiag,lhsDiag,rhsDiag

           if (ierrslap == 2) then
              print *,    'WARNING: Solver may not have converged fully. Treat solution with caution'
              write (1,*) 'WARNING: Solver may not have converged fully. Treat solution with caution'
           endif

        enddo
     endif             ! newton or not

! ----------------------------------------------------------------------
! calculate final energy from energy func (which does not use Q)
! ----------------------------------------------------------------------

     en=energy(enreg)
     
     write (1,*) 'Final Energy       =',en
     write (*,*) 'Final Energy       =',en
     
     write (1,'(a8,a13)') 'Region','Energy'
     do i=1,nregel
        write (1,'(i8,e13.5)') i,enreg(i)
     enddo

     if (.false.) then   ! energy will never be correct
        enQ=0
        do ip=1,lenQ
           iQ1=iQ(ip)
           jQ1=jQ(ip)
           
           enQ=enQ+0.5*vec(iQ1)*Qval(ip)*vec(jQ1)
        enddo
        
        do i=1,ndof
           enQ=enQ-vec(i)*RR(i)
        enddo
        
        write (1,*) 'Final energy from Q=',enQ
        write (*,*) 'Final energy from Q=',enQ
     endif

! ----------------------------------------------------------------------
! Write final snapshot and export to paraview if needed
! ----------------------------------------------------------------------

     call cpu_time(t2)
     write (1,*) 'Time elaspsed',t2-t1,' seconds'
     close (1)   ! .zls file

     write (*,*) 'Writing final state to .zou'

     term=1
     currfile=trim(fname)//'.zou'
     call initfile(2,currfile,term,'unknown',zou_format)
       
     if (zou_format=="ASCII") then
        write (2,*)
        write (2,*) '(final) istep=',nstep
     else
        write (2) -2,nstep,nstep*tstep
     endif

     if (nodecheck=="YES") then
        call nodecheck_snapshot(2)
     else
        call snapshot(2)
     endif

     if (residual /= "NONE") then

        term=1
        call initfile(3,trim(fname)//'.resid',term,'unknown',zou_format)

        if (zou_format=="ASCII") then
           write (3,*)
           write (3,*) '(final) istep=',nstep,' t=',nstep*tstep
        else
           write (3) -2,nstep,nstep*tstep
        endif
        
        if (residual=='UNKNOWNS') then
           call resid_snapshot(3)
        else
           call resid_snapshot2(3)
        endif

        close (3)
     endif

     if (key_export == 1) then
        call pvsnapshot(10,0)
     endif

     close (2)

     call freesolvers
  endif     ! transient or static

! ----------------------------------------------------------------------
! Deallocate all arrays
! ----------------------------------------------------------------------

  call freematrix
  call userfree

  deallocate (rnode,ireg,iregup,Qval,iQ,jQ)
  deallocate (enreg,vec,vecred,veclookup,RR,iunk,iunkred,kfac)
  deallocate (irowst,irowed,us,dus,rs)
  if (newton=='YES') deallocate(Qval0,vecp,dvec,Nres)

  print *,'Zinc completed successfully'
  print *,'Simulation: ',trim(fname)

end program zinc
