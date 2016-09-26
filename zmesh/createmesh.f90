module createmesh

! ----------------------------------------------------------------------
! call processMIN to process MIN file.
! uses global and allocates rnode etc using alloc (from global)
! sets rnode etc in global
! Writes .MLS and .MTF
! User must supply external createerror(s), createmessage(s) for error reporting
!
! ifort error: local integer silently supersedes global [not an error]
! ----------------------------------------------------------------------

  use snap

  implicit none

  save
  private

  public :: processMIN

  integer, allocatable :: fixed(:,:,:)

  real, allocatable :: xspring(:),yspring(:),zspring(:)
  real, allocatable :: spring000(:,:,:)
  real, allocatable :: spring100(:,:,:)
  real, allocatable :: spring110(:,:,:)
  real, allocatable :: spring010(:,:,:)
  real, allocatable :: force(:,:,:,:)

  real, allocatable :: xmeshline(:,:),ymeshline(:,:),zmeshline(:,:),xmesh(:),ymesh(:),zmesh(:)

  type(Tpart), allocatable :: parts(:)
  integer nparts

  integer :: nstep       ! relaxation steps (also no. of fixinv steps)
  integer, parameter :: tstep=1        ! delta-t
  real, parameter :: konst=0.1         ! only tstep*konst is important: should be about 0.1
  logical :: fix    ! Try to fix inverted elements?
  logical :: lrelax ! Relax after each surface fitting?

contains

  function processMIN()

    use global

    integer processMIN
    character(1000) line
    character(5) line5
    character(20) dum
    real x1,x2,dx,y1,y2,dy,z1,z2,dz
    integer nx,ny,nz

    integer nxmeshline,nymeshline,nzmeshline,nmeshline,i,j,k,ios
    integer allocerr,allocerrs(10)

! ----------------------------------------------------------------------
! Read global part of min file
! ----------------------------------------------------------------------

    processMIN=1      ! error. Early return == error

    call createmessage('Processing mesh...')
    
    open (1,file=MinName,status='old',iostat=ios)

    if (ios.ne.0) then
       call createerror('Could not open:'//trim(MinName))
       close (1)
       return
    endif

    nstep=10               ! default
    fix=.true.
    lrelax=.true.

    nxmeshline=0
    nymeshline=0
    nzmeshline=0
    do
       read (1,'(a)',iostat=ios) line
       
       if (ios.ne.0) then
          call createerror('Unexpected end of file in GLOBAL section')
          close (1)
          return
       endif

       call tolower(line)
       line=adjustl(line)
       line5=line(1:5)
       
       if (line(1:1).eq.'*'.or.len_trim(line).eq.0) cycle
       
       if (line(1:6).eq.'global') then
                
       else if (line(1:6).eq.'format') then

       else if (line(1:6).eq.'smooth') then
          read (line,*,iostat=ios) dum,nstep
          if (ios.ne.0) then
             call createerror('Incomprehensible smooth command')
             close (1)
             return
          endif

       else if (line(1:3).eq.'fix') then
          read (line,*,iostat=ios) dum,fix
          if (ios.ne.0) then
             call createerror('Incomprehensible fix command')
             close (1)
             return
          endif

       else if (line(1:5).eq.'relax') then
          read (line,*,iostat=ios) dum,lrelax
          if (ios.ne.0) then
             call createerror('Incomprehensible relax command')
             close (1)
             return
          endif          

       else if (line(1:3).eq.'end') then
          exit
       else if (line5.eq.'xmesh'.or.line5.eq.'ymesh'.or.line5.eq.'zmesh') then
          nmeshline=0
          do
             read (1,'(a)',iostat=ios) line
             
             if (ios.ne.0) then
                call createerror('Unexpected end of file reading '//line5//' section')
                close (1)
                return
             endif
             
             call tolower(line)
             line=adjustl(line)
             if (line(1:3).eq.'end') exit
             
             nmeshline=nmeshline+1
          enddo
                
          do i=1,nmeshline+1
             backspace(1)
          enddo
          
          if (line5.eq.'xmesh') then
             nxmeshline=nmeshline

             if (allocated(xmeshline)) deallocate(xmeshline)
             allocate (xmeshline(nxmeshline,3),stat=allocerr)
             if (allocerr.ne.0) then
                call createerror('Could not allocate xmeshline')
                close (1)
                return
             endif
             
             do i=1,nmeshline
                read (1,*,iostat=ios) (xmeshline(i,j),j=1,3)
                
                if (ios.ne.0) then
                   call createerror('Incomprehensible number in xmesh section')
                   close(1)
                   return
                endif
             enddo
             read (1,*)
          else if (line5.eq.'ymesh') then
             nymeshline=nmeshline

             if (allocated(ymeshline)) deallocate(ymeshline)
             allocate (ymeshline(nymeshline,3),stat=allocerr)
             if (allocerr.ne.0) then
                call createerror('Could not allocate ymeshline')
                close (1)
                return
             endif

             do i=1,nmeshline
                read (1,*,iostat=ios) (ymeshline(i,j),j=1,3)
                
                if (ios.ne.0) then
                   call createerror('Incomprehensible number in ymesh section')
                   close (1)
                   return
                endif
                
             enddo
             read (1,*)
          else
             nzmeshline=nmeshline

             if (allocated(zmeshline)) deallocate(zmeshline)
             allocate (zmeshline(nzmeshline,3),stat=allocerr)
             if (allocerr.ne.0) then
                call createerror('Could not allocate zmeshline')
                close (1)
                return
             endif
             
             do i=1,nmeshline
                read (1,*,iostat=ios) (zmeshline(i,j),j=1,3)
                
                if (ios.ne.0) then
                   call createerror('Incomprehensible number in zmesh section')
                   close (1)
                   return
                endif
                
             enddo
             read (1,*)
          endif
       endif
       
    enddo

    if (nxmeshline.eq.0) then
       call createerror('Did not find xmesh block in global')
       close (1)
       return
    endif
    
    if (nymeshline.eq.0) then
       call createerror('Did not find ymesh block in global')
       close (1)
       return
    endif
    
    if (nzmeshline.eq.0) then
       call createerror('Did not find zmesh block in global')
       close (1)
       return
    endif
    
! ----------------------------------------------------------------------
! Set xmesh
! ----------------------------------------------------------------------

    imax=0
    do i=1,nxmeshline
       x1=xmeshline(i,1)
       x2=xmeshline(i,2)
       dx=xmeshline(i,3)
       
       nx=nint((x2-x1)/dx)
       
       imax=imax+nx
    enddo
    
    if (allocated(xmesh)) deallocate(xmesh)
    allocate(xmesh(0:imax),stat=allocerr)
    if (allocerr.ne.0) then
       call createerror('Could not allocate xmesh')
       close (1)
       return
    endif

    call lmesh(xmeshline,xmesh)
    
! ----------------------------------------------------------------------
! Set ymesh
! ----------------------------------------------------------------------

    jmax=0
    do i=1,nymeshline
       y1=ymeshline(i,1)
       y2=ymeshline(i,2)
       dy=ymeshline(i,3)
       ny=nint((y2-y1)/dy)
       
       jmax=jmax+ny
    enddo

    if (allocated(ymesh)) deallocate(ymesh)
    allocate(ymesh(0:jmax),stat=allocerr)
    if (allocerr.ne.0) then
       call createerror('Could not allocate ymesh')
       close (1)
       return
    endif    

    call lmesh(ymeshline,ymesh)

! ----------------------------------------------------------------------
! Set zmesh
! ----------------------------------------------------------------------

    kmax=0
    do i=1,nzmeshline
       z1=zmeshline(i,1)
       z2=zmeshline(i,2)
       dz=zmeshline(i,3)
       nz=nint((z2-z1)/dz)
       
       kmax=kmax+nz
    enddo
    
    if (allocated(zmesh)) deallocate(zmesh)
    allocate(zmesh(0:kmax),stat=allocerr)
    if (allocerr.ne.0) then
       call createerror('Could not allocate zmesh')
       close (1)
       return
    endif
    
    call lmesh(zmeshline,zmesh)

! ----------------------------------------------------------------------
! allocate storage based on ijkmax
! free and alloc are in global
! ----------------------------------------------------------------------

    if (allocated(rnode)) call free
    call alloc

    if (allocated(fixed)) deallocate(fixed)

    if (allocated(xspring)) deallocate(xspring)
    if (allocated(yspring)) deallocate(yspring)
    if (allocated(zspring)) deallocate(zspring)

    if (allocated(spring000)) deallocate(spring000)
    if (allocated(spring100)) deallocate(spring100)
    if (allocated(spring110)) deallocate(spring110)
    if (allocated(spring010)) deallocate(spring010)

    if (allocated(force)) deallocate(force)

    allocate (fixed(0:imax,0:jmax,0:kmax),stat=allocerrs(1))

    allocate (xspring(0:imax),stat=allocerrs(2))
    allocate (yspring(0:jmax),stat=allocerrs(3))
    allocate (zspring(0:kmax),stat=allocerrs(4))

    allocate (spring000(0:imax,0:jmax,0:kmax),stat=allocerrs(5))
    allocate (spring100(0:imax,0:jmax,0:kmax),stat=allocerrs(6))
    allocate (spring110(0:imax,0:jmax,0:kmax),stat=allocerrs(7))
    allocate (spring010(0:imax,0:jmax,0:kmax),stat=allocerrs(8))

    allocate (force(0:imax,0:jmax,0:kmax,3),stat=allocerrs(9))

    do i=1,9
       if (allocerrs(i).ne.0) then
          call createerror('Could not allocate spring, force matrices')
          close (1)
          return
       endif
    enddo

! ----------------------------------------------------------------------
! Fill logical mesh
! ----------------------------------------------------------------------

    do i=0,imax
       do j=0,jmax
          do k=0,kmax
             rnode(i,j,k,1)=xmesh(i)
             rnode(i,j,k,2)=ymesh(j)
             rnode(i,j,k,3)=zmesh(k)
          enddo
       enddo
    enddo

! ----------------------------------------------------------------------
! Set spring natural lengths
! ----------------------------------------------------------------------

    xspring=0; yspring=0; zspring=0

    do i=0,imax-1
       xspring(i)=abs(xmesh(i+1)-xmesh(i))
    enddo

    do j=0,jmax-1
       yspring(j)=abs(ymesh(j+1)-ymesh(j))
    enddo
    
    do k=0,kmax-1
       zspring(k)=abs(zmesh(k+1)-zmesh(k))
    enddo
    
    do i=0,imax-1
       do j=0,jmax-1
          do k=0,kmax-1
             spring000(i,j,k)=mag(rnode(i  ,j,  k,:)-rnode(i+1,j+1,k+1,:))
             spring100(i,j,k)=mag(rnode(i+1,j,  k,:)-rnode(i,  j+1,k+1,:))
             spring110(i,j,k)=mag(rnode(i+1,j+1,k,:)-rnode(i,  j,  k+1,:))
             spring010(i,j,k)=mag(rnode(i  ,j+1,k,:)-rnode(i+1,j,  k+1,:))
          enddo
       enddo
    enddo

! ----------------------------------------------------------------------
! Find out how many parts there are. Allocate parts
! ----------------------------------------------------------------------

    if (allocated(parts)) then
       do i=1,size(parts,1)
          call freepart(parts(i))             ! would not be needed with TR15581
       enddo

       deallocate(parts)
    endif
          
    nparts=0
    do
       read (1,'(a)',iostat=ios) line
       
       if (ios.ne.0) exit
       
       call tolower(line)
       line=adjustl(line)
       
       if (line(1:1).eq.'*'.or.len_trim(line).eq.0) cycle
       
       if (line(1:4).eq.'part') nparts=nparts+1
    enddo
    
    allocate(parts(nparts),stat=allocerr)
    if (allocerr.ne.0) then
       call createerror('Could not allocate parts array')
       close (1)
       return
    endif

    rewind(1)
          
! ----------------------------------------------------------------------
! Read parts info from MIN file. Close MIN file
! ----------------------------------------------------------------------

    do i=1,nparts
       write (dum,'(i2.2)') i
       call createmessage('Reading part: '//trim(dum))
       if (getpart(parts(i),1).ne.0) then
          call createerror('Failed to read part')
          close (1)
          return
       endif
    enddo
    
    close (1)
          
! ----------------------------------------------------------------------
! Open and write MLS file
! ----------------------------------------------------------------------
    
    open (2,file=trim(stemname)//'.mls',status='unknown',iostat=ios)
    if (ios.ne.0) then
       call createerror('Could not open file:'//trim(stemname)//'.mls')
       return
    endif
    
    write (2,'(/a10,2a13)') 'i','xmesh','xspring'
    do i=0,imax
       write (2,'(i10,2e13.5)') i,xmesh(i),xspring(i)
    enddo
    
    write (2,'(/a10,2a13)') 'j','ymesh','yspring'
    do j=0,jmax
       write (2,'(i10,2e13.5)') j,ymesh(j),yspring(j)
    enddo
    
    write (2,'(/a10,2a13)') 'k','zmesh','zspring'
    do k=0,kmax
       write (2,'(i10,2e13.5)') k,zmesh(k),zspring(k)
    enddo
    
    do i=1,nparts
       call wrtpart(parts(i),2)
    enddo

! ----------------------------------------------------------------------
! Show body diagonals (debug only)
! ----------------------------------------------------------------------

!!$    do i=0,imax-1
!!$       do j=0,jmax-1
!!$          do k=0,kmax-1
!!$             write (2,'(3i5,4e13.5)') i,j,k,spring000(i,j,k),spring100(i,j,k),&
!!$                  spring010(i,j,k),spring110(i,j,k)
!!$          enddo
!!$       enddo
!!$    enddo

    close (2)

! ----------------------------------------------------------------------
! Generate mesh
! ----------------------------------------------------------------------
    
    call elements

    if (any(iregup(0:imax-1,0:jmax-1,0:kmax-1).eq.0)) then
       call createerror('An Element has undefined region number')
       return
    endif
    
    if (any(iregnd.eq.0)) then
       call createerror('A Node has undefined node number')
       return
    endif
    
    if (any(iprtup(0:imax-1,0:jmax-1,0:kmax-1).eq.0)) then
       call createerror('An Element has undefined part number')
       return
    endif
    
    if (any(iprtnd.eq.0)) then
       call createerror('A Node has undefined part number')
       return
    endif
          
    call snapnodes
    call coat
    call openparts

    if (fix) then
       call createmessage('Fixing inverted elements...')
       if (fixinv(.false.).gt.0) then
          call createmessage('Inverted elements remain, trying again...')
          if (fixinv(.true.).gt.0) then
             call createmessage('Some elements distorted: FE results may be innaccurate')
          endif
       endif
    endif
    
! ----------------------------------------------------------------------
! Set negative iregnd to positive
! ----------------------------------------------------------------------

    iregnd=abs(iregnd)

! ----------------------------------------------------------------------
! write out mesh
! ----------------------------------------------------------------------

    call createmessage('Writing file:')
    call createmessage(trim(stemname)//'.mtf')

    open (1,file=trim(stemname)//'.mtf',status='unknown',iostat=ios)
    if (ios.ne.0) then
       call createerror('Could not open file:'//trim(stemname)//'.mtf')
       return
    endif

    write (1,*) 'Written by Zmesh'
    write (1,*) 'iMAX: ',imax
    write (1,*) 'jMAX: ',jmax
    write (1,*) 'kMAX: ',kmax
    
    write (1,'(/3a10,3a13,2a10/)') 'i','j','k','x','y','z','regno','regup'
          
    do k=0,kmax
       do j=0,jmax
          do i=0,imax
             write (1,'(3i10,3e13.5,2i10)') i,j,k,rnode(i,j,k,:),iregnd(i,j,k),iregup(i,j,k)
          enddo
       enddo
    enddo
    
    close (1)

    processMIN=0              ! successful

  end function processMIN

! ######################################################################

  subroutine createfree

    integer i

    if (allocated(xmeshline)) then
       deallocate(xmeshline,ymeshline,zmeshline,xmesh,ymesh,zmesh)
    endif
    
    if (allocated(parts)) then
       do i=1,nparts
          call freepart(parts(i))
       enddo
       
       deallocate(parts)
    endif

    if (allocated(xspring)) then
       deallocate (xspring,yspring,zspring,spring000,spring100,spring110,spring010)
       deallocate (fixed,force)
    endif
  end subroutine createfree

! ######################################################################

  subroutine elements

!    use kernel32
!    use user32
    use global

    type(Tpart) part
    integer ipart,i,j,k,ioff,joff,koff,ret
    real rc(3)
    character(2) cpart

    iregup=0
    iregnd=0
    iprtnd=0
    iprtup=0

    do ipart=1,nparts

       part=parts(ipart)

       if (.not.isfilled(part)) cycle
       
       write (cpart,'(i2.2)') ipart
       call createmessage('Processing part: '//cpart)
              
       do i=0,imax-1
          do j=0,jmax-1
             do k=0,kmax-1
                
                rc=0
                do ioff=0,1
                   do joff=0,1
                      do koff=0,1
                         rc=rc+rnode(i+ioff,j+joff,k+koff,:)
                      enddo
                   enddo
                enddo
                
                rc=rc/8
                
                if (inpart(rc,part)) then
                   iregup(i,j,k)=part%region
                   iprtup(i,j,k)=ipart

                   do ioff=0,1
                      do joff=0,1
                         do koff=0,1
                            call setregnd(i+ioff,j+joff,k+koff,part%region,part%protect)
!                            iregnd(i+ioff,j+joff,k+koff)=part%region
                            iprtnd(i+ioff,j+joff,k+koff)=ipart
                         enddo
                      enddo
                   enddo

                endif
             enddo
          enddo
       enddo
    enddo
    
  end subroutine elements

! ######################################################################

  subroutine setregnd(i,j,k,region,protect)

! ----------------------------------------------------------------------
! Attempt to set node region number. If it is negative, we cannot
! change it (protected). IF this part is protected, set node reg to negative.
! ----------------------------------------------------------------------

    use global, only : iregnd

    integer i,j,k,region
    logical protect

    if (iregnd(i,j,k).lt.0) then
!       call createmessage('protected node')
       return
    endif

    if (protect) then
       iregnd(i,j,k)=-region
    else
       iregnd(i,j,k)=region
    endif
  end subroutine setregnd

! ######################################################################

  function bestpoint(i,j,k,r,part,flag,what,region)
    
! ----------------------------------------------------------------------
! Find the "best" point near to node point r on surface of part
! flag=1: on surface or edge; 2 on edge only
! i,j,k is the index to the node. 
! normally best point is just snappart but if (i,j,k) correspond to side
! or edge of simulation, we must change one or more of the components
! so that the box integrity is maintained.
! ----------------------------------------------------------------------

    use global, only : imax,jmax,kmax

    real r(3),bestpoint(3),maskpos(3),r2(3)
    type(Tpart) part
    integer i,j,k,flag,what,region,mask(3),iter

    mask=1
    maskpos=0

    if (i.eq.0) then
       mask(1)=0
       maskpos(1)=xmesh(0)
    else if (i.eq.imax) then
       mask(1)=0
       maskpos(1)=xmesh(imax)
    endif
    
    if (j.eq.0) then
       mask(2)=0
       maskpos(2)=ymesh(0)
    else if (j.eq.jmax) then
       mask(2)=0
       maskpos(2)=ymesh(jmax)
    endif
    
    if (k.eq.0) then
       mask(3)=0
       maskpos(3)=zmesh(0)
    else if (k.eq.kmax) then
       mask(3)=0
       maskpos(3)=zmesh(kmax)
    endif

    if (all(mask.eq.1)) then
       bestpoint=snappart(r,part,flag,what,region)
    else
       bestpoint=r
       do iter=1,5
          r2=snappart(bestpoint,part,flag,what,region)*mask+maskpos
          bestpoint=0.9*r2+0.1*bestpoint
       enddo
    endif
  end function bestpoint
  
! ######################################################################

  subroutine snapnodes

! ----------------------------------------------------------------------
! Snap nodes
! ----------------------------------------------------------------------

!    use kernel32
!    use user32
    use global

    type(Tpart) part

    integer ipart,ret,i,j,k,istep,what,region
    real r(3),factor,r2(3),dist(3),rdf
    logical doedge
    character(2) cpart
    character(20) cdum

    fixed=0

    do ipart=1,nparts

       part=parts(ipart)

       if (.not.isfilled(part)) cycle

       write (cpart,'(i2.2)') ipart
       
       if (part%nsurf.gt.0) then

          call createmessage('Surface fitting part: '//cpart)
          
          do i=0,imax
             do j=0,jmax
                do k=0,kmax

                   r=rnode(i,j,k,:)
                   
                   if (elligible(ipart,part,i,j,k,doedge,factor)) then
                      
                      if (.not.doedge) then
                      
!                         r2=snappart(r,part,1,what,region)
                         r2=bestpoint(i,j,k,r,part,1,what,region)

                         rnode(i,j,k,:)=(1-factor)*r+factor*r2
                         
                      else
                         
                         dist=localdist(i,j,k)                         
!                         r2=snappart(r,part,2,what,region)
                         r2=bestpoint(i,j,k,r,part,2,what,region)
                         
                         if (what.eq.1) then
                            rnode(i,j,k,:)=(1-factor)*r+factor*r2
                            
                         else if (what.ge.2) then
                            if (all(abs(r-r2).lt.dist)) then
                               rnode(i,j,k,:)=(1-factor)*r+factor*r2
                            else
!                               r2=snappart(r,part,1,what,region)
                               r2=bestpoint(i,j,k,r,part,1,what,region)
                               rnode(i,j,k,:)=(1-factor)*r+factor*r2
                            endif
                         endif
                       
                      endif    ! doedge
                      
                      fixed(i,j,k)=1
                   endif       ! elligible node?
                   
                enddo
             enddo
          enddo
          
          if (lrelax) then
             do istep=1,nstep
                call relax(.true.,rdf)   ! use diagonal springs
             enddo
          endif

       endif
    enddo
  end subroutine snapnodes

! ######################################################################

  function exists(i,j,k)

! ----------------------------------------------------------------------
! Element defined by lower left corner exists
! ----------------------------------------------------------------------

    use global, only : imax,jmax,kmax

    logical exists
    integer i,j,k

    if ( i.ge.0.and.i.lt.imax.and.&
         j.ge.0.and.j.lt.jmax.and.&
         k.ge.0.and.k.lt.kmax) then
       exists=.true.
    else
       exists=.false.
    endif
  end function exists

! ######################################################################

  function elligible(ipart,part,i,j,k,doedge,factor)

    use global

    type(Tpart) part

    logical elligible,doedge
    integer i,j,k,surftype,surfto,isurf,ioff,joff,koff,ipart
    real factor,surftol
    logical prtfound,region,edge

    elligible=.false.

    prtfound=.false.
    do ioff=-1,0
       do joff=-1,0
          do koff=-1,0
             if (exists(i+ioff,j+joff,k+koff)) then
                if (iprtup(i+ioff,j+joff,k+koff).eq.ipart) then
                   prtfound=.true.
                endif
             endif
          enddo
       enddo
    enddo

    if (.not.prtfound) return          ! elligible = .false.
    
    doedge=.false.
    factor=0
    
    do isurf=1,part%nsurf
       surftype=part%surftype(isurf)
       surfto=part%surfto(isurf)
       surftol=part%surftol(isurf)
       
       if (surftype.eq.3.or.surftype.eq.4) then
          edge=.true.
       else
          edge=.false.
       endif
       
       if (surftype.eq.1.or.surftype.eq.3) then
          region=.true.
       else
          region=.false.
       endif
       
       do ioff=-1,0
          do joff=-1,0
             do koff=-1,0

                if (exists(i+ioff,j+joff,k+koff)) then
                   if (region.and.iregup(i+ioff,j+joff,k+koff).eq.surfto) then
                      elligible=.true.
                      if (edge) doedge=.true.
                      factor=max(surftol,factor)
                      return
                   else if (.not.region.and.iprtup(i+ioff,j+joff,k+koff).eq.surfto) then
                      elligible=.true.
                      if (edge) doedge=.true.
                      factor=max(surftol,factor)
                      return
                   endif
                endif
             enddo
          enddo
       enddo
       
    enddo

! elligible = .false.

  end function elligible

! ######################################################################

  function fixinv(diag)

! ----------------------------------------------------------------------
! Check if any element is inverted and fix
! (elements are not sorted as in Zinc: we assume we can just take
! them out of the grid)
! if diag=.true. use diagonal springs in relaxation
! ----------------------------------------------------------------------

    use global
    logical diag
    integer fixinv,ipass,itot,dethist(8),i,j,k,n
    character(100) mess
    real rdf

    do ipass=1,nstep
       
       itot=0
       dethist=0
       
       do i=0,imax-1
          do j=0,jmax-1
             do k=0,kmax-1
                
                n=nnegdet(i,j,k)
                
                if (n.gt.0) then
                   
                   dethist(n)=dethist(n)+1
                   
                   itot=itot+1
                   
                   fixed(i,j,k)=0
                   fixed(i+1,j,k)=0
                   fixed(i,j+1,k)=0
                   fixed(i+1,j+1,k)=0
        
                   fixed(i,j,k+1)=0
                   fixed(i+1,j,k+1)=0
                   fixed(i,j+1,k+1)=0
                   fixed(i+1,j+1,k+1)=0
                   
                endif
                
             enddo
          enddo
       enddo
       
       write (mess,'(a,i5,a,i5,a)') 'pass ',ipass,' : ',itot,' inverted elements'
       call createmessage(mess)

       if (itot.eq.0) exit
       
       call relax(diag,rdf)
       
    enddo
    
    fixinv=itot

  end function fixinv

! ######################################################################

  subroutine lmesh(meshline,mesh)
    
    real meshline(:,:),mesh(0:)
    real z1,z2,dz
    integer ind,i,nmeshline,nz,j

    nmeshline=size(meshline,1)
    
    ind=0
    do i=1,nmeshline
       
       z1=meshline(i,1)
       z2=meshline(i,2)
       dz=meshline(i,3)
       nz=nint((z2-z1)/dz)
       
       do j=1,nz
          mesh(ind)=z1+(j-1)*dz
          ind=ind+1
       enddo
    enddo
    
    mesh(ind)=meshline(nmeshline,2)
  end subroutine lmesh

! ######################################################################

  subroutine relax(diag,rdf)

    use global

    logical diag
    real rdf,r1(3),r2(3),r12(3),length,r12hat(3),f1(3),length0,lengthdiag
    integer i,j,k,mask(3)
    
    force=0

! ----------------------------------------------------------------------
! Straight springs
! ----------------------------------------------------------------------

    do i=0,imax-1
       length0=xspring(i)

       do k=0,kmax
          do j=0,jmax

             r1=rnode(i,j,k,:)
             r2=rnode(i+1,j,k,:)
          
             r12=r2-r1
             
             length=mag(r12)

             r12hat=r12/length
             
             f1=konst*(length-length0)*r12hat

             force(i,j,k,:)=force(i,j,k,:)+f1
             force(i+1,j,k,:)=force(i+1,j,k,:)-f1
          enddo
       enddo
    enddo

    do j=0,jmax-1
       length0=yspring(j)

       do k=0,kmax
          do i=0,imax

             r1=rnode(i,j,k,:)
             r2=rnode(i,j+1,k,:)
          
             r12=r2-r1
             
             length=mag(r12)

             
             r12hat=r12/length
             
             f1=konst*(length-length0)*r12hat
             
             force(i,j,k,:)=force(i,j,k,:)+f1
             force(i,j+1,k,:)=force(i,j+1,k,:)-f1
          enddo
       enddo
    enddo

    do k=0,kmax-1
       length0=zspring(k)

       do j=0,jmax
          do i=0,imax

             r1=rnode(i,j,k,:)
             r2=rnode(i,j,k+1,:)
          
             r12=r2-r1
             
             length=mag(r12)

             r12hat=r12/length
             
             f1=konst*(length-length0)*r12hat
             
             force(i,j,k,:)=force(i,j,k,:)+f1
             force(i,j,k+1,:)=force(i,j,k+1,:)-f1
          enddo
       enddo
    enddo

! ----------------------------------------------------------------------
! Add body diagonal springs
! ----------------------------------------------------------------------

    if (diag) then

    do i=0,imax-1
       do j=0,jmax-1
          do k=0,kmax-1
             
             r1=rnode(i,j,k,:)
             r2=rnode(i+1,j+1,k+1,:)
          
             r12=r2-r1
             length=mag(r12)
             r12hat=r12/length

             lengthdiag=spring000(i,j,k)

             f1=konst*(length-lengthdiag)*r12hat

             force(i,j,k,:)=force(i,j,k,:)+f1
             force(i+1,j+1,k+1,:)=force(i+1,j+1,k+1,:)-f1

! ----------------------------------------------------------------------

             r1=rnode(i+1,j,k,:)
             r2=rnode(i,j+1,k+1,:)
          
             r12=r2-r1
             length=mag(r12)
             r12hat=r12/length

             lengthdiag=spring100(i,j,k)

             f1=konst*(length-lengthdiag)*r12hat
             
             force(i+1,j,k,:)=force(i+1,j,k,:)+f1
             force(i,j+1,k+1,:)=force(i,j+1,k+1,:)-f1

! ----------------------------------------------------------------------

             r1=rnode(i+1,j+1,k,:)
             r2=rnode(i,j,k+1,:)
          
             r12=r2-r1
             length=mag(r12)
             r12hat=r12/length

             lengthdiag=spring110(i,j,k)

             f1=konst*(length-lengthdiag)*r12hat

             force(i+1,j+1,k,:)=force(i+1,j+1,k,:)+f1
             force(i,j,k+1,:)=force(i,j,k+1,:)-f1

! ----------------------------------------------------------------------

             r1=rnode(i,j+1,k,:)
             r2=rnode(i+1,j,k+1,:)
          
             r12=r2-r1
             length=mag(r12)
             r12hat=r12/length

             lengthdiag=spring010(i,j,k)

             f1=konst*(length-lengthdiag)*r12hat

             force(i,j+1,k,:)=force(i,j+1,k,:)+f1
             force(i+1,j,k+1,:)=force(i+1,j,k+1,:)-f1

          enddo
       enddo
    enddo

    endif

! ----------------------------------------------------------------------
! Zero forces on edges atoms
! ----------------------------------------------------------------------

    do i=0,imax
       do j=0,jmax
          do k=0,kmax

             mask=[1,1,1]
             if (i.eq.0.or.i.eq.imax) mask(1)=0
             if (j.eq.0.or.j.eq.jmax) mask(2)=0
             if (k.eq.0.or.k.eq.kmax) mask(3)=0

             force(i,j,k,:)=force(i,j,k,:)*mask
          enddo
       enddo
    enddo

! ----------------------------------------------------------------------
! Update node positions
! ----------------------------------------------------------------------

    rdf=0
    do i=0,imax
       do j=0,jmax
          do k=0,kmax

             if (fixed(i,j,k).eq.0) then
                rnode(i,j,k,:)=rnode(i,j,k,:)+tstep*force(i,j,k,:)
                rdf=rdf+mag(tstep*force(i,j,k,:))
             endif
          enddo
       enddo
    enddo

!    rdf=rdf/nnod
  end subroutine relax

! ######################################################################

  function nnegdet(i,j,k)

    use global, only : rnode
    integer i,j,k,icorn,nnegdet

    real, parameter :: xic(8)= [-1,1,1,-1,-1,1,1,-1]
    real, parameter :: etac(8)=[-1,-1,1,1,-1,-1,1,1]
    real, parameter :: muc(8)= [-1,-1,-1,-1,1,1,1,1]
    
    real mat(3,3),drdxi(3),drdeta(3),drdmu(3),det,xi,eta,mu
    real, dimension(3) :: r1,r2,r3,r4,r5,r6,r7,r8

    r1=rnode(i,j,k,:)
    r2=rnode(i+1,j,k,:)
    r3=rnode(i+1,j+1,k,:)
    r4=rnode(i,j+1,k,:)
    
    r5=rnode(i,j,k+1,:)
    r6=rnode(i+1,j,k+1,:)
    r7=rnode(i+1,j+1,k+1,:)
    r8=rnode(i,j+1,k+1,:)

    nnegdet=0
        
    do icorn=1,8
       xi=xic(icorn)
       eta=etac(icorn)
       mu=muc(icorn)

       drdxi= 0.125*(-r1*(1-eta)*(1-mu)+r2*(1-eta)*(1-mu)+r3*(1+eta)*(1-mu)-r4*(1+eta)*(1-mu) &
                     -r5*(1-eta)*(1+mu)+r6*(1-eta)*(1+mu)+r7*(1+eta)*(1+mu)-r8*(1+eta)*(1+mu))

       drdeta=0.125*(-r1*(1-xi )*(1-mu)-r2*(1+xi )*(1-mu)+r3*(1+xi )*(1-mu)+r4*(1-xi )*(1-mu) &
                     -r5*(1-xi )*(1+mu)-r6*(1+xi )*(1+mu)+r7*(1+xi )*(1+mu)+r8*(1-xi )*(1+mu))

       drdmu= 0.125*(-r1*(1-xi)*(1-eta)-r2*(1+xi)*(1-eta)-r3*(1+xi)*(1+eta)-r4*(1-xi)*(1+eta) &
                     +r5*(1-xi)*(1-eta)+r6*(1+xi)*(1-eta)+r7*(1+xi)*(1+eta)+r8*(1-xi)*(1+eta))

       mat(1,:)=drdxi
       mat(2,:)=drdeta
       mat(3,:)=drdmu
        
       det=  mat(1,1)*(mat(2,2)*mat(3,3)-mat(3,2)*mat(2,3)) &
            -mat(1,2)*(mat(2,1)*mat(3,3)-mat(3,1)*mat(2,3)) &
            +mat(1,3)*(mat(2,1)*mat(3,2)-mat(3,1)*mat(2,2))

       if (det.lt.0) nnegdet=nnegdet+1
    enddo
  end function nnegdet

! ######################################################################

  subroutine coat
    
    use global

    type(Tpart) part
    integer ipart,i,j,k,ioff,joff,koff,newreg,icoat

    logical prtfound,regfound
    character(2) cdum

    do ipart=1,nparts

       part=parts(ipart)

       if (part%ncoat.gt.0) then

          write (cdum,'(i2.2)') ipart
          call createmessage('Coating part: '//cdum)

          do i=0,imax
             do j=0,jmax
                do k=0,kmax
                   
                   prtfound=.false.
                   do ioff=-1,0
                      do joff=-1,0
                         do koff=-1,0
                            if (exists(i+ioff,j+joff,k+koff)) then
                               if (iprtup(i+ioff,j+joff,k+koff).eq.ipart) then
                                  prtfound=.true.
                               endif
                            endif
                         enddo
                      enddo
                   enddo
                 
                   regfound=.false.
                   do icoat=1,part%ncoat
                      do ioff=-1,0
                         do joff=-1,0
                            do koff=-1,0
                               
                               if (exists(i+ioff,j+joff,k+koff)) then
                                  if (part%coatdata(icoat,1).eq.iregup(i+ioff,j+joff,k+koff)) then
                                     regfound=.true.
                                     newreg=part%coatdata(icoat,2)
                                  endif
                               endif
                            enddo
                         enddo
                      enddo
                   enddo
                 
                   if (regfound.and.prtfound) then
                      iregnd(i,j,k)=newreg          ! COAT ignores protection
                   endif
                enddo
             enddo
          enddo
       endif
    enddo
  end subroutine coat

! ######################################################################

  subroutine openparts

    use global

    type(Tpart) part
    integer ipart,i,j,k,imin,jmin,kmin,what,region
    real r(3),r0(3),dist(3),dist2,distmin2,r2(3)
    logical first
    character(2) cdum
    character(200) string

    do ipart=1,nparts

       part=parts(ipart)

       if (isfilled(part)) cycle

       write (cdum,'(i2.2)') ipart
       call createmessage('Processing OPEN part: '//cdum)

       if (index(part%type,'bound').ne.0) then
        
          if (index(part%type,'x').ne.0) then
           
             if (index(part%type,'up').ne.0) then
                do j=0,jmax
                   do k=0,kmax
!                      iregnd(imax,j,k)=part%region
                      call setregnd(imax,j,k,part%region,.false.)
                   enddo
                enddo
             else
                do j=0,jmax
                   do k=0,kmax
!                      iregnd(0,j,k)=part%region
                      call setregnd(0,j,k,part%region,.false.)
                   enddo
                enddo
             endif

          else if (index(part%type,'y').ne.0) then
             
             if (index(part%type,'up').ne.0) then
                do i=0,imax
                   do k=0,kmax
!                      iregnd(i,jmax,k)=part%region
                      call setregnd(i,jmax,k,part%region,.false.)
                   enddo
                enddo
             else
                do i=0,imax
                   do k=0,kmax
!                      iregnd(i,0,k)=part%region
                      call setregnd(i,0,k,part%region,.false.)
                   enddo
                enddo
             endif
           
          else

             if (index(part%type,'up').ne.0) then
                do i=0,imax
                   do j=0,jmax
!                      iregnd(i,j,kmax)=part%region
                      call setregnd(i,j,kmax,part%region,.false.)
                   enddo
                enddo
             else
                do i=0,imax
                   do j=0,jmax
!                      iregnd(i,j,0)=part%region
                      call setregnd(i,j,0,part%region,.false.)
                   enddo
                enddo
             endif
          endif

       else if (part%type.eq.'point') then

          r0=part%shift(1:3)
          
          first=.true.
        
          do i=0,imax
             do j=0,jmax
                do k=0,kmax
                   r=rnode(i,j,k,:)
                   dist2=(r(1)-r0(1))**2+(r(2)-r0(2))**2+(r(3)-r0(3))**2
                 
                   if (first) then
                      first=.false.
                      distmin2=dist2
                      imin=i; jmin=j; kmin=k
                   else if (dist2.lt.distmin2) then
                      distmin2=dist2
                      imin=i; jmin=j; kmin=k
                   endif
                enddo
             enddo
          enddo
          
!          rnode(imin,jmin,kmin,:)=r0
!          iregnd(imin,jmin,kmin)=part%region
          call setregnd(imin,jmin,kmin,part%region,.false.)
          fixed(imin,jmin,kmin)=1
           
       else if (is1d(part)) then
          do i=0,imax
             do j=0,jmax
                do k=0,kmax
                   r=rnode(i,j,k,:)
                   dist=localdist(i,j,k)
                   
                   r2=snappart(r,part,2,what,region)
                   if (all(abs(r-r2).lt.dist)) then
!                    rnode(i,j,k,:)=(1-factor)*r+r2*factor
!                      iregnd(i,j,k)=part%region
                      call setregnd(i,j,k,part%region,.false.)
                      fixed(i,j,k)=1
                   endif
                enddo
             enddo
          enddo
       else                  ! 2-D
          do i=0,imax
             do j=0,jmax
                do k=0,kmax
                   r=rnode(i,j,k,:)
                   dist=localdist(i,j,k)*0.2
                   
                   r2=snappart(r,part,2,what,region)       ! try with edge

! ----------------------------------------------------------------------
!                   if (i.eq.imax.and.j.eq.jmax.and.k.eq.kmax) then
!                      write (string,*) ipart,mag(r2-r),dist,what
!                      call createmessage(string)
!                   endif
! ----------------------------------------------------------------------

                   if (all(abs(r-r2).lt.dist)) then
                      !                    rnode(i,j,k,:)=(1-factor)*r+r2*factor
!                      iregnd(i,j,k)=part%region
                      call setregnd(i,j,k,part%region,.false.)
                      fixed(i,j,k)=1
                   else if (what.ne.1) then
                      r2=snappart(r,part,1,what,region)    ! try again with surface
                      if (all(abs(r2-r).lt.dist)) then
                         !                       rnode(i,j,k,:)=(1-factor)*r+r2*factor
!                         iregnd(i,j,k)=part%region
                         call setregnd(i,j,k,part%region,.false.)
                         fixed(i,j,k)=1
                      endif
                   endif
                   
                enddo
             enddo
          enddo
       endif

    enddo
  end subroutine openparts

! ######################################################################

  function localdist(i,j,k)

    use global

    real localdist(3)
    integer i,j,k,ioff,joff,koff
    real drmax(3),rn(3),r(3)
    logical isin

    r=rnode(i,j,k,:)

    drmax=0
    do ioff=-1,1
       do joff=-1,1
          do koff=-1,1
             if (.not.(ioff.eq.0.and.joff.eq.0.and.koff.eq.0)) then

                isin=i+ioff.ge.0.and.i+ioff.le.imax.and.&
                     j+joff.ge.0.and.j+joff.le.jmax.and.&
                     k+koff.ge.0.and.k+koff.le.kmax

                if (isin) then
                   rn=rnode(i+ioff,j+joff,k+koff,:)

                   drmax=max(abs(rn-r),drmax)
!                   drmax=max(rn-r,drmax)
                endif
             endif
          enddo
       enddo
    enddo
    
!    localdist=mag(drmax)/2
    localdist=drmax/2
  end function localdist

! In future, should change this so that localdist returns a vector
! requirement: all(r-r2.lt.localdist) must be TRUE
! DONE: but concerned about 0.2 magic number

! Add alternative fixinv based on iteratively moving nodes to 
! increase jacobian of each node (in each element)

! ######################################################################

  function mag(r)
    real r(3),mag
    mag=sqrt(r(1)**2+r(2)**2+r(3)**2)
  end function mag

end module createmesh
