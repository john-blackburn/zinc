! Part of Zinc FE package. Author: John Blackburn

program zpp

! ----------------------------------------------------------------------
! output series of line/plane scans from zinc output.
! key_plot=1
! linescan 100 0.0 0.0 -0.0125 0.0 0.0 0.0125 = Vx^2+Vy^2+Vz^2
! planescan 1 0.0 -1 1 -1 1 20 20 = V^2+1 &
! = Vx*{1=eps1,2=eps2}. Here, {1=eps1,2=eps2} is
!
! if (region.eq.1) then
!    eps1
! else if (region.eq.2) then
!    eps2
! endif
!
! (could consider moving readcon and readmtf to iofile module)
! ----------------------------------------------------------------------

  use common
  use iofile
  use strings, only : lowercase
  use scan
  use evaluate
  use matrices, only : readmatrix,usersetup,freematrix,userfree,lookind

  implicit none

  type(Texpr), allocatable :: exprlist(:)
  
  character cdum*20

  integer nexpr,ixpr,ind,key_debug
  
  character(1000) line,lefts,rights,vals,output
  character(20) temp,sfrom,sto
  character(2) sxpr,inorms

  integer istatus,length,lfname,i,j,k,term,ios,itot,key_plot
  integer allocerr(6),ic,jc,kc,mreg,mregup,Nline
  integer Na,Nb,ii,i1,j1,k1,ii1,inorm,ifrom,ito,iregion,ipnorm,ieq

  logical ftol_set,key_plot_set,inbrak
  double precision x,y,z,x1,y1,z1,x2,y2,z2,pnorm,a1,a2,b1,b2
  integer numbrak,maxperbrak,numperbrak

  integer, allocatable :: keyexpr(:),iparexpr(:,:)
  double precision, allocatable :: parexpr(:,:)

! ----------------------------------------------------------------------
! Compiled settings, not available to user
! ----------------------------------------------------------------------

  key_debug=0

! ----------------------------------------------------------------------
! get filename
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

  write (*,*) 'This is '//trim(name)//' (zpp)'
  print *,'Processing files beginning:'
  print *,trim(fname)

  date=' '
  time=' '
  call date_and_time(date,time)

! ----------------------------------------------------------------------
! Read CAFQG matrices from ZIN/CON file. Allocate and populate
! appropriate arrays. Also, set global ZIN variables.
! Set US, DUS and DEFPARAM constants and default variables
! ----------------------------------------------------------------------

  call readmatrix
  call defparam('NX',0.0)
  call defparam('NY',0.0)
  call defparam('NZ',0.0)

  ng_xi=3; ng_eta=3; ng_mu=3   ! only 3d for now
  ndim=3
  
! ----------------------------------------------------------------------
!     Read from MTF file, multiply positions by factor "scale"
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.mtf',term,'old','ASCII')

  itot=0
  do
     call getline(1,ios,line)
     if (ios.ne.0) call prterr('Unexpected end of file while searching for imax,jmax,kmax')

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
!     Allocate arrays u,rnode,ireg,iregup,labur,labdur
! ----------------------------------------------------------------------

  allocate (u(0:imax,0:jmax,0:kmax,nvar),stat=allocerr(1))
  allocate (rnode(0:imax,0:jmax,0:kmax,3),stat=allocerr(2))
  allocate (ireg(0:imax,0:jmax,0:kmax),stat=allocerr(3))
  allocate (iregup(0:imax,0:jmax,0:kmax),stat=allocerr(4))
  
  do i=1,4
     if (allocerr(i).ne.0) then
        write (temp,*) i
        call prterr('Failed to allocate array '//trim(temp))
     endif
  enddo

! ----------------------------------------------------------------------
! Read geometry data. scale up
! ----------------------------------------------------------------------

  do k=0,kmax
     do j=0,jmax
        do i=0,imax

           do
              call getline(1,ios,line)
              if (index(line,'REGNO').eq.0.and.index(line,'===').eq.0) exit
           enddo

           read (line,*,iostat=ios) ic,jc,kc,x,y,z,mreg,mregup

           if (ios.ne.0) call prterr('Incomprehensible mesh file line')

           if (.not.(i.eq.ic.and.j.eq.jc.and.k.eq.kc)) then
              call prterr('Unexpected node: vary i first, then j, then k')
           endif
           
           if ( .not.(mreg  .ge.1.and.mreg  .le.nregnd &
                .and. mregup.ge.0.and.mregup.le.nregel)) then
              call prterr('Region number out of range in .mtf file')
           endif
           
           rnode(i,j,k,1)=x*scale
           rnode(i,j,k,2)=y*scale
           rnode(i,j,k,3)=z*scale
           ireg(i,j,k)=mreg
           iregup(i,j,k)=mregup
        enddo
     enddo
  enddo
  
  close (1)
  write (*,*) 'Done.'

! ----------------------------------------------------------------------
! Pre-read .zpp file just to find how many expressions we need to store
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.zpp',term,'old','ASCII')

  nexpr=0
  do
     call getline(1,ios,line)
     if (ios.ne.0) exit

     if (index(line,'LINESCAN').ne.0.or.index(line,'PLANESCAN').ne.0 &
          .or.index(line,'SURFINT').ne.0.or.index(line,'VOLINT').ne.0 &
          .or.index(line,'FDCHK').ne.0) then
        nexpr=nexpr+1
     endif
  enddo

  close (1)

  allocate(keyexpr(nexpr),parexpr(nexpr,6),iparexpr(nexpr,6),stat=allocerr(1))
  allocate(exprlist(nexpr),stat=allocerr(2))

  do i=1,2
     if (allocerr(i).ne.0) then
        write (temp,*) i
        call prterr('Failed to allocate array '//trim(temp))
     endif
  enddo

! ----------------------------------------------------------------------
! Read .zpp again to find how many brackets are needed for each expression
! and the maximum number of entries per bracket for each line.
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.zpp',term,'old','ASCII')

  ind=0
  do
     call getline(1,ios,line)
     if (ios.ne.0) exit

     if (index(line,'LINESCAN').ne.0.or.index(line,'PLANESCAN').ne.0.or. &
          index(line,'SURFINT').ne.0.or.index(line,'VOLINT').ne.0) then
        i=spliteq(line,lefts,rights)

        if (i.eq.0) call prterr('Expecting = in scan line')

        ind=ind+1

        numbrak=0
        maxperbrak=0

        numperbrak=0
        inbrak=.false.

        do j=1,len_trim(rights)
           if (rights(j:j).eq.'{') then
              inbrak=.true.
              numbrak=numbrak+1
              
           else if (rights(j:j).eq.'}') then
              numperbrak=numperbrak+1
              maxperbrak=max(numperbrak,maxperbrak)
              numperbrak=0
              inbrak=.false.
              
           else if (inbrak.and.rights(j:j).eq.',') then
              numperbrak=numperbrak+1
           endif
        enddo

        exprlist(ind)%nbrak=numbrak

!        print *,'bracket storage',ind,numbrak,maxperbrak

        if (numbrak.gt.0) then
           allocate(exprlist(ind)%nperbrak(numbrak),        stat=allocerr(1))

           allocate(exprlist(ind)%sbrak(numbrak,maxperbrak),stat=allocerr(2))
           allocate(exprlist(ind)% brak(numbrak,maxperbrak),stat=allocerr(3))
           allocate(exprlist(ind)%ibrak(numbrak,maxperbrak),stat=allocerr(4))
           allocate(exprlist(ind)%tbrak(numbrak,maxperbrak),stat=allocerr(5))
           allocate(exprlist(ind)%rbrak(numbrak,maxperbrak),stat=allocerr(6))

           do i=1,6
              if (allocerr(i).ne.0) call prterr('Failed to allocate {...} storage')
           enddo
        endif
     
     endif
  enddo

  close (1)

! ----------------------------------------------------------------------
! Read .zpp file. Scale up scan geometry.
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.zpp',term,'old','ASCII')

  write (*,*) trim(fname)//'.zpp'

  ind=0

  key_plot_set=.false.
  ftol_set=.false.

  key_plot=0
  ftol=1e-3
  
  do
     call getline(1,ios,line)

     if (ios.ne.0) exit

     if (index(line,'LINESCAN').ne.0) then
        i=spliteq(line,lefts,vals)

        if (i.eq.0) call prterr('Expected =')

        read (lefts,*,iostat=ios) cdum,nline,x1,y1,z1,x2,y2,z2
        if (ios.ne.0) call prterr('Could not interpret linescan command')

        ind=ind+1
        keyexpr(ind)=1
        parexpr(ind,1:6)=[x1,y1,z1,x2,y2,z2]*scale
        iparexpr(ind,1)=nline

        call getexpr(ind,vals,.true.,.false.,.false.)   ! fill exprlist(ind). deriv, norm, istep

     else if (index(line,'PLANESCAN').ne.0) then

        i=spliteq(line,lefts,vals)

        if (i.eq.0) call prterr('Expected =')

        read (lefts,*,iostat=ios) cdum,inorm,pnorm,a1,a2,b1,b2,Na,Nb
        if (ios.ne.0) call prterr('Could not interpret planescan command')

        ind=ind+1
        keyexpr(ind)=2
        parexpr(ind,1:5)=[pnorm,a1,a2,b1,b2]*scale
        iparexpr(ind,1:3)=[inorm,Na,Nb]

        call getexpr(ind,vals,.true.,.false.,.false.)  ! deriv, norm, istep

     else if (index(line,'SURFINT').ne.0) then
        
        i=spliteq(line,lefts,vals)
        if (i.eq.0) call prterr('Expected =')

        do i=1,len_trim(lefts)
           if (lefts(i:i).eq.'-') lefts(i:i)=' '
        enddo

        read (lefts,*,iostat=ios) cdum,sfrom,sto
        if (ios.ne.0) call prterr('Could not interpret surfint command')

        ifrom=lookind(sfrom)
        ito=lookind(sto)

        ind=ind+1
        keyexpr(ind)=3
        iparexpr(ind,1:2)=[ifrom,ito]

        call getexpr(ind,vals,.true.,.true.,.false.)                  ! deriv, norm, istep

     else if (index(line,'VOLINT').ne.0) then

        i=spliteq(line,lefts,vals)
        if (i.eq.0) call prterr('Expected =')

        read (lefts,*,iostat=ios) cdum,iregion
        if (ios.ne.0) call prterr('Could not interpret VOLINT command')

        ind=ind+1
        keyexpr(ind)=4
        iparexpr(ind,1)=iregion

        call getexpr(ind,vals,.true.,.false.,.false.)

     else if (index(line,'FDCHK').ne.0) then

        ! FDCHK 1 23 2;  check plane i=23, equation 2

        read (line,*,iostat=ios) cdum,inorms,ipnorm,ieq
        if (ios.ne.0) call prterr('Could not interpret FDCHK command')

        if (inorms.eq.'I') then
           inorm=1
        else if (inorms.eq.'J') then
           inorm=2
        else if (inorms.eq.'K') then
           inorm=3
        else
           call prterr('First param for FDCHK must be I, J or K, got: '//inorms)
        endif

        ind=ind+1
        keyexpr(ind)=5
        iparexpr(ind,1:3)=[inorm,ipnorm,ieq]

        exprlist(ind)%sexpr=line

!        print *,'found FDCHK'

     else if (index(line,'=').ne.0) then
        i=spliteq(line,lefts,rights)

        if (lefts.eq.'GRAPH') then
           if (rights.eq.'EPS') then
              key_plot=1
           else if (rights.eq.'EMF') then
              key_plot=2
           else if (rights.eq.'NONE') then
              key_plot=0
           else
              call prterr('GRAPH must be EPS, EMF or NONE')
           endif

           key_plot_set=.true.
        else if (lefts.eq.'FTOL') then
           read (rights,*,err=911) ftol
           ftol_set=.true.
        else
           call prterr('unknown command')
        endif
     endif
  enddo

  close (1)

! ----------------------------------------------------------------------
! Load the DLL if necessary. Set up Cfun... and maybe scanfun
! ----------------------------------------------------------------------

  call usersetup(any(exprlist%iexpr.eq.3))

! ----------------------------------------------------------------------
! write out expression etc for each scan line
! ----------------------------------------------------------------------

  if (key_debug.eq.1) then
     do ixpr=1,nexpr
        print *,'Expression: ',ixpr
        
        print *,'keyexpr',keyexpr(ixpr)
        print *,'parexpr',parexpr(ixpr,:)
        print *,'iparexpr',iparexpr(ixpr,:)
        
        print *,'sexpr',trim(exprlist(ixpr)%sexpr)
        print *,'expr',exprlist(ixpr)%expr
        print *,'iexpr',exprlist(ixpr)%iexpr
        print *,'nbrak',exprlist(ixpr)%nbrak
     enddo
  endif

! ----------------------------------------------------------------------
! Read in solutions from .zou. Process scan lines
! produce .out files for gnuplot plotting
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.zou',term,'old',zou_format)

  do ixpr=1,nexpr

     write (sxpr,'(i2.2)') ixpr

     term=1
     call initfile(2,trim(fname)//trim(sxpr)//'.out',term,'unknown','ASCII')
     write (2,'(a,a)') '# ',trim(exprlist(ixpr)%sexpr)

     rewind(1)
     if (zou_format=='BINARY') read (1)

     do
        if (zou_format=='ASCII') then
           call getline(1,ios,line)
        else
           read (1,iostat=ios) i
        endif

        if (ios.ne.0) exit

        if ( (zou_format=='ASCII'.and.index(line,'ISTEP')/=0).or. &
             (zou_format=='BINARY'.and.i<0)) then

           write (2,'(a)') '#'//trim(line)
        
           do i=0,imax
              do j=0,jmax
                 do k=0,kmax
                    do ii=1,nvar

                       if (zou_format=='ASCII') then
                          call getline(1,ios,line)
                          if (ios.ne.0) call prterr('Unexpected end of file')
                          read (line,*,iostat=ios) i1,j1,k1,ii1,u(i,j,k,ii)
                       else
                          read (1,iostat=ios) i1,j1,k1,ii1,u(i,j,k,ii)
                       endif

                       if (ios.ne.0) call prterr('Incomprehensible line in file')

                       if (.not.(i.eq.i1.and.j.eq.j1.and.k.eq.k1.and.ii.eq.ii1)) then
                          call prterr('Bad ZOU file. Not compatible with input files.')
                       endif

                    enddo
                 enddo
              enddo
           enddo

           if (keyexpr(ixpr).eq.1) then
              x1=parexpr(ixpr,1)
              y1=parexpr(ixpr,2)
              z1=parexpr(ixpr,3)
              x2=parexpr(ixpr,4)
              y2=parexpr(ixpr,5)
              z2=parexpr(ixpr,6)

              Nline=iparexpr(ixpr,1)

              write (2,'(a,3e13.5,a,3e13.5,a,i5,a)') '# LINESCAN (',x1,y1,z1,') to (',x2,y2,z2,'): ',Nline,' steps'
              write (2,'(5a13,5a6)') '#    distance','x','y','z','value','ifail','i','j','k','reg'
              call linescan(2,Nline,x1,y1,z1,x2,y2,z2,exprlist(ixpr))
           else if (keyexpr(ixpr).eq.2) then
              pnorm=parexpr(ixpr,1)
              a1=parexpr(ixpr,2)
              a2=parexpr(ixpr,3)
              b1=parexpr(ixpr,4)
              b2=parexpr(ixpr,5)

              inorm=iparexpr(ixpr,1)
              Na=iparexpr(ixpr,2)
              Nb=iparexpr(ixpr,3)

              if (inorm.eq.1) then
                 write (2,'(a,e13.5,a,2e13.5,a,2e13.5,a)') '# PLANESCAN x=',pnorm,'; y=[',a1,a2,']; z=[',b1,b2,']'
              else if (inorm.eq.2) then
                 write (2,'(a,e13.5,a,2e13.5,a,2e13.5,a)') '# PLANESCAN y=',pnorm,'; x=[',a1,a2,']; z=[',b1,b2,']'
              else
                 write (2,'(a,e13.5,a,2e13.5,a,2e13.5,a)') '# PLANESCAN z=',pnorm,'; x=[',a1,a2,']; y=[',b1,b2,']'
              endif

              write (2,'(4a13,5a6)') '#           x','y','z','value','ifail','i','j','k','reg'
              call planescan(2,inorm,pnorm,a1,a2,b1,b2,Na,Nb,exprlist(ixpr))
           else if (keyexpr(ixpr).eq.3) then
              ifrom=iparexpr(ixpr,1)
              ito=iparexpr(ixpr,2)

              write (2,'(a,i5,a,i5)') '# SURFINT: from',ifrom,' to',ito
              write (2,*) 'integral=',surfint(ifrom,ito,exprlist(ixpr))

           else if (keyexpr(ixpr).eq.4) then
              iregion=iparexpr(ixpr,1)

              write (2,'(a,i5)') '# VOLINT: region',iregion
              write (2,*) 'integral=',volint(iregion,exprlist(ixpr))

           else if (keyexpr(ixpr).eq.5) then
              inorm=iparexpr(ixpr,1)
              ipnorm=iparexpr(ixpr,2)
              ieq=iparexpr(ixpr,3)
              
              if (inorm.eq.1) then
                 write (2,'(a,i5,a,i5)') '# FDCHK: i=',ipnorm,' Equation: ',ieq
              else if (inorm.eq.2) then
                 write (2,'(a,i5,a,i5)') '# FDCHK: j=',ipnorm,' Equation: ',ieq
              else if (inorm.eq.3) then
                 write (2,'(a,i5,a,i5)') '# FDCHK: k=',ipnorm,' Equation: ',ieq
              endif

              call fdchk(2,inorm,ipnorm,ieq)

           endif

           write (2,*)
           write (2,*)

        endif
     enddo

     close (2)

  enddo

  close (1)

! ----------------------------------------------------------------------
! Output .gnu file for gnuplot
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.gnu',term,'unknown','ASCII')

  if (key_plot.eq.1) then
     write (1,*) 'set term post eps'
  else
     write (1,*) 'set term emf'
  endif

!  write (1,*) 'cd '//trim(fname)
  write (1,*) 'set nokey'
  write (1,*) 'set hidden'
  write (1,*) 'set view 30,30'
  write (1,*) 'set format x "%h"'
  write (1,*) 'set format y "%h"'
  write (1,*) 'set format z "%h"'
  write (1,*)

  do ixpr=1,nexpr

     write (sxpr,'(i2.2)') ixpr

     if (key_plot.eq.1) then
        write (output,'(a)') 'set output '''//trim(fname)//sxpr//'.eps'''
     else
        write (output,'(a)') 'set output '''//trim(fname)//sxpr//'.emf'''
     endif

     if (keyexpr(ixpr).eq.1) then                   ! LINESCAN
        write (1,'(a)') trim(output)
        write (1,'(a,3e13.5,a,3e13.5,a)') &
             'set title "(',parexpr(ixpr,1:3),') to (',parexpr(ixpr,4:6),')"'

        write (1,'(a)') 'set xlabel "scan distance"'
        write (1,'(a)') 'set ylabel "'//trim(exprlist(ixpr)%sexpr)//'" noenhanced'

        write (1,'(a)') 'plot '''//trim(fname)//sxpr//'.out'' u 1:5 w l'

     else if (keyexpr(ixpr).eq.2) then              ! PLANESCAN
        
        write (1,'(a)') trim(output)
        write (1,'(a)') 'set zlabel "'//trim(exprlist(ixpr)%sexpr)//'" offset -5,0 rotate by 90 noenhanced'

        if (iparexpr(ixpr,1).eq.1) then
           write (1,'(a,e13.5,a)') 'set title "x=',parexpr(ixpr,1),'"'
           write (1,'(a)') 'set xlabel "y"'
           write (1,'(a)') 'set ylabel "z"'
           write (1,'(a)') 'splot '''//trim(fname)//sxpr//'.out'' u 2:3:4 w l'

        else if (iparexpr(ixpr,1).eq.2) then
           write (1,'(a,e13.5,a)') 'set title "y=',parexpr(ixpr,1),'"'
           write (1,'(a)') 'set xlabel "x"'
           write (1,'(a)') 'set ylabel "z"'
           write (1,'(a)') 'splot '''//trim(fname)//sxpr//'.out'' u 1:3:4 w l'

        else if (iparexpr(ixpr,1).eq.3) then
           write (1,'(a,e13.5,a)') 'set title "z=',parexpr(ixpr,1),'"'
           write (1,'(a)') 'set xlabel "x"'
           write (1,'(a)') 'set ylabel "y"'
           write (1,'(a)') 'splot '''//trim(fname)//sxpr//'.out'' u 1:2:4 w l'
        endif

     endif
     write (1,*)
  enddo

  close (1)

! ----------------------------------------------------------------------
! If key_plot not zero, call gnuplot
! ----------------------------------------------------------------------

  if (key_plot.ne.0) then
     print '(a)','gnuplot "'//trim(fname)//'.gnu"'
     call system('gnuplot "'//trim(fname)//'.gnu"')
  endif

! ----------------------------------------------------------------------
! free memory. Free DLL (if attached)
! ----------------------------------------------------------------------

  call freematrix
  call userfree
  deallocate (us,dus,u,rnode,ireg,iregup,keyexpr,parexpr,iparexpr)

  do ixpr=1,nexpr
     if (exprlist(ixpr)%nbrak.ne.0) then
!        print *,'deallocate',ixpr
        deallocate (exprlist(ixpr)%nperbrak, exprlist(ixpr)%sbrak, exprlist(ixpr)%rbrak, &
             exprlist(ixpr)%brak, exprlist(ixpr)%ibrak, exprlist(ixpr)%tbrak)
     endif
  enddo
  
  deallocate (exprlist)

  print *,'Program zpp completed successfully'
  
  call exit(0)

911 call prterr('Unable to decipher number/expression in file')
  
contains

  subroutine getexpr(ind,vals,derivOK,normOK,istepOK)

! ----------------------------------------------------------------------
! fill exprlist(ind) based on vals
! note that %nbrak is already set
! planescan 1 0.0 -1 1 -1 1 20 20 = V^2+1
! OR = Vx*{1=eps1,2=eps2}.
! ----------------------------------------------------------------------

    logical derivOK,normOK,istepOK
    character(*), intent(in) :: vals
    integer, intent(in) :: ind

    integer ival,ii,ibraket,ist,ied,i,j,jst,ibrakexp,ieq,ireg
    double precision val
    type(Tlist) list

    character(1000) string,expression,brakstr,lefts,rights
    character(2) temp
    logical found

    allocate (list%token(numtok))

    ival=0
    do ii=1,nvar
       if (us(ii).eq.adjustl(vals)) then
          exprlist(ind)%iexpr=-ii
          exprlist(ind)%sexpr=vals
          return
       endif
    enddo

    string=vals
    ibraket=0

    do
       ist=0
       ied=0
       do i=1,len_trim(string)
          if (string(i:i).eq.'{') ist=i
          if (string(i:i).eq.'}') then
             ied=i
             exit
          endif
       enddo
       
       if (ist.ne.0.and.ied.eq.0.or.ist.eq.0.and.ied.ne.0) call prterr('Unmatched {}')
       
       if (ist.ne.0.and.ied.ne.0) then
          ibraket=ibraket+1
          write (temp,'(i2.2)') ibraket
          
          brakstr='#'//trim(string(ist+1:ied-1))//','
          
          jst=1
          ibrakexp=0
          do j=1,len_trim(brakstr)
             if (brakstr(j:j).eq.',') then
                expression=brakstr(jst+1:j-1)
                ibrakexp=ibrakexp+1

                ieq=spliteq(expression,lefts,rights)
                if (ieq.eq.0) call prterr('Expecting "=" within {..} expression')

                read (lefts,*,iostat=ios) ireg
                if (ios.ne.0) call prterr('Incomprehensible region number in {..} expression')
                
                if (ireg.lt.1.or.ireg.gt.nregel) call prterr('Region number out of range in {..} expression')

                exprlist(ind)%rbrak(ibraket,ibrakexp)=ireg
                exprlist(ind)%sbrak(ibraket,ibrakexp)=rights
                
                found=.false.
                do ii=1,nvar
                   if (us(ii).eq.rights) then
                      exprlist(ind)%ibrak(ibraket,ibrakexp)=-ii
                      found=.true.
                   endif
                enddo
                   
                if (.not.found) then
                   call parsval(rights,val,ival,list,derivOK,normOK,istepOK)
                      
                   exprlist(ind)% brak(ibraket,ibrakexp)=val
                   exprlist(ind)%ibrak(ibraket,ibrakexp)=ival
                   exprlist(ind)%tbrak(ibraket,ibrakexp)=list
                endif

                jst=j
             endif
          enddo

          exprlist(ind)%nperbrak(ibraket)=ibrakexp

          string=string(:ist-1)//'regexp'//temp//string(ied+1:)
          call defparam('regexp'//temp,0)
       else
          exit
       endif
    enddo

    call parsval(string,val,ival,list,derivOK,normOK,istepOK)
    
    exprlist(ind)%expr=val
    exprlist(ind)%sexpr=vals
    exprlist(ind)%texpr=list
    exprlist(ind)%iexpr=ival          ! will only use this if no {}

    deallocate (list%token)

!    print *,'getexpr:',val,ival
  end subroutine getexpr
  
end program zpp
