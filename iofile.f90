! Part of Zinc FE package. Author: John Blackburn

module iofile

! subroutine intifile
! subroutine Ctest
! subroutine aatest
! subroutine fftest
! subroutine parsval
! subroutine getline
! subroutine prterr
! subroutine spliteq
! subroutine range
! subroutine drange
! subroutine shapshot
! subroutine resid_snapshot
! pvshapshot
! abbrev
! Cindex, aindex, findex

! uses evalexpr and tokenize

  implicit none

  integer, private :: lineno
  character(1000), private :: currfile

contains

subroutine initfile(iunit,filename,term,status,fmt)

! ----------------------------------------------------------------------
! open file if it exists to unit iunit.
! if file does not exist and term=1, stop program
! if file does not exist and term=0, issue warning and set term=-1
! set lineno to zero and currfile to filename
! if status='unknown' etc, write strapline
! fmt is BINARY of ASCII depending on format of file
! ----------------------------------------------------------------------

  use common, only : name,date,time
  implicit none
  
  integer iunit,term,ios
  character(*) filename,status
  character(*) fmt

  lineno=0
  currfile=filename

  print *,'Opening file:'
  print '(5x,a)',trim(filename)

  if (fmt=="ASCII") then
     open (iunit,file=filename,status=status,iostat=ios,form='formatted')
  else
     open (iunit,file=filename,status=status,iostat=ios,form='unformatted')
  endif

  if (ios.ne.0) then
     if (term.eq.1) then
        print *,'Error: failed to open:'
        print *,trim(filename)
        call exit(1)
     else
        print *,'Warning: failed to open OPTIONAL file'
        print *,trim(filename)
        print *,'Proceeding anyway'
        term=-1
     endif
  endif

  if (status.eq.'unknown'.or.status.eq.'replace'.or.status.eq.'new') then

     if (fmt=='BINARY') then
        write (iunit) '# Prepared by ',name, &
             '. Date: ',date(7:8),'/',date(5:6),'/',date(1:4), &
             '. Time: ',time(1:2),':',time(3:4),':',time(5:6)
     else
        write (iunit,'(14a)') '# Prepared by ',name, &
             '. Date: ',date(7:8),'/',date(5:6),'/',date(1:4), &
             '. Time: ',time(1:2),':',time(3:4),':',time(5:6)
     endif

  endif

end subroutine initfile

! ######################################################################

subroutine Ctest(ii,jj,kk,ll,nvar)

  use common, only : NDIM

  implicit none
  logical ok

  integer ii,jj,kk,ll,nvar

  ok=  ii.ge.1.and.ii.le.nvar.and.jj.ge.1.and.jj.le.NDIM.and. &
       kk.ge.1.and.kk.le.nvar.and.ll.ge.1.and.ll.le.NDIM
      
  if (.not.ok) then
     call prterr('"C" specification indices out of range')
  endif
  
end subroutine Ctest

! ######################################################################

subroutine aatest(ii,jj,nvar)
  implicit none
  logical ok

  integer ii,jj,nvar
  
  ok=ii.ge.1.and.ii.le.nvar.and.jj.ge.1.and.jj.le.nvar
  
  if (.not.ok) then
     call prterr ('"a" specification indices out of range')
  endif
  
end subroutine aatest

! ######################################################################

subroutine fftest(ii,nvar)
  implicit none
  logical ok
  
  integer ii,nvar

  ok=ii.ge.1.and.ii.le.nvar
  
  if (.not.ok) then
     call prterr('"f" specification indices out of range')
  endif
  
end subroutine fftest

! ######################################################################

subroutine parsval(vals,val,ival,list,derivOK,normOK,istepOK)

! ----------------------------------------------------------------------
! given vals, evaluate to number if possible and store in val.
! Otherwise prepare list using tokenize function.
! parsval checks that expressions are legitimate and fails otherwise
! set ival to 0 (const), 1 (contains xyz), 2 (contains u1 etc)
! 3 starts with $. In this case, function will be called
! note that starting blanks have already been removed
! ----------------------------------------------------------------------

  use common, only : us,dus,rs,restart,scoord,nvar,nvar_rst,NDIM
  use evaluate

  implicit none
  
  logical derivOK,normOK,istepOK
  character, intent(in) :: vals*(*)
  double precision, intent(out) :: val
  integer, intent(out) :: ival
  type(Tlist) :: list

  character(len(vals)) :: temp
  type(item) :: tok
  integer ifound,ivar,j,icp,isp

  if (vals(1:1).eq.'$') then
     val=0
     ival=3
  else

! ----------------------------------------------------------------------
! Check if expression contains forbidden tokens
! ----------------------------------------------------------------------

     ifound=0
     temp=vals
     do
        call get_next_token(temp,tok,icp,isp)
        if (tok%type.eq.'E') exit

        do ivar=1,nvar
           do j=1,NDIM
              if (tok%char.eq.dus(ivar,j).and.(.not.derivOK)) then
                 call prterr('This expression cannot contain derivatives')
              endif
           enddo
        enddo

        if ((.not.normOK).and.(tok%char.eq.'NX'.or.tok%char.eq.'NY'.or.tok%char.eq.'NZ')) then
           call prterr('This expression cannot contain unit normal (NX, NY, NZ)')
        endif

        if (tok%char.eq.'ISTEP'.and.(.not.istepOK)) then
           call prterr('This expression cannot contain ISTEP')
        endif
     enddo

! ----------------------------------------------------------------------
! Check if expression makes sense
! ----------------------------------------------------------------------

     call evalexpr(vals,val)
     if ( ierr.ne.0 .and.ierr.ne. 8.and.ierr.ne.9.and. &
          ierr.ne.10.and.ierr.ne.14) then
        call prterr('Error parsing expression')
     endif

! ----------------------------------------------------------------------
! IF expression contains variables and derivatives: non-linear, type 2
! ----------------------------------------------------------------------

     ifound=0
     temp=vals
     do
        call get_next_token(temp,tok,icp,isp)
        if (tok%type.eq.'E') exit

        if (restart/='NONE') then
           do ivar=1,nvar_rst
              if (tok%char.eq.rs(ivar)) ifound=1
           enddo
        endif

        do ivar=1,nvar
           if (tok%char.eq.us(ivar)) ifound=1
           do j=1,NDIM
              if (tok%char.eq.dus(ivar,j)) ifound=1
           enddo
        enddo
     enddo

     if (ifound.eq.1) then
        ival=2
        val=0
        call tokenize(vals,list)
        return
     endif

! ----------------------------------------------------------------------
! else if contains X, Y, Z, NX, NY, NZ, type 1.
! ----------------------------------------------------------------------

     ifound=0
     temp=vals
     do
        call get_next_token(temp,tok,icp,isp)
        if (tok%type.eq.'E') exit

        do j=1,NDIM
           if (tok%char.eq.scoord(j)) ifound=1
        enddo

        if (tok%char.eq.'NX') ifound=1
        if (tok%char.eq.'NY') ifound=1
        if (tok%char.eq.'NZ') ifound=1
     enddo

     if (ifound.eq.1) then
        ival=1
        val=0
        call tokenize(vals,list)
        return
     endif

! ----------------------------------------------------------------------
! no vars and no xyz, norm, so this is contant which we have already
! evaluated to val. type 0
! ----------------------------------------------------------------------

     ival=0

  endif

end subroutine parsval

! ######################################################################

subroutine getline(iunit,ios,line)

! ----------------------------------------------------------------------
! read the next non-blank, non-comment line from iunit
! (return with ios negative if end of file encountered)
! If line ends with & read next line and concatenate.
! Return total line in line
!
! Also, updates lineno and refers to currfile, both global
! ----------------------------------------------------------------------

  use strings
  implicit none

  integer iunit,ios
  character line*(*)
  character(len(line)) temp
  integer i,icont

  line=' '
  icont=0

  do
     lineno=lineno+1
     read (iunit,'(a)',iostat=ios) temp

     if (ios.gt.0) then
        print *,'Error reading from file:'
        print *,trim(currfile)
        print *,'line no:',lineno
        call exit(1)
     else if (ios.lt.0) then
        if (icont.eq.1) then
           print *,'Expecting continuation of line'
           print *,trim(currfile)
           print *,'line no:',lineno
           call exit(1)
        endif
        return
     endif

     temp=uppercase(temp)
     temp=adjustl(temp)
     
     i=index(temp,'!')
  
     if (i.gt.1) temp=temp(:i-1)

     if (temp(1:1).ne.'!'.and.len_trim(temp).ne.0) then
        i=index(temp,'&')
        if (i.eq.0) then
           line=trim(line)//trim(temp)
           return
        else
           icont=1
           temp=temp(:i-1)
           line=trim(line)//trim(temp)
        endif
     endif

  enddo
  
end subroutine getline

! ######################################################################

subroutine prterr(string)

  implicit none

  character string*(*)

  print '(a)',' ERROR: '//trim(string)
  print '(a,a)',' File: ',trim(currfile)
  print '(a,i5)',' Line: ',lineno

  call exit(1)
end subroutine prterr

! ######################################################################

function spliteq(line,left,right)
  implicit none
  integer spliteq
  character(len=*) line,left,right

  spliteq=index(line,'=')
  if (spliteq.ne.0) then
     left=line(:spliteq-1)
     right=line(spliteq+1:)
     
     left=adjustl(left)
     right=adjustl(right)
  endif

end function spliteq

! ######################################################################

subroutine range(ivar,name,imin,imax)
  implicit none
  integer ivar,imin,imax
  character name*(*)

  if (.not.(ivar.ge.imin.and.ivar.le.imax)) then
     print *,'Variable out of range:',trim(name)
     print *,'Should be ',imin,'to',imax
     call exit(1)
  endif
end subroutine range

! ######################################################################

subroutine drange(var,name,dmin,dmax)
  implicit none
  double precision var,dmin,dmax
  character name*(*)

  if (.not.(var.ge.dmin.and.var.le.dmax)) then
     print *,'Variable out of range:',trim(name)
     print *,'Should be ',dmin,'to',dmax
     call exit(1)
  endif
end subroutine drange

! ######################################################################

subroutine snapshot(iunit)

! ----------------------------------------------------------------------
! write current system state to iunit which must already be open
! Note: this routine does NOT write istep=xxx etc
! Write in ASCII or Binary format depending on zou_format
! ----------------------------------------------------------------------

  use common
  use indexQ
  implicit none

  integer iunit
  integer i,j,k,ialpha,ii

  if (zou_format.eq.'ASCII') then

     do i=0,imax
        do j=0,jmax
           do k=0,kmax
              do ii=1,nvar
                 ialpha=indQ(i,j,k,ii)
                 if (NDIM==3) then
                    write (iunit,'(4i5,e15.6)') i,j,k,ii,vec(ialpha)
                 else if (NDIM==2) then
                    write (iunit,'(3i5,e15.6)') i,j,ii,vec(ialpha)
                 else
                    write (iunit,'(2i5,e15.6)') i,ii,vec(ialpha)
                 endif
              enddo
           enddo
        enddo
     enddo

  else

     do i=0,imax
        do j=0,jmax
           do k=0,kmax
              do ii=1,nvar
                 ialpha=indQ(i,j,k,ii)
                 if (NDIM==3) then
                    write (iunit) i,j,k,ii,vec(ialpha)
                 else if (NDIM==2) then
                    write (iunit) i,j,ii,vec(ialpha)
                 else
                    write (iunit) i,ii,vec(ialpha)
                 endif
              enddo
           enddo
        enddo
     enddo

  endif

end subroutine snapshot

! ######################################################################

subroutine resid_snapshot2(iunit)          ! PIVOTS

  use common
  use indexQ

  integer iunit,i,j,k,ii,ialpha,ip,jQx,nunk(nvar),ivec(NDIM)
  double precision lhsI,rhsI,lhsav(nvar),rhsav(nvar),vecx
  logical found
  character(20) fmt

  lhsav=0
  rhsav=0
  nunk=0

  do i=0,imax
     do j=0,jmax
        do k=0,kmax
           do ii=1,nvar

              ialpha=indQ(i,j,k,ii)
              if (iunk(ialpha).eq.1) then
                 if (removeFixed=="YES") ialpha=veclookup(ialpha)

                 lhsI=0
                 rhsI=RR(ialpha)

                 do ip=irowst(ialpha),irowed(ialpha)
                    jQx=jQ(ip)
                    if (removeFixed=="NO") then
                       vecx=vec(jQx)
                    else
                       vecx=vecred(jQx)
                    endif

                    if (jQx.eq.ialpha) then
                       lhsI=Qval(ip)*vecx
                    else
                       rhsI=rhsI-Qval(ip)*vecx
                    endif
                 enddo

                 lhsav(ii)=lhsav(ii)+lhsI
                 rhsav(ii)=rhsav(ii)+rhsI
                 nunk(ii)=nunk(ii)+1
              endif
              
           enddo
        enddo
     enddo
  enddo

!  print *,lhsav
!  print *,rhsav
!  print *,nunk

  lhsav=lhsav/nunk
  rhsav=rhsav/nunk

  do i=0,imax
     do j=0,jmax
        do k=0,kmax
           do ii=1,nvar

              ialpha=indQ(i,j,k,ii)
              if (iunk(ialpha).eq.1) then
                 if (removeFixed=="YES") ialpha=veclookup(ialpha)

                 lhsI=0
                 rhsI=RR(ialpha)

                 found=.false.
                 do ip=irowst(ialpha),irowed(ialpha)
                    jQx=jQ(ip)
                    if (removeFixed=="NO") then
                       vecx=vec(jQx)
                    else
                       vecx=vecred(jQx)
                    endif

                    if (jQx.eq.ialpha) then
                       if (found) call prterr('Double pivot')
                       lhsI=Qval(ip)*vecx
                       found=.true.
                    else
                       rhsI=rhsI-Qval(ip)*vecx
                    endif
                 enddo

                 if (.not.found) call prterr('Pivot was not found in resid_snapshot2')
              else
                 lhsI=lhsav(ii)
                 rhsI=rhsav(ii)
              endif

              if (NDIM==3) then
                 ivec=[i,j,k]
                 fmt='(4i5,e15.6)'
              else if (NDIM==2) then
                 ivec=[i,j]
                 fmt='(3i5,e15.6)'
              else
                 ivec=[i]
                 fmt='(2i5,e15.6)'
              endif

              if (zou_format=='ASCII') then
                 write (iunit,fmt) ivec,3*(ii-1)+1,lhsI
                 write (iunit,fmt) ivec,3*(ii-1)+2,rhsI
                 write (iunit,fmt) ivec,3*(ii-1)+3,lhsI-rhsI
              else
                 write (iunit) ivec,3*(ii-1)+1,lhsI
                 write (iunit) ivec,3*(ii-1)+2,rhsI
                 write (iunit) ivec,3*(ii-1)+3,lhsI-rhsI
              endif

           enddo
        enddo
     enddo
  enddo

end subroutine resid_snapshot2

! ######################################################################

subroutine resid_snapshot(iunit)      ! UNKNOWNS

  use common
  use indexQ

  integer iunit,i,j,k,ii,ialpha,ip,jQx,nunk(nvar),ivec(NDIM)
  double precision lhsI,rhsI,lhsav(nvar),rhsav(nvar)
  character(20) fmt

  if (key_sim /= 2) then
     call prterr('RESID=UNKNOWNS should only be used for TRANSIENT')
  endif

  lhsav=0
  rhsav=0
  nunk=0

  do i=0,imax
     do j=0,jmax
        do k=0,kmax
           do ii=1,nvar

              ialpha=indQ(i,j,k,ii)
              if (iunk(ialpha).eq.1) then

                 lhsI=0
                 rhsI=RR(ialpha)
                 
                 do ip=irowst(ialpha),irowed(ialpha)
                    jQx=jQ(ip)
                    if (iunk(jQx)==1) then
                       lhsI=lhsI+Qval(ip)*vec(jQx)
                    else
                       rhsI=rhsI-Qval(ip)*vec(jQx)
                    endif
                 enddo

                 lhsav(ii)=lhsav(ii)+lhsI
                 rhsav(ii)=rhsav(ii)+rhsI
                 nunk(ii)=nunk(ii)+1
              endif

           enddo
        enddo
     enddo
  enddo

  lhsav=lhsav/nunk
  rhsav=rhsav/nunk

  do i=0,imax
     do j=0,jmax
        do k=0,kmax
           do ii=1,nvar

              ialpha=indQ(i,j,k,ii)

              if (iunk(ialpha).eq.1) then

                 lhsI=0
                 rhsI=RR(ialpha)

                 do ip=irowst(ialpha),irowed(ialpha)
                    jQx=jQ(ip)
                    if (iunk(jQx)==1) then
                       lhsI=lhsI+Qval(ip)*vec(jQx)
                    else
                       rhsI=rhsI-Qval(ip)*vec(jQx)
                    endif
                 enddo
              else
                 lhsI=lhsav(ii)
                 rhsI=rhsav(ii)
              endif

              if (NDIM==3) then
                 ivec=[i,j,k]
                 fmt='(4i5,e15.6)'
              else if (NDIM==2) then
                 ivec=[i,j]
                 fmt='(3i5,e15.6)'
              else
                 ivec=[i]
                 fmt='(2i5,e15.6)'
              endif

              ! eg V.lhs, V.rhs, V.resid, T.lhs, T.rhs, T.resid
              
              if (zou_format=='ASCII') then
                 write (iunit,fmt) ivec,3*(ii-1)+1,lhsI
                 write (iunit,fmt) ivec,3*(ii-1)+2,rhsI
                 write (iunit,fmt) ivec,3*(ii-1)+3,lhsI-rhsI
              else
                 write (iunit) ivec,3*(ii-1)+1,lhsI
                 write (iunit) ivec,3*(ii-1)+2,rhsI
                 write (iunit) ivec,3*(ii-1)+3,lhsI-rhsI
              endif
           
           enddo
        enddo
     enddo
  enddo

end subroutine resid_snapshot

! ######################################################################

subroutine pvsnapshot(iunit,tind)

! ----------------------------------------------------------------------
! Write a Paraview snapshot to unit iunit. Open and close the file
! if tind>0, write fname(tind).vtk. If tind=0 write to fname.vtk.
! ----------------------------------------------------------------------

  use common
  use indexQ
  implicit none

  integer iunit,tind,i,j,k,ii,ialpha,term
  character(30) temp
  character(1000) currfile

  if (tind.eq.0) then
     currfile=trim(fname)//'.vtk'
  else
     write (temp,*) tind
     currfile=trim(fname)//trim(adjustl(temp))//'.vtk'
  endif

  term=1
  call initfile(iunit,currfile,term,'unknown','ASCII')

  write (iunit,'(a)') '# vtk DataFile Version 3.0'
  write (iunit,'(a,i5,a,e13.5)') 'istep=',istep,' t=',istep*tstep
  write (iunit,'(a)') 'ASCII'
  write (iunit,'(a)') 'DATASET STRUCTURED_GRID'
  write (iunit,'(a,3i10)') 'DIMENSIONS',imax+1,jmax+1,kmax+1
  write (iunit,'(a,i10,a)') 'POINTS',nnod,' float'
  
  do i=0,imax
     do j=0,jmax
        do k=0,kmax
           write (iunit,'(3e13.5)') rnode(i,j,k,1),rnode(i,j,k,2),rnode(i,j,k,3)
        enddo
     enddo
  enddo
  
  write (iunit,'(a,i10)') 'POINT_DATA',nnod

  do ii=1,nvar
     
     write (temp,*) ii
     write (iunit,'(a)') 'SCALARS U'//trim(adjustl(temp))//' float 1'
     write (iunit,'(a)') 'LOOKUP_TABLE default'
     
     do i=0,imax
        do j=0,jmax
           do k=0,kmax
              ialpha=indQ(i,j,k,ii)
              write (iunit,'(e15.6)') vec(ialpha)
           enddo
        enddo
     enddo
  enddo
  
  close (iunit)
end subroutine pvsnapshot

! ######################################################################

function abbrev(string,n)

! ----------------------------------------------------------------------
! return string of length n corresponding to "string" perhaps abbreviated
! return either all of string, if len_trim(string)<=n
! OR abbreviated version with ... to represent abbreviation
! eg abracadabra (len_trim=11) with n=9 is abr...bra
! ----------------------------------------------------------------------

  character(*) string
  character(n) abbrev
  integer n,ls,lleft,lright

  ls=len_trim(string)

  if (ls.le.n) then
     abbrev=adjustl(trim(string))
  else
     lleft=(n-3)/2
     lright=n-lleft-3

     abbrev=adjustl(string(1:lleft)//'...'//string(ls-lright+1:ls))
  endif

end function abbrev

! ######################################################################

subroutine Cindex(string,ii,jj,kk,ll,matrix)

  use common, only : nvar,us,scoord,NDIM

  character(*) string
  integer ii,jj,kk,ll

  character(30) iis,jjs,kks,lls
  integer ivar,ic,ios
  logical found,matrix

  matrix=.false.

  read (string,*,iostat=ios) iis,jjs,kks,lls

  if (ios /= 0) then
     read (string,*,iostat=ios) iis,kks
     if (ios /= 0) call prterr('Incomprehensible C indices')
     matrix=.true.
     jj=0
     ll=0
  endif

! ----------------------------------------------------------------------
! ii
! ----------------------------------------------------------------------

  found=.false.
  do ivar=1,nvar
     if (iis.eq.us(ivar)) then
        ii=ivar
        found=.true.
     endif
  enddo

  if (.not.found) then
     read (iis,*,err=999) ii
  endif

! ----------------------------------------------------------------------
! jj
! ----------------------------------------------------------------------
  
  if (.not.matrix) then
     found=.false.
     do ic=1,NDIM
        if (jjs.eq.scoord(ic)) then
           jj=ic
           found=.true.
        endif
     enddo
     
     if (.not.found) then
        read (jjs,*,err=999) jj
     endif
  endif

! ----------------------------------------------------------------------
! kk
! ----------------------------------------------------------------------

  found=.false.
  do ivar=1,nvar
     if (kks.eq.us(ivar)) then
        kk=ivar
        found=.true.
     endif
  enddo

  if (.not.found) then
     read (kks,*,err=999) kk
  endif

! ----------------------------------------------------------------------
! ll
! ----------------------------------------------------------------------

  if (.not.matrix) then
     found=.false.
     do ic=1,NDIM
        if (lls.eq.scoord(ic)) then
           ll=ic
           found=.true.
        endif
     enddo
     
     if (.not.found) then
        read (lls,*,err=999) ll
     endif
  endif

  return

999 call prterr('Incomprehensible C index')
end subroutine Cindex

! ######################################################################

subroutine aindex(string,ii,kk)

  use common, only : nvar,us

  character(*) string
  integer ii,kk

  character(30) iis,kks
  integer ivar
  logical found

  read (string,*) iis,kks

! ----------------------------------------------------------------------
! ii
! ----------------------------------------------------------------------

  found=.false.
  do ivar=1,nvar
     if (iis.eq.us(ivar)) then
        ii=ivar
        found=.true.
     endif
  enddo

  if (.not.found) then
     read (iis,*,err=999) ii
  endif

! ----------------------------------------------------------------------
! kk
! ----------------------------------------------------------------------

  found=.false.
  do ivar=1,nvar
     if (kks.eq.us(ivar)) then
        kk=ivar
        found=.true.
     endif
  enddo

  if (.not.found) then
     read (kks,*,err=999) kk
  endif

  return

999 call prterr('Incomprehensible a index')
end subroutine Aindex

! ######################################################################

subroutine findex(string,ii)

  use common, only : nvar,us

  character(*) string
  integer ii,ivar
  logical found

  found=.false.
  do ivar=1,nvar
     if (string.eq.us(ivar)) then
        ii=ivar
        found=.true.
     endif
  enddo

  if (.not.found) then
     read (string,*,err=999) ii
  endif

  return

999 call prterr('Incomprehensible index')

end subroutine findex

end module iofile
