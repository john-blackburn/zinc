program zmeshcon

  use createmesh
  use global, only : minname,stemname

  integer ret,istatus,length,i

  call get_command_argument(1,stemname,length,istatus)

  if (istatus.gt.0) then
     print *,'Unable to retrieve command line argument!'
     print *,'Usage: zmesh file. This processes file.min and creates file.mtf'
     call exit(1)
  else if (istatus.eq.-1) then
     print *,'Filename too long, max 1000 characters'
     call exit(1)
  endif

  if (len_trim(stemname).eq.0) then
     print *,' Error: no input file specified!'
     call exit(1)
  endif

  do i=len_trim(stemname),1,-1
     if (stemname(i:i+3).eq.'.min') then
        stemname=stemname(:i-1)
        exit
     endif
  enddo

  minname=trim(stemname)//'.min'

! read in MIN. process. write out MTF
! allocates several variables in global

  ret=processMIN()

end program zmeshcon

! ######################################################################

subroutine createmessage(s)

  character(*) s
  print *,trim(s)

end subroutine createmessage

subroutine createerror(s)

  character(*) s
  print *,trim(s)
  stop

end subroutine createerror

subroutine snaperror(s)

  character(*) s
  print *,trim(s)

end subroutine snaperror
