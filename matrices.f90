! Part of Zinc FE package. Author: John Blackburn

module matrices

! ----------------------------------------------------------------------
! Set up and evalute c,a,f,q,g
! Used by zinc, qpr, zpp and scan (eventually zmesh?)
! Gives access to c,a,f,q,g possibly by calling user functions cfun etc
! Uses explicit DLL loading for these functions
! The pointers below should be initialised by calling usersetup
! CONDITIONAL COMPILATION: For static linking, define staticzinc (CPP)
! ----------------------------------------------------------------------

  ! define staticzinc for static linkage
  
  use common, only : EXPRLEN
  use evaluate, only : Tlist
  use ISO_C_BINDING, only : C_INTPTR_T

  implicit none

  save
!  private (would like to but...)

! These are the procedures
  public CCval,aaval,ffval,qqval,ggval
  public usersetup,userfree,lookup,lookind,addtable
  public getelinfo,getfaceinfo,fixnodes,readmatrix,wrtmatrix,freematrix

! Various data below is used. Hence all public for now

#ifdef STATICZINC

! Static interfaces

  interface
     function Cfun(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character label*(*)
       integer nvar,istep,imax,jmax,kmax
       double precision Cfun,x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function Cfun
     
     function afun(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character label*(*)
       integer nvar,istep,imax,jmax,kmax
       double precision afun,x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function afun

     function ffun(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character label*(*)  
       integer nvar,istep,imax,jmax,kmax
       double precision x,y,z,ffun,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function ffun
     
     function qfun(label,x,y,z,nx,ny,nz,ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character(*) label
       double precision qfun,x,y,z,nx,ny,nz,ur(nvar)
       integer nvar,istep,imax,jmax,kmax
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function qfun
     
     function gfun(label,x,y,z,nx,ny,nz,ur,nvar,istep, &
       ireg,iregup,rnode,vec,imax,jmax,kmax)

       character(*) label
       double precision gfun,x,y,z,nx,ny,nz,ur(nvar)
       integer nvar,istep,imax,jmax,kmax
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function gfun

     function BCfun(label,x,y,z,ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character(*) label
       double precision BCfun,x,y,z,ur(nvar)
       integer nvar,istep,imax,jmax,kmax
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function BCfun

     function scanfun(label,x,y,z,nx,ny,nz,ur,dur,nvar)
       character(*) label
       double precision scanfun,x,y,z,nx,ny,nz,ur(nvar),dur(nvar,*)
       integer nvar
     end function scanfun

! ----------------------------------------------------------------------

     subroutine dCfun_du(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax,dCdu)

       character(*) label
       integer nvar,istep,imax,jmax,kmax
       double precision dCdu(nvar),x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end subroutine dCfun_du

     subroutine dCfun_ddu(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax,dCddu)

       character(*) label
       integer nvar,istep,imax,jmax,kmax
       double precision dCddu(nvar,*),x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end subroutine dCfun_ddu

     subroutine dffun_du(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax,dfdu)

       character(*) label
       integer nvar,istep,imax,jmax,kmax
       double precision dfdu(nvar),x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end subroutine dffun_du

     subroutine dffun_ddu(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax,dfddu)

       character(*) label
       integer nvar,istep,imax,jmax,kmax
       double precision dfddu(nvar,*),x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end subroutine dffun_ddu
  end interface

#else

! interfaces for dynamic linking of DLL

  abstract interface
     function Cfunint(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character label*(*)
       integer nvar,istep,imax,jmax,kmax
       double precision Cfunint,x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function Cfunint
     
     function afunint(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character label*(*)
       integer nvar,istep,imax,jmax,kmax
       double precision afunint,x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function afunint

     function ffunint(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character label*(*)  
       integer nvar,istep,imax,jmax,kmax
       double precision x,y,z,ffunint,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function ffunint
     
     function qfunint(label,x,y,z,nx,ny,nz,ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character(*) label
       double precision qfunint,x,y,z,nx,ny,nz,ur(nvar)
       integer nvar,istep,imax,jmax,kmax
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function qfunint
     
     function gfunint(label,x,y,z,nx,ny,nz,ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)
       
       character(*) label
       double precision gfunint,x,y,z,nx,ny,nz,ur(nvar)
       integer nvar,istep,imax,jmax,kmax
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function gfunint

     function BCfunint(label,x,y,z,ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

       character(*) label
       double precision BCfunint,x,y,z,ur(nvar)
       integer nvar,istep,imax,jmax,kmax
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end function BCfunint

     function scanfunint(label,x,y,z,nx,ny,nz,ur,dur,nvar)

       character(*) label
       double precision scanfunint,x,y,z,nx,ny,nz,ur(nvar),dur(nvar,*)
       integer nvar
     end function scanfunint

! ----------------------------------------------------------------------

     subroutine dCfun_duint(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax,dCdu)

       character(*) label
       integer nvar,istep,imax,jmax,kmax
       double precision dCdu(nvar),x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end subroutine dCfun_duint

     subroutine dCfun_dduint(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax,dCddu)

       character(*) label
       integer nvar,istep,imax,jmax,kmax
       double precision dCddu(nvar,*),x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end subroutine dCfun_dduint

     subroutine dffun_duint(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax,dfdu)

       character(*) label
       integer nvar,istep,imax,jmax,kmax
       double precision dfdu(nvar),x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end subroutine dffun_duint

     subroutine dffun_dduint(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax,dfddu)

       character(*) label
       integer nvar,istep,imax,jmax,kmax
       double precision dfddu(nvar,*),x,y,z,ur(nvar),dur(nvar,*)
       integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
       double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)
     end subroutine dffun_dduint

!     function dCfun_duint(token,ur,dur)
!       character(*) token
!       double precision ur(:),dur(:,:),dCfun_duint(size(ur))
!     end function dCfun_duint
     
!     function dCfun_dduint(token,ur,dur)
!       character(*) token
!       double precision ur(:),dur(:,:),dCfun_dduint(size(ur),3)
!     end function dCfun_dduint
     
!     function dffun_duint(token,ur,dur)
!       character(*) token
!       double precision ur(:),dur(:,:),dffun_duint(size(ur))
!     end function dffun_duint
     
!     function dffun_dduint(token,ur,dur)
!       character(*) token
!       double precision ur(:),dur(:,:),dffun_dduint(size(ur),3)
!     end function dffun_dduint
  end interface

  procedure(cfunint), pointer :: cfun
  procedure(afunint), pointer :: afun
  procedure(ffunint), pointer :: ffun
  procedure(qfunint), pointer :: qfun
  procedure(gfunint), pointer :: gfun
  procedure(BCfunint), pointer :: BCfun
  procedure(scanfunint), pointer :: scanfun

  procedure(dCfun_duint),  pointer :: dCfun_du
  procedure(dCfun_dduint), pointer :: dCfun_ddu
  procedure(dffun_duint),  pointer :: dffun_du
  procedure(dffun_dduint), pointer :: dffun_ddu

  integer(C_INTPTR_T) :: hmodule=0

#endif

! ----------------------------------------------------------------------
! C, a, f, q, g storage
! ----------------------------------------------------------------------

  double precision, allocatable ::  CC(:,:,:,:,:)
  integer, allocatable          :: iCC(:,:,:,:,:)
  character(EXPRLEN), allocatable  :: sCC(:,:,:,:,:)
  type(Tlist), allocatable      :: tCC(:,:,:,:,:)
  
  double precision, allocatable ::  aa(:,:,:)
  integer, allocatable          :: iaa(:,:,:)
  character(EXPRLEN), allocatable  :: saa(:,:,:)
  type(Tlist), allocatable      :: taa(:,:,:)

  double precision, allocatable ::  ff(:,:)
  integer, allocatable          :: iff(:,:)
  character(EXPRLEN), allocatable  :: sff(:,:)
  type(Tlist), allocatable      :: tff(:,:)

  double precision, allocatable ::  qq(:,:,:)
  integer, allocatable          :: iqq(:,:,:)
  character(EXPRLEN), allocatable  :: sqq(:,:,:)
  type(Tlist), allocatable      :: tqq(:,:,:)

  double precision, allocatable ::  gg(:,:)
  integer, allocatable          :: igg(:,:)
  character(EXPRLEN), allocatable  :: sgg(:,:)
  type(Tlist), allocatable      :: tgg(:,:)

  integer, allocatable :: table(:,:)
  integer ntable

  logical, allocatable :: regCNL(:),regANL(:),regFNL(:),regQNL(:),regGNL(:),regfixNL(:)
  logical isQNL,isRNL,isfixNL,isNL

  integer, allocatable :: iregfix(:,:)
  double precision, allocatable ::  regfix(:,:)
  character(EXPRLEN), allocatable  :: sregfix(:,:)
  type(Tlist), allocatable      :: tregfix(:,:)

  double precision, allocatable :: init(:)
  integer, allocatable :: initi(:)
  character(EXPRLEN), allocatable :: inits(:)
  type(Tlist), allocatable :: initt(:)


  logical, private :: nvar_set, key_db_set, omega_set, key_sim_set, ndim_set
  logical, private :: nreg_set, restart_set, scale_set
  logical, private :: key_Q_set, nstride_set, labels_set, zou_format_set, rst_format_set
  logical, private :: nstep_set, tstep_set, nstrideNL_set, key_export_set

  logical, private :: idbmin_set,idbmax_set
  logical, private :: jdbmin_set,jdbmax_set
  logical, private :: kdbmin_set,kdbmax_set
  logical, private :: tol_set,itmax_set,residual_set,nodecheck_set,newton_set,removeFixed_set

contains

  subroutine usersetup(scansetup)

! ----------------------------------------------------------------------
! If any of (C,a,f,q,g,BC) have tokens or scansetup, load DLL
! set up appropriate functions
! if scansetup is true, set up scanfun also (for ZPP)
! Conditional compilation: if STATICZINC, this function is blank
! ----------------------------------------------------------------------

    use ISO_C_BINDING
    use common, only : fname,newton
    use iofile, only : prterr
    logical scansetup

#ifndef STATICZINC

    interface
       function LoadLibrary(lpFileName) bind(C,name='LoadLibraryA')
         use ISO_C_BINDING
         implicit none
         !GCC$ ATTRIBUTES STDCALL :: LoadLibrary
         integer(C_INTPTR_T) LoadLibrary
         character(kind=C_CHAR) lpFileName(*)
       end function LoadLibrary
       
       function GetProcAddress(hModule, lpProcName) bind(C,name='GetProcAddress')
         use ISO_C_BINDING
         implicit none
         !GCC$ ATTRIBUTES STDCALL :: GetProcAddress
         type(C_FUNPTR) GetProcAddress
         integer(C_INTPTR_T), value :: hModule
         character(kind=C_CHAR) lpProcName(*)
       end function GetProcAddress
    end interface

    type(C_FUNPTR) FunAddress

    logical CCtoken,aatoken,fftoken,qqtoken,ggtoken,BCtoken

    CCtoken=any(iCC.eq.3)
    aatoken=any(iaa.eq.3)
    fftoken=any(iff.eq.3)
    qqtoken=any(iqq.eq.3)
    ggtoken=any(igg.eq.3)
    BCtoken=any(iregfix.eq.4)    ! Note "4"

    if (CCtoken.or.aatoken.or.fftoken.or.qqtoken.or.ggtoken.or.BCtoken.or.scansetup) then
       
       print *,'Loading DLL file: ',trim(fname)//'.dll'
       
       hmodule=LoadLibrary(trim(fname)//'.dll'//achar(0))
       
       if (hmodule.eq.0) call prterr('DLL not found')
       
       if (CCtoken) then
          FunAddress=GetProcAddress(hmodule,'cfun'//achar(0))
          if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find cfun')
          call C_F_PROCPOINTER(FunAddress,cfun)

          if (newton=='YES') then
             FunAddress=GetProcAddress(hmodule,'dcfun_du'//achar(0))
             if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find dCfun_du')
             call C_F_PROCPOINTER(FunAddress,dcfun_du)

             FunAddress=GetProcAddress(hmodule,'dcfun_ddu'//achar(0))
             if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find dCfun_ddu')
             call C_F_PROCPOINTER(FunAddress,dcfun_ddu)
          endif
       endif
       
       if (aatoken) then
          FunAddress=GetProcAddress(hmodule,'afun'//achar(0))
          if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find afun')
          call C_F_PROCPOINTER(FunAddress,afun)
       endif

       if (fftoken) then
          FunAddress=GetProcAddress(hmodule,'ffun'//achar(0))
          if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find ffun')
          call C_F_PROCPOINTER(FunAddress,ffun)

          if (newton=='YES') then
             FunAddress=GetProcAddress(hmodule,'dffun_du'//achar(0))
             if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find dffun_du')
             call C_F_PROCPOINTER(FunAddress,dffun_du)

             FunAddress=GetProcAddress(hmodule,'dffun_ddu'//achar(0))
             if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find dffun_ddu')
             call C_F_PROCPOINTER(FunAddress,dffun_ddu)
          endif
       endif

       if (qqtoken) then
          FunAddress=GetProcAddress(hmodule,'qfun'//achar(0))
          if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find qfun')
          call C_F_PROCPOINTER(FunAddress,qfun)
       endif

       if (ggtoken) then
          FunAddress=GetProcAddress(hmodule,'gfun'//achar(0))
          if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find gfun')
          call C_F_PROCPOINTER(FunAddress,gfun)
       endif
       
       if (BCtoken) then                                
          FunAddress=GetProcAddress(hmodule,'bcfun'//achar(0))
          if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find BCfun')
          call C_F_PROCPOINTER(FunAddress,BCfun)
       endif
       
       if (scansetup) then
          FunAddress=GetProcAddress(hmodule,'scanfun'//achar(0))
          if (.NOT.C_ASSOCIATED(FunAddress)) call prterr('Could not find scanfun')
          call C_F_PROCPOINTER(FunAddress,scanfun)
       endif
    endif

#endif 

  end subroutine usersetup

! ######################################################################

  subroutine userfree

#ifndef STATICZINC
    
    interface
       function FreeLibrary(hLibModule) bind(C,name='FreeLibrary')
         use ISO_C_BINDING
         implicit none
         !GCC$ ATTRIBUTES STDCALL :: FreeLibrary
         integer(C_INT) FreeLibrary
         integer(C_INTPTR_T), value :: hLibModule
       end function FreeLibrary
    end interface

    integer ret

    if (hmodule.ne.0) then
       ret=FreeLibrary(hmodule)
       if (ret.eq.0) print *,'Warning: failed to release DLL library'
    endif

#endif

  end subroutine userfree

! ######################################################################

function CCval(irege,ii,jj,kk,ll,ur,dur,rc,iCCx)

! ----------------------------------------------------------------------
! Get C(ii,jj,kk,ll) for reg irege
! In cases 2 and 3, this will depend on ur, dur, rc
! These can be calculated using getelinfo
! ----------------------------------------------------------------------

  use evaluate
  use common

  integer irege,ii,jj,kk,ll
  double precision CCval,ur(nvar),dur(nvar,NDIM),rc(NDIM),x,y,z
  integer iCCx,ivar,j
                    
  iCCx=iCC(irege,ii,jj,kk,ll)

               x=rc(1)
  if (NDIM>=2) y=rc(2)
  if (NDIM==3) z=rc(3)

  if (iCCx.eq.0) then
     CCval=CC(irege,ii,jj,kk,ll)
  else if (iCCx.eq.1) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)
     
!     call evalexpr(sCCx,CCval)
     call evaltoken(tCC(irege,ii,jj,kk,ll),CCval)

  else if (iCCx.eq.2) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)
     
     do ivar=1,nvar
        call defparam(us(ivar),ur(ivar))
        do j=1,NDIM
           call defparam(dus(ivar,j),dur(ivar,j))
        enddo
     enddo
     
!     call evalexpr(sCCx,CCval)
     call evaltoken(tCC(irege,ii,jj,kk,ll),CCval)

  else if (iCCx.eq.3) then
     CCval= Cfun(sCC(irege,ii,jj,kk,ll),x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)
  endif
  
end function CCval

! ######################################################################

function aaval(irege,ii,kk,ur,dur,rc,iaax)

  use evaluate
  use common

  integer irege,ii,kk
  double precision aaval,ur(nvar),dur(nvar,NDIM),rc(NDIM),x,y,z
  integer iaax,ivar,j

  iaax=iaa(irege,ii,kk)

               x=rc(1)
  if (NDIM>=2) y=rc(2)
  if (NDIM==3) z=rc(3)

  if (iaax.eq.0) then
     aaval=aa(irege,ii,kk)
  else if (iaax.eq.1) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)
     
!     call evalexpr(saax,aaval)
     call evaltoken(taa(irege,ii,kk),aaval)
     
  else if (iaax.eq.2) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)
     
     do ivar=1,nvar
        call defparam(us(ivar),ur(ivar))
        do j=1,3
           call defparam(dus(ivar,j),dur(ivar,j))
        enddo
     enddo
     
!     call evalexpr(saax,aaval)
     call evaltoken(taa(irege,ii,kk),aaval)

  else if (iaax.eq.3) then
     aaval= afun(saa(irege,ii,kk),x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)
  endif

end function aaval

! ######################################################################

function ffval(irege,ii,ur,dur,rc)

  use evaluate
  use common

  integer irege,ii
  double precision ffval,ur(nvar),dur(nvar,NDIM),rc(NDIM),x,y,z
  integer iffx,ivar,j

  iffx=iff(irege,ii)

               x=rc(1)
  if (NDIM>=2) y=rc(2)
  if (NDIM==3) z=rc(3)

  if (iffx.eq.0) then
     ffval=ff(irege,ii)
  else if (iffx.eq.1) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)

!     call evalexpr(sffx,ffval)
     call evaltoken(tff(irege,ii),ffval)
     
  else if (iffx.eq.2) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)
     
     do ivar=1,nvar
        call defparam(us(ivar),ur(ivar))
        do j=1,3
           call defparam(dus(ivar,j),dur(ivar,j))
        enddo
     enddo
     
!     call evalexpr(sffx,ffval)
     call evaltoken(tff(irege,ii),ffval)
     
  else if (iffx.eq.3) then
     ffval= ffun(sff(irege,ii),x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)
  endif

end function ffval

! ######################################################################

function qqval(irege1,irege2,ii,kk,ur,rc,norm,iqqx)

  use evaluate
  use common

  integer irege1,irege2,ii,kk
  double precision qqval,ur(nvar),rc(NDIM),norm(3),x,y,z
  integer iqqx,ivar,ind
                    
  ind=lookup(irege1,irege2)

  if (ind.eq.0) then
     iqqx=-1   ! not found
     qqval=0
     return
  endif

  iqqx=iqq(ind,ii,kk)

               x=rc(1)
  if (NDIM>=2) y=rc(2)
  if (NDIM==3) z=rc(3)

  if (iqqx.eq.0) then
     qqval=qq(ind,ii,kk)
  else if (iqqx.eq.1) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)

     call defparam('NX',norm(1))
     call defparam('NY',norm(2))
     call defparam('NZ',norm(3))

!     call evalexpr(sqqx,qqval)
     call evaltoken(tqq(ind,ii,kk),qqval)

  else if (iqqx.eq.2) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)

     call defparam('NX',norm(1))
     call defparam('NY',norm(2))
     call defparam('NZ',norm(3))
     
     do ivar=1,nvar
        call defparam(us(ivar),ur(ivar))
     enddo
     
!     call evalexpr(sqqx,qqval)
     call evaltoken(tqq(ind,ii,kk),qqval)

  else if (iqqx.eq.3) then
     qqval=qfun(sqq(ind,ii,kk),x,y,z,norm(1),norm(2),norm(3),ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)
  endif
  
end function qqval

! ######################################################################

function ggval(irege1,irege2,ii,ur,rc,norm)

  use evaluate
  use common

  implicit none

  integer irege1,irege2,ii
  double precision ggval,ur(nvar),rc(NDIM),norm(3),x,y,z
  integer iggx,ivar,ind
                    
  ind=lookup(irege1,irege2)

  if (ind.eq.0) then
     ggval=0
     return
  endif

  iggx=igg(ind,ii)

               x=rc(1)
  if (NDIM>=2) y=rc(2)
  if (NDIM==3) z=rc(3)

  if (iggx.eq.0) then
     ggval=gg(ind,ii)
  else if (iggx.eq.1) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)
     
     call defparam('NX',norm(1))
     call defparam('NY',norm(2))
     call defparam('NZ',norm(3))

!     call evalexpr(sggx,ggval)
     call evaltoken(tgg(ind,ii),ggval)

  else if (iggx.eq.2) then
                  call defparam('X',x)
     if (NDIM>=2) call defparam('Y',y)
     if (NDIM==3) call defparam('Z',z)
     call defparam('ISTEP',istep)
     
     call defparam('NX',norm(1))
     call defparam('NY',norm(2))
     call defparam('NZ',norm(3))

     do ivar=1,nvar
        call defparam(us(ivar),ur(ivar))
     enddo
     
!     call evalexpr(sggx,ggval)
     call evaltoken(tgg(ind,ii),ggval)

  else if (iggx.eq.3) then
     ggval=gfun(sgg(ind,ii),x,y,z,norm(1),norm(2),norm(3),ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)
  endif
  
end function ggval

! ######################################################################

function lookup(ireg1,ireg2)

! given region numbers of elements on either side of face. Get index number of face

  integer lookup,ireg1,ireg2,i

  do i=1,ntable
     if ((ireg1.eq.table(i,1).and.ireg2.eq.table(i,2)).or.&
         (ireg1.eq.table(i,2).and.ireg2.eq.table(i,1))) then
        lookup=i
        return
     endif
  enddo

  lookup=0
end function lookup

! ######################################################################

function lookind(string)

! convert user input to integer. user can enter [xyz]max/min

  use iofile, only : prterr

  character(*) string
  integer lookind

  if (string.eq.'XMIN') then
     lookind=-1
  else if (string.eq.'XMAX') then
     lookind=-2
  else if (string.eq.'YMIN') then
     lookind=-3
  else if (string.eq.'YMAX') then
     lookind=-4
  else if (string.eq.'ZMIN') then
     lookind=-5
  else if (string.eq.'ZMAX') then
     lookind=-6
  else
     read (string,*,err=999) lookind
  endif

  return

999 call prterr('Incomprehensible surface specification')

end function lookind

! ######################################################################

subroutine addtable(ireg1,ireg2)

  integer ireg1,ireg2,i

  do i=1,ntable
     if ((ireg1.eq.table(i,1).and.ireg2.eq.table(i,2)).or.&
         (ireg1.eq.table(i,2).and.ireg2.eq.table(i,1))) then
!        call prterr('Duplicate surface specification')
        return
     endif
  enddo

  ntable=ntable+1
  table(ntable,1)=ireg1
  table(ntable,2)=ireg2
end subroutine addtable

! ######################################################################

subroutine getelinfo(itabsrt,rtabsrt,ur,dur)

! ----------------------------------------------------------------------
! Only called for non-linear systems
! note that element centre is not here since this may be needed
! for linear problems which depend on space
! ----------------------------------------------------------------------

  use common, only : vec,nvar,zero,NDIM
  use geom
  use shape
  use indexQ

  integer itabsrt(2**NDIM,NDIM)
  double precision rtabsrt(2**NDIM,NDIM),ur(nvar),dur(nvar,NDIM)
  
  integer ii,jj,l,inode,jnode,knode
  double precision ul(2**NDIM,nvar),jac(NDIM,NDIM),invjac(NDIM,NDIM),jacdet,dN(NDIM)

  if (NDIM==3) then
     do l=1,8
        inode=itabsrt(l,1)
        jnode=itabsrt(l,2)
        knode=itabsrt(l,3)
        
        do ii=1,nvar
           ul(l,ii)=vec(indQ(inode,jnode,knode,ii))
        enddo
     enddo
  else if (NDIM==2) then
     do l=1,4
        inode=itabsrt(l,1)
        jnode=itabsrt(l,2)
        
        do ii=1,nvar
           ul(l,ii)=vec(indQ(inode,jnode,0,ii))
        enddo
     enddo
  else
     do l=1,2
        inode=itabsrt(l,1)
        
        do ii=1,nvar
           ul(l,ii)=vec(indQ(inode,0,0,ii))
        enddo
     enddo
  endif

  call jacobian(zero,zero,zero,rtabsrt,jac,invjac,jacdet)
           
  ur=0
  dur=0
           
  do l=1,2**NDIM
     call getdN(l,zero,zero,zero,invjac,dN)
     
     do ii=1,nvar
        ur(ii)=ur(ii)+ul(l,ii)*Nl(l,zero,zero,zero)
     enddo
                  
     do ii=1,nvar
        do jj=1,NDIM
           dur(ii,jj)=dur(ii,jj)+ul(l,ii)*dN(jj)
        enddo
     enddo
  enddo

end subroutine getelinfo

! ######################################################################

subroutine getfaceinfo(irege1,irege2,centre1,centre2,itabsrt,rtabsrt,ur,rc,norm)

! ----------------------------------------------------------------------
! Called for all systems (with Neuman boundaries) including
! linear systems. itabsrt, rtabsrt are (8,3) but only first (4,3) is used
! Bit unnecessary in linear systems but faces not called for every element
! ----------------------------------------------------------------------

  use common, only : vec,nvar,zero
  use indexQ
  use shape
  use util
  use iofile, only : prterr

  integer irege1,irege2,itabsrt(8,3)
  double precision centre1(3),centre2(3),rtabsrt(8,3),ur(nvar),rc(3),norm(3)

  integer l,inode,jnode,knode,ii,ind
  double precision ul(4,nvar)
  double precision norm1(3),norm2(3),norm3(3),norm4(3),test(3),vec1(3),vec2(3)
  
  rc(1)=sum(rtabsrt(1:4,1))/4
  rc(2)=sum(rtabsrt(1:4,2))/4
  rc(3)=sum(rtabsrt(1:4,3))/4
  
  do l=1,4
     inode=itabsrt(l,1)
     jnode=itabsrt(l,2)
     knode=itabsrt(l,3)
     
     do ii=1,nvar
        ul(l,ii)=vec(indQ(inode,jnode,knode,ii))
     enddo
  enddo

  ur=0
  do l=1,4
     do ii=1,nvar
        ur(ii)=ur(ii)+ul(l,ii)*N2d(l,zero,zero)
     enddo
  enddo

  ind=lookup(irege1,irege2)

  if (ind.eq.0) then
     call prterr('getfaceinfo: illegal face')
  endif

  vec1=rtabsrt(2,:)-rtabsrt(1,:)
  vec2=rtabsrt(4,:)-rtabsrt(1,:)
  norm1=cross(vec1,vec2)
  
  vec1=rtabsrt(3,:)-rtabsrt(2,:)
  vec2=rtabsrt(1,:)-rtabsrt(2,:)
  norm2=cross(vec1,vec2)
  
  vec1=rtabsrt(4,:)-rtabsrt(3,:)
  vec2=rtabsrt(2,:)-rtabsrt(3,:)
  norm3=cross(vec1,vec2)
  
  vec1=rtabsrt(1,:)-rtabsrt(4,:)
  vec2=rtabsrt(3,:)-rtabsrt(4,:)
  norm4=cross(vec1,vec2)
  
  if (dot_product(norm1,norm2).lt.0.or.dot_product(norm1,norm3).lt.0.or.dot_product(norm1,norm4).lt.0.or.&
      dot_product(norm2,norm3).lt.0.or.dot_product(norm2,norm4).lt.0.or.dot_product(norm3,norm4).lt.0) then
     call prterr('incompatible normal vectors')
  endif
  
  norm=(norm1+norm2+norm3+norm4)/4
  norm=norm/mag(norm)
  
  if (irege2.eq.table(ind,2)) then
     test=centre2-centre1
  else
     test=centre1-centre2
  endif
  
  if (dot_product(norm,test).lt.0) norm=-norm

end subroutine getfaceinfo

! ######################################################################

subroutine fixnodes

! ----------------------------------------------------------------------
! Set nodes which are fixed to rightful values. If the fixed nodes
! depend on initial state of others nodes these should be set before
! calling. Note that this does not mark the nodes fixed in iunk
! ----------------------------------------------------------------------

  use common
  use evaluate
  use indexQ
  
  integer ii,i,j,k,ireg1,jj
  double precision ur(nvar),x,y,z,val
  
  do ii=1,nvar
     do i=0,imax
        do j=0,jmax
           do k=0,kmax
              
              x=rnode(i,j,k,1)
              if (NDIM>=2) y=rnode(i,j,k,2)
              if (NDIM==3) z=rnode(i,j,k,3)
              
              ireg1=ireg(i,j,k)
              if (iregfix(ireg1,ii).eq.1) then
                 vec(indQ(i,j,k,ii))=regfix(ireg1,ii)
                 
              else if (iregfix(ireg1,ii).eq.2) then
                 call defparam('X',x)
                 if (NDIM>=2) call defparam('Y',y)
                 if (NDIM==3) call defparam('Z',z)
                 call defparam('ISTEP',istep)
                 
                 call evaltoken(tregfix(ireg1,ii),val)

                 if (ierr.ne.0) then
                    print *,'ERROR: Could not resolve initial expression for variable ',ii
                    call exit(1)
                 endif
                 
                 vec(indQ(i,j,k,ii))=val
              else if (iregfix(ireg1,ii).eq.3) then
                 call defparam('X',x)
                 if (NDIM>=2) call defparam('Y',y)
                 if (NDIM==3) call defparam('Z',z)
                 call defparam('ISTEP',istep)
                 
                 do jj=1,nvar
                    ur(jj)=vec(indQ(i,j,k,jj))
                 enddo
                 
                 do jj=1,nvar
                    call defparam(us(jj),ur(jj))
                 enddo
                 
                 call evaltoken(tregfix(ireg1,ii),val)
                 vec(indQ(i,j,k,ii))=val
                 
              else if (iregfix(ireg1,ii).eq.4) then

                 do jj=1,nvar
                    ur(jj)=vec(indQ(i,j,k,jj))
                 enddo

                 vec(indQ(i,j,k,ii))=BCfun(sregfix(ireg1,ii),x,y,z,ur,nvar,istep, &
                      ireg,iregup,rnode,vec,imax,jmax,kmax)
              endif
           enddo
        enddo
     enddo
  enddo
end subroutine fixnodes

! ######################################################################

subroutine readmatrix

! ----------------------------------------------------------------------
! Read ZIN and CON files. Allocate and fill matrix related arrays
! do sanity checks. Set up us,dus,scoord
! also sets node fix info (Dirichlet boundaries) and init
! Set all global ZIN variables
! Set up constants and u1(xyz) etc [or user names] using DEFPARAM
! opens own files unit 1 (which must be available)
! ----------------------------------------------------------------------

  use common
  use evaluate
  use iofile

  integer, allocatable :: isetel(:),isetnd(:)
  integer allocerr(30),term,i,ios,ivar,j,ireg1,ireg2,ii,jj,kk,ll,irg,ind
  integer ivals,ival,nvals,itable,k,iiold,ist,ied
  character(1000) lefts,rights,vals,line
  character cdum*20,cvar*20,keyword*20
  character(20) temp,sreg1,sreg2
  type(Tlist) list
  logical nodeinit,found,submatrix,endok
  
  double precision val

  allocate(list%token(numtok))  ! local

! ----------------------------------------------------------------------
! Pre-read .zin file purely to allow matrix allocation
! we just need NDIM, nvar, nreg and ntable. Also nregel and nregnd
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.zin',term,'old','ASCII')

  nvar_set=.false.; key_db_set=.false.; omega_set=.false.
  nreg_set=.false.; restart_set=.false.; scale_set=.false.
  key_Q_set=.false.; nstride_set=.false.
  nstep_set=.false.; tstep_set=.false.; nstrideNL_set=.false.
  key_export_set=.false. ; labels_set=.false.; zou_format_set=.false.; rst_format_set=.false.
  tol_set=.false.; itmax_set=.false.; residual_set=.false.; nodecheck_set=.false.; newton_set=.false.
  nvar_set=.false.; ndim_set=.false.

  idbmin_set=.false.
  idbmax_set=.false.
  jdbmin_set=.false.
  jdbmax_set=.false.
  kdbmin_set=.false.
  kdbmax_set=.false.

  ntable=0
  nregel=0
  nregnd=0
  do
     call getline(1,ios,line)
!     print *,trim(line)

     if (ios.lt.0) exit
     line=adjustl(line)

     if (index(line,'=').ne.0) then
        
        i=spliteq(line,lefts,rights)
        
        if (i.ne.0) then
           
           if (lefts.eq.'NVAR') then
              read (rights,*,iostat=ios) nvar
              if (ios.ne.0) call prterr('Could not understand "nvar" in .zin file')
              nvar_set=.true.
           else if (lefts.eq.'NDIM') then
              read (rights,*,iostat=ios) NDIM
              if (ios.ne.0) call prterr('Could not understand "ndim" in .zin file')
              ndim_set=.true.
           endif
        endif
     else if (line(1:7)=='SURFACE') then
        ntable=ntable+1
     else if (line(1:6)=="REGION") then
        read (line,*,err=911) cdum,irg,keyword
        if (keyword=="ELEMENTS") then
           nregel=max(nregel,irg)
        else
           nregnd=max(nregnd,irg)
        endif
     endif
  enddo

  close (1)

  if (.not.nvar_set) then
     call prterr('nvar not found in .zin file')
  endif

  if (.not.ndim_set) then
     print *,'ndim not found in .zin file, assuming 3D system'
     ndim=3
  endif

  if (NDIM<1.or.NDIM>3) then
     call prterr('NDIM must be 1, 2 or 3')
  endif

  print *,nregel,' Element regions found'
  print *,nregnd,' Node regions found'

! ----------------------------------------------------------------------
! Allocate arrays iset,C,aa,ff. allocate the table for surfaces.
! ----------------------------------------------------------------------

  allocate (isetel(nregel),isetnd(nregnd),stat=allocerr(1))  ! locals
  
  allocate ( CC(nregel,nvar,NDIM,nvar,NDIM),stat=allocerr(2))
  allocate (iCC(nregel,nvar,NDIM,nvar,NDIM),stat=allocerr(3))
  allocate (sCC(nregel,nvar,NDIM,nvar,NDIM),stat=allocerr(4))
  allocate (tCC(nregel,nvar,NDIM,nvar,NDIM),stat=allocerr(5))

  allocate ( aa(nregel,nvar,nvar),stat=allocerr(6))
  allocate (iaa(nregel,nvar,nvar),stat=allocerr(7))
  allocate (saa(nregel,nvar,nvar),stat=allocerr(8))
  allocate (taa(nregel,nvar,nvar),stat=allocerr(9))
  
  allocate ( ff(nregel,nvar),stat=allocerr(10))
  allocate (iff(nregel,nvar),stat=allocerr(11))
  allocate (sff(nregel,nvar),stat=allocerr(12))
  allocate (tff(nregel,nvar),stat=allocerr(13))

  allocate (iregfix(nregnd,nvar),regfix(nregnd,nvar),&
       sregfix(nregnd,nvar),tregfix(nregnd,nvar),stat=allocerr(14))
  
  allocate (init(nvar),initi(nvar),inits(nvar),initt(nvar),stat=allocerr(15))

  allocate (us(nvar),dus(nvar,NDIM),stat=allocerr(16))

  allocate ( qq(ntable,nvar,nvar),stat=allocerr(17))
  allocate (iqq(ntable,nvar,nvar),stat=allocerr(18))
  allocate (sqq(ntable,nvar,nvar),stat=allocerr(19))
  allocate (tqq(ntable,nvar,nvar),stat=allocerr(20))

  allocate ( gg(ntable,nvar),stat=allocerr(21))
  allocate (igg(ntable,nvar),stat=allocerr(22))
  allocate (sgg(ntable,nvar),stat=allocerr(23))
  allocate (tgg(ntable,nvar),stat=allocerr(24))

  allocate (table(ntable,2),stat=allocerr(25))

  allocate (regQNL(ntable),stat=allocerr(26))
  allocate (regGNL(ntable),stat=allocerr(27))

  do i=1,27
     if (allocerr(i).ne.0) then
        write (temp,*) i
        call prterr('Failed to allocate array:'//trim(temp))
     endif
  enddo

! ----------------------------------------------------------------------
!  Read .zin file again to find LABELS and set up the surface table
!  Also discover RESTART and RST_FORMAT
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.zin',term,'old','ASCII')

  ntable=0
  restart='NONE'
  rst_format='BINARY'
  do
     call getline(1,ios,line)

     if (ios.lt.0) exit

     if (index(line,'=').ne.0) then

        i=spliteq(line,lefts,rights)

        if (i.ne.0) then
           
           if (lefts.eq.'LABELS') then
              read (rights,*,iostat=ios) (us(ii),ii=1,nvar)
              if (ios.ne.0) call prterr('Could not understand "labels" in .zin file')
              labels_set=.true.
           else if (lefts.eq.'RESTART') then
              restart=adjustl(rights)
              restart_set=.true.
           else if (lefts.eq.'RST_FORMAT') then
              if (rights /= 'ASCII' .and. rights /= 'BINARY') then
                 call prterr('RST_FORMAT must be "ASCII" or "BINARY"')
              endif
              
              rst_format=trim(rights)
              rst_format_set=.true.
           endif
        endif
     else if (index(line,'SURFACE').ne.0) then
        do i=1,len_trim(line)
           if (line(i:i).eq.'-') line(i:i)=' '
        enddo

        read (line,*,err=911) cdum,sreg1,sreg2
        ireg1=lookind(sreg1)
        ireg2=lookind(sreg2)

        call addtable(ireg1,ireg2)
     endif

  enddo

  close (1)

! ----------------------------------------------------------------------
! Pre-read restart file, if specified, to find nvar_rst
! ----------------------------------------------------------------------

  if (restart /= 'NONE') then
     term=1
     call initfile(1,restart,term,'old',rst_format)

     if (rst_format=='ASCII') then ! Ignore strapline
        read (1,*)
        read (1,*)
        read (1,*)
     else
        read (1) 
        read (1)
     endif

     iiold=0
     nvar_rst=0
     do
        if (rst_format=='ASCII') then
           read (1,*,iostat=ios) i,j,k,ii
        else
           read (1,iostat=ios) i,j,k,ii
        endif

!        print *,ii,iiold

        if (ios.ne.0) call prterr('Problem reading .rst file')
        
        nvar_rst=max(nvar_rst,ii)
        
        if (ii<=iiold) exit
        iiold=ii
     enddo

     print *,nvar_rst,' variables found in restart file'
  endif

! ----------------------------------------------------------------------
! prepare strings for vars and derivatives u1, u2, u1x etc
! ----------------------------------------------------------------------

  allocate (rs(nvar_rst),stat=allocerr(1))
  if (allocerr(1)/=0) call prterr('Could not allocate rs')

  scoord=['X','Y','Z']

  if (.not.labels_set) then
     do ivar=1,nvar
        write (temp,*) ivar
        temp=adjustl(temp)
        us(ivar)='U'//temp
     enddo
  endif

  if (restart /= 'NONE') then
     do ivar=1,nvar_rst
        write (temp,*) ivar
        rs(ivar)='R'//adjustl(temp)
     enddo
  endif

  do ivar=1,nvar
     do j=1,NDIM
        dus(ivar,j)=trim(us(ivar))//scoord(j)
     enddo
  enddo

! ----------------------------------------------------------------------
! set region defaults, nodes not fixed and elements in each reg 
! have C=0, a=0, f=0. Unspecified region components are linear
! with zero fixed value (matrix operations below)
! ----------------------------------------------------------------------

  isetel=0; isetnd=0
  regfix=0
  iregfix=0
  sregfix='-'

  CC=0
  iCC=0
  sCC='-'
  
  aa=0
  iaa=0
  saa='-'
  
  ff=0
  iff=0
  sff='-'

  qq=0
  iqq=0
  sqq='-'

  gg=0
  igg=0
  sgg='-'

  init=0
  initi=-99
  inits='-'
  
! ----------------------------------------------------------------------
! Now read in .con file
! ----------------------------------------------------------------------

  term=0
  call initfile(1,trim(fname)//'.con',term,'old','ASCII')

  if (term.eq.-1) then
     print *,'Warning: constants file not found - there will be no constants!'
  else

     do
        call getline(1,ios,line)

        if (ios.lt.0) exit

        i=spliteq(line,lefts,rights)
        
        if (i.eq.0) then
           call prterr('Incomprehensible line in .con file')
        else
           
           if (lefts.eq.'I'.or.lefts.eq.'PI'.or.lefts.eq.'REGION') then
              call prterr('Reserved variable names: I, PI, REGION')
           endif
           
           do ivar=1,nvar
              if (lefts.eq.us(ivar)) then
                 call prterr(trim(us(ivar))//' may not be reset')
              endif
              do j=1,NDIM
                 
                 if (lefts.eq.dus(ivar,j)) then
                    call prterr(trim(dus(ivar,j))//' may not be reset')
                 endif
                 
              enddo
           enddo
           
           do j=1,NDIM
              if (lefts.eq.scoord(j)) then
                 call prterr('x,y,z may not be reset')
              endif
           enddo

           if (lefts.eq.'NX'.or.lefts.eq.'NY'.or.lefts.eq.'NZ') then
              call prterr('NX, NY, NZ may not be reset')
           endif

           call defparam(lefts,rights)
           
           if (ierr.ne.0) call prterr('Error parsing string')
        endif

     enddo
  endif

  close (1)

!  call listvar

! ----------------------------------------------------------------------
! Add u1, u2, u1x,..., x, y, z to variable list and set to zero
! This is needed for parsval which will fail if expression does not
! make sense
! ----------------------------------------------------------------------

  call defparam('X',0)
  call defparam('Y',0)
  call defparam('Z',0)

  call defparam('NX',0)
  call defparam('NY',0)
  call defparam('NZ',0)

  do ivar=1,nvar
     call defparam(us(ivar),0)
     do j=1,NDIM
        call defparam(dus(ivar,j),0)
     enddo
  enddo

  if (restart /= 'NONE') then
     do ivar=1,nvar_rst
        call defparam(rs(ivar),0)
     enddo
  endif

  call defparam('ISTEP',0)

! ----------------------------------------------------------------------
! Input data fname.zin. Note regfix only important if iregfix=1
! ----------------------------------------------------------------------

  term=1
  call initfile(1,trim(fname)//'.zin',term,'old','ASCII')

! ----------------------------------------------------------------------
! Now, read control variables and region data
! set defaults
! ----------------------------------------------------------------------

  key_db=0     ! default values
  nstep=10
  omega=0.1
  scale=1
  key_Q=1
  nstride=10
  nstrideNL=10
  key_export=0
  itmax=1000
  tol=1e-6
  zou_format='BINARY'
  residual='NONE'
  nodecheck='NO'
  newton='NO'
  removeFixed='YES'

!  print *,'starting zin loop'

  endok=.false.
  do
     call getline(1,ios,line)

     line=adjustl(line)
     print *,trim(line)

     if (ios.ne.0) exit

     if (index(line,'=').ne.0) then

! ----------------------------------------------------------------------
! Control parameters, param = value
! ----------------------------------------------------------------------

        i=spliteq(line,lefts,rights)

        if (lefts.eq.'NVAR') then
           ! already done
        else if (lefts.eq.'NDIM') then
           ! already done
        else if (lefts.eq.'KEY_DB') then
           read (rights,*,err=911) key_db
           key_db_set=.true.

        else if (lefts.eq.'OMEGA') then
           read (rights,*,err=911) omega
           omega_set=.true.

        else if (lefts.eq.'NREG') then
           ! ignore

        else if (lefts.eq.'RESTART') then
           ! already done

        else if (lefts.eq.'SCALE') then
           read (rights,*,err=911) scale
           scale_set=.true.

        else if (lefts.eq.'KEY_Q') then
           read (rights,*,err=911) key_Q
           key_Q_set=.true.

        else if (lefts.eq.'SOLVER') then

           if (rights=='SOR_SYM') then
              key_sim=0
           else if (rights=='SOR') then
              key_sim=1
           else if (rights=='TRANSIENT') then
              key_sim=2
           else if (rights=='GMRES_DIAG') then
              key_sim=3
           else if (rights=='GMRES_LU') then
              key_sim=4
           else if (rights=='BCG_DIAG') then
              key_sim=5
           else if (rights=='BCG_LU') then
              key_sim=6
           else if (rights=='CONJ_DIAG') then
              key_sim=7
           else if (rights=='CONJ_CHOLESKY') then
              key_sim=8
           else if (rights=='UMFPACK') then
              key_sim=9
           else
              call prterr('SOLVER must be: TRANSIENT, SOR[_SYM], GMRES_[DIAG|LU], BGC_[DIAG|LU], CONJ_[DIAG|CHOLESKY], UMFPACK')
           endif
              
           key_sim_set=.true.

        else if (lefts.eq.'NSTRIDE') then
           read (rights,*,err=911) nstride
           nstride_set=.true.

        else if (lefts.eq.'NSTRIDENL') then
           read (rights,*,err=911) nstrideNL
           nstrideNL_set=.true.

        else if (lefts.eq.'NSTEP') then
           read (rights,*,err=911) nstep
           nstep_set=.true.

        else if (lefts.eq.'TSTEP') then
           read (rights,*,err=911) tstep
           tstep_set=.true.

        else if (lefts.eq.'EXPORT') then
           if (rights=='NONE') then
              key_export=0
           else if (rights=='PARAVIEW') then
              key_export=1
           else
              call prterr('EXPORT must be "NONE" or "PARAVIEW"')
           endif

           key_export_set=.true.

        else if (lefts.eq.'ZOU_FORMAT') then
           if (rights /= 'ASCII' .and. rights /= 'BINARY') then
              call prterr('ZOU_FORMAT must be "ASCII" or "BINARY"')
           endif

           zou_format=trim(rights)
           zou_format_set=.true.

        else if (lefts.eq.'RST_FORMAT') then
           ! already done

        else if (lefts.eq.'RESIDUAL') then
           if (rights /= 'PIVOTS' .and. rights /= 'UNKNOWNS' .and. rights/='NONE') then
              call prterr('RESIDUAL must be PIVOTS, UNKNOWNS or NONE')
           endif

           residual=trim(rights)
           residual_set=.true.

        else if (lefts.eq.'NODECHECK') then
           if (rights /= 'YES'.and.rights /= 'NO') then
              call prterr('NODECHECK must be YES or NO')
           endif

           nodecheck=trim(rights)
           nodecheck_set=.true.

        else if (lefts.eq.'LABELS') then
           ! already done
        else if (lefts.eq.'TOL') then
           read (rights,*,err=911) tol
           tol_set=.true.

        else if (lefts.eq.'ITMAX') then
           read (rights,*,err=911) itmax
           itmax_set=.true.
           

        else if (lefts.eq.'NEWTON') then
           if (rights /='YES'.and.rights /= 'NO') then
              call prterr('NEWTON must be YES or NO')
           endif

           newton=trim(rights)
           newton_set=.true.

        else if (lefts.eq.'REMOVEFIXED') then
           if (rights /='YES'.and.rights /= 'NO') then
              call prterr('REMOVEFIXED must be YES or NO')
           endif

           removeFixed=trim(rights)
           removeFixed_set=.true.
        else
           call prterr('Unknown control parameter: '//trim(lefts))
        endif

! ----------------------------------------------------------------------
! initialisation. init. needs to be present for each var
! ----------------------------------------------------------------------

     else if (line(1:4)=='INIT') then

        do ivar=1,nvar
           call getline(1,ios,line)
           if (ios.ne.0) call prterr('Unexpected end of file')

           line=adjustl(line)
           if (line(1:3)=="END") call prterr('Not enough statements in INIT block. Each variable must be initialised')

           i=spliteq(line,lefts,vals)
           if (i.eq.0) call prterr('Expecting =')
           
           call findex(lefts,ii)

           if (vals.eq.'INTERP') then
              val=0
              ival=-1
           else
              call parsval(vals,val,ival,list,.false.,.false.,.false.)    ! deriv, norm, istep
              if (ival>2) call prterr('Non-linearity not allowed in initialisation expressions')
              if (ival>1.and.restart=='NONE') call prterr('Restart found in "init" but no restart file nominated')
           endif
           
           inits(ii)=vals
           init(ii) =val
           initi(ii)=ival

           if (ival==1.or.ival==2) then
              allocate(initt(ii)%token(numtok))
              initt(ii)=list
           endif
           
           if (len_trim(vals)>EXPRLEN.and.key_db>0) then
              print *,'WARNING: initt will be truncated'
           endif

        enddo

        endok=.true.   ! if next non-blank line is END that is ok

! ----------------------------------------------------------------------
! region statements, either elements or nodes
! region 1 elements [13 values] C
! region 2 nodes [13 values]
! ----------------------------------------------------------------------

     else if (line(1:6)=='REGION') then

        if (index(line,'ELEMENTS').ne.0) then
           read (line,*,iostat=ios) cdum,irg,keyword,nvals,cdum,cvar

           if (ios /= 0) then
              read (line,*,iostat=ios) cdum,irg,keyword,cvar
              if (ios/=0) call prterr('Incomprehensible REGION..ELEMENTS line')
              nvals=999
           endif

           isetel(irg)=1

           if (.not.(irg.ge.1.and.irg.le.nregel)) then
              call prterr('ELEMENT Region number out of bounds')
           endif

!           write (*,'(4a13)') 'irg','keyword','nvals','cvar'
!           write (*,'(i13,a13,i13,a13)') irg,trim(keyword),nvals,trim(cvar)
        else if (index(line,'NODES').ne.0) then
           read (line,*,iostat=ios) cdum,irg,keyword,nvals

           if (ios/=0) then
              read (line,*,iostat=ios) cdum,irg,keyword
              if (ios/=0) call prterr('Incomprehensible REGION..NODES line')
              nvals=999
           endif

           isetnd(irg)=1

           if (.not.(irg.ge.1.and.irg.le.nregnd)) then
              call prterr('NODE Region number out of bounds')
           endif

!           write (*,'(3a13)') 'irg','keyword','nvals'
!           write (*,'(i13,a13,i13)') irg,trim(keyword),nvals
        else
           call prterr('Region command must contain "elements" or "nodes"')
        endif
                
! ----------------------------------------------------------------------
! elements
! ----------------------------------------------------------------------

        if (keyword.eq.'ELEMENTS') then
           if (cvar.eq.'C') then

              do ivals=1,nvals
                 call getline(1,ios,line)

                 if (ios.ne.0) then
                    call prterr('Unexpected end of file')
                 endif

                 line=adjustl(line)
                 if (line(1:3)=="END") exit

                 i=spliteq(line,lefts,vals)
                 if (i.eq.0) call prterr('Expecting =')

                 call Cindex(lefts,ii,jj,kk,ll,submatrix)

                 ! vals = [1,2,3, 4,5,6, 7,8,9]
                 if (submatrix) then

                    if (ii<1.or.ii>nvar.or.kk<1.or.kk>nvar) then
                       call prterr('C specification indices out of range')
                    endif

                    if (vals(1:1)/='[') then
                       call prterr('Expecting submatrix of form [1,2,3...,9]')
                    endif

                    ist=2
                    ied=2
                    do jj=1,NDIM
                       do ll=1,NDIM
                          found=.false.
                          do i=ist,len_trim(vals)
                             if (vals(i:i)==','.or.vals(i:i)==']') then
                                ied=i-1
                                found=.true.
                                exit
                             endif
                          enddo

                          if (.not.found) call prterr('Unexpected end of line while reading C submatrix')

                          call parsval(vals(ist:ied),val,ival,list,.true.,.false.,.true.)

                          CC(irg,ii,jj,kk,ll)=val
                          iCC(irg,ii,jj,kk,ll)=ival
                          sCC(irg,ii,jj,kk,ll)=vals(ist:ied)

                          if (ival==1.or.ival==2) then
                             allocate(tCC(irg,ii,jj,kk,ll)%token(numtok))
                             tCC(irg,ii,jj,kk,ll)=list
                          endif
                          
                          ist=ied+2
                       enddo
                    enddo

                 else
                                
                    call Ctest(ii,jj,kk,ll,nvar)
                    call parsval(vals,val,ival,list,.true.,.false.,.true.)   ! deriv, norm, istep

                    CC(irg,ii,jj,kk,ll)=val
                    iCC(irg,ii,jj,kk,ll)=ival
                    sCC(irg,ii,jj,kk,ll)=vals
                    
                    if (ival==1.or.ival==2) then
                       allocate(tCC(irg,ii,jj,kk,ll)%token(numtok))
                       tCC(irg,ii,jj,kk,ll)=list
                    endif
                    
                    if (len_trim(vals)>EXPRLEN.and.key_db>0) then
                       print *,'WARNING: sCC will be truncated'
                    endif

                 endif
              enddo
              
           else if (cvar.eq.'A') then
              do ivals=1,nvals
                 call getline(1,ios,line)

                 if (ios.ne.0) then
                    call prterr('Unexpected end of file')
                 endif

                 line=adjustl(line)
                 if (line(1:3)=="END") exit

                 i=spliteq(line,lefts,vals)
                 if (i.eq.0) call prterr('Expecting =')

                 call aindex(lefts,ii,jj)
                 call aatest(ii,jj,nvar)
                 call parsval(vals,val,ival,list,.true.,.false.,.true.)   ! deriv, norm, istep

                 aa(irg,ii,jj)=val
                 iaa(irg,ii,jj)=ival
                 saa(irg,ii,jj)=vals

                 if (ival==1.or.ival==2) then
                    allocate(taa(irg,ii,jj)%token(numtok))
                    taa(irg,ii,jj)=list
                 endif

                 if (len_trim(vals)>EXPRLEN.and.key_db>0) then
                    print *,'WARNING: saa will be truncated'
                 endif

!                 write (*,'(2i13,e13.5)') ii,jj,val
              enddo
              
           else if (cvar.eq.'F') then
              do ivals=1,nvals
                 call getline(1,ios,line)
                 if (ios.ne.0) call prterr('Unexpected end of file')

                 line=adjustl(line)
                 if (line(1:3)=="END") exit

                 i=spliteq(line,lefts,vals)
                 if (i.eq.0) call prterr('Expecting =')

                 call findex(lefts,ii)
                 call fftest(ii,nvar)
                 call parsval(vals,val,ival,list,.true.,.false.,.true.)   ! deriv, norm, istep

                 ff(irg,ii)=val
                 iff(irg,ii)=ival
                 sff(irg,ii)=vals

                 if (ival==1.or.ival==2) then
                    allocate(tff(irg,ii)%token(numtok))
                    tff(irg,ii)=list
                 endif

                 if (len_trim(vals)>EXPRLEN.and.key_db>0) then
                    print *,'WARNING: sff will be truncated'
                 endif
                 
!                 write (*,'(i13,e13.5)') ii,val
              enddo
           else
              call prterr('Unknown material property variable')
           endif

! ----------------------------------------------------------------------
! Nodes
! ----------------------------------------------------------------------

        else if (keyword.eq.'NODES') then

           do ivals=1,nvals
              call getline(1,ios,line)
              if (ios.ne.0) call prterr('Unexpected end of file')

              line=adjustl(line)
              if (line(1:3)=="END") exit

              i=spliteq(line,lefts,vals)
              if (i.eq.0) call prterr('Expecting =')

              nodeinit=.false.
              if (lefts(1:4)=='INIT') then
                 nodeinit=.true.
                 lefts=lefts(5:)
                 lefts=adjustl(lefts)
              endif

              call findex(lefts,ii)
              call fftest(ii,nvar)
              call parsval(vals,val,ival,list,.false.,.false.,.true.)     ! deriv, norm, istep

              regfix(irg,ii)=val

              if (nodeinit) then
                 iregfix(irg,ii)=-(ival+1)       ! -1 = const init
                 if (ival>1) call prterr('Initialisation expression must be const or vary only with XYZ')
              else
                 iregfix(irg,ii)=  ival+1        ! note: 0 = not fixed or init. 1=const etc
              endif

              sregfix(irg,ii)=vals

              if (ival==1.or.ival==2) then
                 allocate (tregfix(irg,ii)%token(numtok))
                 tregfix(irg,ii)=list
              endif
                 
!              write (*,'(i13,e13.5)') ii,val
           enddo
        else
           call prterr('Unknown keyword in Region command')
        endif

! ----------------------------------------------------------------------
! surface 1 2 [3 values] q
! surface 1-2 [3 values] g
! surface 1-xmax [3 values] q
! ----------------------------------------------------------------------

     else if (line(1:7)=='SURFACE') then
        do i=1,len_trim(line)
           if (line(i:i).eq.'-') line(i:i)=' '
        enddo

        print *,trim(line)

        read (line,*,iostat=ios) cdum,sreg1,sreg2,nvals,cdum,cvar

        if (ios /= 0) then
           read (line,*,iostat=ios) cdum,sreg1,sreg2,cvar
           if (ios /=0 ) call prterr('Incomprehensible SURFACE line')
           nvals=999
        endif

        ireg1=lookind(sreg1)
        ireg2=lookind(sreg2)

        ind=lookup(ireg1,ireg2)

        if (ind.eq.0) call prterr('Illegal surface')

        if (cvar.eq.'Q') then
           do ivals=1,nvals
              call getline(1,ios,line)

              if (ios.ne.0) then
                 call prterr('Unexpected end of file')
                 call exit(1)
              endif

              line=adjustl(line)
              if (line(1:3)=="END") exit

              i=spliteq(line,lefts,vals)
              if (i.eq.0) call prterr('Expecting =')

              call aindex(lefts,ii,jj)
              call aatest(ii,jj,nvar)
              call parsval(vals,val,ival,list,.false.,.true.,.true.)    ! deriv,norm,istep

              qq(ind,ii,jj)=val
              iqq(ind,ii,jj)=ival
              sqq(ind,ii,jj)=vals

              if (ival==1.or.ival==2) then
                 allocate(tqq(ind,ii,jj)%token(numtok))
                 tqq(ind,ii,jj)=list
              endif
           enddo
        else if (cvar.eq.'G') then
           do ivals=1,nvals
              call getline(1,ios,line)

              if (ios.ne.0) then
                 call prterr('Unexpected end of file')
                 call exit(1)
              endif

              line=adjustl(line)
              if (line(1:3)=="END") exit

              i=spliteq(line,lefts,vals)
              if (i.eq.0) call prterr('Expecting =')

              call findex(lefts,ii)
              call fftest(ii,nvar)
              call parsval(vals,val,ival,list,.false.,.true.,.true.)   ! deriv, norm, istep

              gg(ind,ii)=val
              igg(ind,ii)=ival
              sgg(ind,ii)=vals

              if (ival==1.or.ival==2) then
                 allocate (tgg(ind,ii)%token(numtok))
                 tgg(ind,ii)=list
              endif
           enddo
        else
           call prterr('Surface variable must be "q" or "g"')
        endif
     else if (line(1:3)=='END') then
        if (endok) then
           endok=.false.
        else
           call prterr('Unexpected END statement outside of block')
        endif
     else
        call prterr('Incomprehensible line. Expecting x=y, REGION, INIT, SURFACE, comment or blank line')
     endif
     
  enddo

  close (1)

! ----------------------------------------------------------------------
! Check all is present and correct
! ----------------------------------------------------------------------

  if (.not.nvar_set) call prterr('NVAR not set')
  if (.not.key_sim_set) call prterr('SOLVER not set')

  if (key_sim.eq.2) then
     if (.not.tstep_set) call prterr('TSTEP not set')
     if (.not.nstep_set) call prterr('NSTEP must be set for transient simulations')
  endif

  if (key_db.gt.1) then
     if (.not.idbmin_set) call prterr('idbmin not set')
     if (.not.idbmax_set) call prterr('idbmax not set')

     if (.not.jdbmin_set) call prterr('jdbmin not set')
     if (.not.jdbmax_set) call prterr('jdbmax not set')

     if (.not.kdbmin_set) call prterr('kdbmin not set')
     if (.not.kdbmax_set) call prterr('kdbmax not set')
  endif

  write (*,*) 'Done.'

! ----------------------------------------------------------------------
! Prepare matrices which specify if matrices are nonlinear
! reg(C,A,F,fix,Q,G)NL
! isQNL, isRNL, isFixNL
! isNL
! ----------------------------------------------------------------------

  allocate (regCNL(nregel),stat=allocerr(1))
  allocate (regANL(nregel),stat=allocerr(2))
  allocate (regFNL(nregel),stat=allocerr(3))
  allocate (regfixNL(nregnd),stat=allocerr(4))

  do i=1,4
     write (temp,*) i
     if (allocerr(i).ne.0) then
        call prterr('NL: Failed to allocate array:'//trim(temp))
     endif
  enddo

  do irg=1,nregel
     regCNL(irg)=any(iCC(irg,:,:,:,:).gt.1)
     regANL(irg)=any(iaa(irg,:,:).gt.1)
     
     regFNL(irg)=any(iff(irg,:).gt.1)
  enddo

  do irg=1,nregnd
     regfixNL(irg)=any(iregfix(irg,:).gt.2)    ! note 2
  enddo

  do itable=1,ntable
     regQNL(itable)=any(iqq(itable,:,:).gt.1)
     regGNL(itable)=any(igg(itable,:).gt.1)
  enddo

  isQNL=any(regCNL).or.any(regANL)
  isRNL=any(regFNL)
  isFixNL=any(regfixNL)

  if (ntable.gt.0) then
     isQNL=isQNL.or.any(regQNL)
     isRNL=isRNL.or.any(regGNL)
  endif

  isNL=isQNL.or.isRNL.or.isFixNL

! ----------------------------------------------------------------------
! sanity checks
! ----------------------------------------------------------------------

  if ((any(regANL).or.any(regQNL).or.any(regQNL).or.any(regGNL).or.any(regfixNL)).and.newton=='YES') then
     call prterr('NEWTON not supported for non-linear A, Q, G matrices or non-linear Dirichlet')
  endif

  if ((any(iCC==2).or.any(iff==2)).and.newton=='YES') then
     call prterr('Inline non-linear expressions not allowed for NEWTON. Use token replacement (DLL) instead')
  endif

  if (key_Q>=2.and.ntable>0) then
     call prterr('Surfaces not supported for key_Q > 2')
  endif

  if (any(aa/=0).and.key_Q>=4) then
     call prterr('A matrix not supported for cuboid elements (key_Q=4)')
  endif

  do ii=1,nvar
     if (initi(ii).eq.-99) call prterr('Variable not initialised')
  enddo

  do i=1,nregel
     if (isetel(i).eq.0) then
        write (*,*) 'Warning: Element region ',i,' not set in ZIN file'
     endif
  enddo

  call range(key_db,'key_db',0,2)
  call drange(omega,'omega',0d0,2d0)

  if (key_db.gt.1) then
     call range(idbmin,'idbmin',0,imax)
     call range(idbmax,'idbmax',0,imax)

     call range(jdbmin,'jdbmin',0,jmax)
     call range(jdbmax,'jdbmax',0,jmax)

     call range(kdbmin,'kdbmin',0,kmax)
     call range(kdbmax,'kdbmax',0,kmax)
  endif

  call range(key_Q,'key_Q',1,5)
  call range(key_sim,'key_sim',0,9)
!  call range(nstride,'nstride',1,nstep)
!  call range(nstrideNL,'nstrideNL',1,nstep)
  call range(key_export,'key_export',0,1)

  deallocate(isetel,isetnd)
  deallocate(list%token)

  return

911 call prterr('Unable to decipher a number/expression in file')

end subroutine readmatrix

! ######################################################################

subroutine wrtmatrix(iunit)

  use common
  use evaluate
  use iofile

  integer iunit,ii,irg,jj,kk,ll,itable,i
  character(20) temp
  double precision val

! ----------------------------------------------------------------------
! write control variables
! ----------------------------------------------------------------------

  write (iunit,'(/a15,a15,a15)') 'Control','Value','Set by user'

  write (iunit,'(a15,i15,l15)') 'nvar',nvar,nvar_set
  write (iunit,'(a15,i15,l15)') 'ndim',ndim,ndim_set

  write (iunit,'(a15,a15,l15)') 'restart',trim(restart),restart_set
  write (iunit,'(a15,a15,l15)') 'newton',trim(newton),newton_set
  write (iunit,'(a15,a15,l15)') 'removefixed',trim(removeFixed),removeFixed_set

  if (key_sim==0) then
     write (iunit,'(a15,a15,l15)') 'solver','SOR_SYM',key_sim_set
  else if (key_sim==1) then
     write (iunit,'(a15,a15,l15)') 'solver','SOR',key_sim_set
  else if (key_sim==2) then
     write (iunit,'(a15,a15,l15)') 'solver','TRANSIENT',key_sim_set
  else if (key_sim==3) then
     write (iunit,'(a15,a15,l15)') 'solver','GMRES_DIAG',key_sim_set
  else if (key_sim==4) then
     write (iunit,'(a15,a15,l15)') 'solver','GMRES_LU',key_sim_set
  else if (key_sim==5) then
     write (iunit,'(a15,a15,l15)') 'solver','BCG_DIAG',key_sim_set
  else if (key_sim==6) then
     write (iunit,'(a15,a15,l15)') 'solver','BCG_LU',key_sim_set
  else if (key_sim==7) then
     write (iunit,'(a15,a15,l15)') 'solver','CONJ_DIAG',key_sim_set
  else if (key_sim==8) then
     write (iunit,'(a15,a15,l15)') 'solver','CONJ_CHOLESKY',key_sim_set
  else if (key_sim==9) then
     write (iunit,'(a15,a15,l15)') 'solver','UMFPACK',key_sim_set
  endif

  write (iunit,'(a15,i15,l15)') 'nstep',nstep,nstep_set
  write (iunit,'(a15,e15.5,l15/)') 'tstep',tstep,tstep_set

  write (iunit,'(a15,a15,l15)') 'zou_format',zou_format,zou_format_set
  write (iunit,'(a15,a15,l15)') 'rst_format',rst_format,rst_format_set
  write (iunit,'(a15,a15,l15)') 'residual',residual,residual_set
  write (iunit,'(a15,a15,l15)') 'nodecheck',nodecheck,nodecheck_set

  write (iunit,'(a15,i15,l15)') 'key_db',key_db,key_db_set
  write (iunit,'(a15,e15.5,l15)') 'omega',omega,omega_set
  write (iunit,'(a15,e15.5,l15)') 'scale',scale,scale_set
  write (iunit,'(a15,i15,l15)') 'key_Q',key_Q,key_Q_set
  write (iunit,'(a15,i15,l15)') 'nstride',nstride,nstride_set
  write (iunit,'(a15,i15,l15)') 'nstrideNL',nstrideNL,nstrideNL_set
  write (iunit,'(a15,i15,l15)') 'itmax',itmax,itmax_set
  write (iunit,'(a15,e15.5,l15)') 'tol',tol,tol_set

  if (key_export==0) then
     write (iunit,'(a15,a15,l15)') 'export','NONE',key_export_set
  else
     write (iunit,'(a15,a15,l15)') 'export','PARAVIEW',key_export_set
  endif

  write (iunit,'(/a15,i15)') '# surfaces',ntable
  write (iunit,'(a15,i15,l15)') 'El. Regions',nregel
  write (iunit,'(a15,i15,l15)') 'Node Regions',nregnd

! ----------------------------------------------------------------------
! write constants
! write matrix details
! ----------------------------------------------------------------------

  write (iunit,'(/a,19x,a13)') 'Symbol','Value'
  do i=1,nparams
     val=real(params(i)%value)
     write (iunit,'(a,a,e13.5)') params(i)%symbol,'=',val
  enddo

! ----------------------------------------------------------------------
! write initial expressions
! ----------------------------------------------------------------------

  write (iunit,*)
  write (iunit,*) 'Initialisation expressions:'
  write (iunit,*) 'Types: 0: number; 1: xyz; 2: non-linear; 3: NL token'
  write (iunit,'(a10,5x,a,10x,a10,a13)') 'Var','Expression','type','value'
  
  do ii=1,nvar
     temp=abbrev(inits(ii),20)
     write (iunit,'(i10,5x,a20,i10,e13.5)') ii,temp,initi(ii),init(ii)
  enddo
  
! ----------------------------------------------------------------------
! Write material properties/pdes data+fixed nodes
! ----------------------------------------------------------------------
  
  do irg=1,nregel
     write (iunit,*)
     write (iunit,*) 'ELEMENT Region ',irg,'=========================================================='
     write (iunit,*) 'Types: 0: number; 1: xyz; 2: non-linear; 3: NL token'
     write (iunit,'(2x,4a5,a13,a8,2x,a)') 'i','j','k','l','value','type','expression'
     do ii=1,nvar
        do kk=1,nvar
           do jj=1,NDIM
              do ll=1,NDIM
                 temp=abbrev(sCC(irg,ii,jj,kk,ll),20)
                 write (iunit,'(a,4i5,e13.5,i8,2x,a)') &
                      'C ',ii,jj,kk,ll,&
                      CC(irg,ii,jj,kk,ll),iCC(irg,ii,jj,kk,ll),temp
              enddo
           enddo
           write (iunit,*)
        enddo
     enddo
     
     write (iunit,'(2x,2a5,a13,a8,2x,a)') 'i','k','value','type','expression'
     do ii=1,nvar
        do kk=1,nvar
           temp=abbrev(saa(irg,ii,kk),20)
           write (iunit,'(a,2i5,e13.5,i8,2x,a)') &
                'a ',ii,kk,aa(irg,ii,kk),iaa(irg,ii,kk),temp
        enddo
     enddo

     write (iunit,'(2x,a5,a13,a8,2x,a)') 'i','value','type','expression'
     do ii=1,nvar
        temp=abbrev(sff(irg,ii),20)
        write (iunit,'(a,i5,e13.5,i8,2x,a)') 'f ',ii,ff(irg,ii),iff(irg,ii),temp
     enddo
  enddo

  do irg=1,nregnd
     write (iunit,*)
     write (iunit,*) 'NODE Region ',irg,'=========================================================='
     write (iunit,*) 'Types: 0: not fixed; 1: number; 2: xyz; 3: non-linear; 4: NL token'
     write (iunit,'(a5,a13,a8,2x,a)') 'i','value','type','expression'
     do ii=1,nvar
        if (iregfix(irg,ii).eq.0) then
           write (iunit,'(i5,a13,i8)') ii,'-',0
        else
           temp=abbrev(sregfix(irg,ii),20)
           write (iunit,'(i5,e13.5,i8,2x,a)') ii,regfix(irg,ii),iregfix(irg,ii),temp
        endif
     enddo
  enddo

  write (iunit,*) '==================END REGION INFO==================='

! ----------------------------------------------------------------------
! Write out surface information
! ----------------------------------------------------------------------

  do itable=1,ntable

     write (iunit,*)
     write (iunit,*) 'SURFACE',itable,'================================='
     write (iunit,*) 'connection:',table(itable,:)

     write (iunit,'(2x,2a5,a13,a8,2x,a)') 'i','k','value','type','expression'
     do ii=1,nvar
        do kk=1,nvar
           temp=abbrev(sqq(itable,ii,kk),20)
           write (iunit,'(a,2i5,e13.5,i8,2x,a)') &
                'q ',ii,kk,qq(itable,ii,kk),iqq(itable,ii,kk),temp
        enddo
     enddo

     write (iunit,'(2x,a5,a13,a8,2x,a)') 'i','value','type','expression'
     do ii=1,nvar
        temp=abbrev(sgg(itable,ii),20)
        write (iunit,'(a,i5,e13.5,i8,2x,a)') 'g ',ii,gg(itable,ii),igg(itable,ii),temp
     enddo
  enddo

! ----------------------------------------------------------------------
! write info about non-linearity
! ----------------------------------------------------------------------

  write (iunit,*)
  write (iunit,*) 'Non linear region settings'
  write (iunit,'(4a5)') 'ireg','CNL','ANL','FNL'
 
  do irg=1,nregel
     write (iunit,'(i5,3l5)') irg,regCNL(irg),regANL(irg),regFNL(irg)
  enddo
  
  if (ntable>0) then
     write (iunit,*) 'Non linear surfaces'
     write (iunit,'(3a5)') 'itable','QNL','GNL'
     do itable=1,ntable
        write (iunit,'(i5,2l5)') itable,regQNL(itable),regGNL(itable)
     enddo
  endif

  write (iunit,*) 'Q contains non-linear terms?',isQNL
  write (iunit,*) 'R contains non-linear terms?',isRNL

end subroutine wrtmatrix

! ######################################################################

subroutine freematrix

  deallocate(CC,iCC,sCC,tCC)
  deallocate(aa,iaa,saa,taa)
  deallocate(ff,iff,sff,tff)
  deallocate(qq,iqq,sqq,tqq)
  deallocate(gg,igg,sgg,tgg)
  deallocate(table,regQNL,regGNL)
  deallocate(regCNL,regANL,regFNL,regfixNL)
  deallocate(iregfix,regfix,sregfix,tregfix)
  deallocate(init,initi,inits,initt)

end subroutine freematrix

! ######################################################################

subroutine nodecheck_snapshot(iunit)

  use common

  integer iunit,i,j,k,ii,ireg1,iel,iregel(8)
  logical different,zerofound

  do i=0,imax
     do j=0,jmax
        do k=0,kmax
           do ii=1,nvar

              ireg1=ireg(i,j,k)

              if (iregfix(ireg1,ii)==1) then   ! fixed literal const
                 if (regfix(ireg1,ii)==0) then
                    call writenode(iunit, i,j,k,ii,1 )          ! deactivated
                 else
                    call writenode(iunit, i,j,k,ii,2 )          ! fixed (non-zero)
                 endif
              else if (iregfix(ireg1,ii)>1) then
                 call writenode(iunit, i,j,k,ii,2 )            ! fixed (non-zero)
              else if (i==0.or.i==imax.or.j==0.or.j==jmax.or.k==0.or.k==kmax) then
                 call writenode(iunit, i,j,k,ii,3 )           ! external Neumann
              else
                 iregel(1)=iregup(i,  j,  k)
                 iregel(2)=iregup(i-1,j,  k)
                 iregel(3)=iregup(i,  j-1,k)
                 iregel(4)=iregup(i,  j,  k-1)
                 iregel(5)=iregup(i-1,j-1,k)
                 iregel(6)=iregup(i-1,j,  k-1)
                 iregel(7)=iregup(i,  j-1,k-1)
                 iregel(8)=iregup(i-1,j-1,k-1)

                 different=.false.
                 do iel=2,8
                    if (.not.equalrow(ii,iregel(iel),iregel(1))) then
                       different=.true.
                       exit
                    endif
                 enddo

                 if (different) then

                    zerofound=.false.
                    do iel=1,8
                       if (zerorow(ii,iregel(iel))) then
                          zerofound=.true.
                          exit
                       endif
                    enddo

                    if (zerofound) then
                       call writenode(iunit, i,j,k,ii,4 )       ! internal Neumann
                    else
                       call writenode(iunit, i,j,k,ii,5 )       ! continuity
                    endif
                 else
                    call writenode(iunit, i,j,k,ii,6 )       ! variable
                 endif
              endif

           enddo
        enddo
     enddo
  enddo
                 
end subroutine nodecheck_snapshot

! ######################################################################

subroutine writenode(iunit,i,j,k,ii,ival)

  use common, only : zou_format

  integer iunit,i,j,k,ii,ival

  if (zou_format=='ASCII') then
     write (iunit,'(5i5)') i,j,k,ii,ival
  else
     write (iunit) i,j,k,ii,ival
  endif

end subroutine writenode

! ######################################################################

function zerorow(ii,ireg1)

  use common, only : nvar,NDIM

  integer ii,ireg1,kk,jj,ll
  logical zerorow

  do kk=1,nvar
     do jj=1,NDIM
        do ll=1,NDIM
           if (.not.zero1(ii,jj,kk,ll,ireg1)) then
              zerorow=.false.
              return
           endif
        enddo
     enddo
  enddo

  zerorow=.true.
end function zerorow

! ######################################################################

function zero1(ii,jj,kk,ll,ireg1)

  integer ii,jj,kk,ll,ireg1
  logical zero1

  if (iCC(ireg1,ii,jj,kk,ll)==0.and.CC(ireg1,ii,jj,kk,ll)==0) then
     zero1=.true.
  else
     zero1=.false.
  endif

end function zero1

! ######################################################################

function equalrow(ii,ireg1,ireg2)

  use common, only : nvar,NDIM

  integer ii,ireg1,ireg2,kk,jj,ll
  logical equalrow

  do kk=1,nvar
     do jj=1,NDIM
        do ll=1,NDIM
           if (.not.equal1(ii,jj,kk,ll,ireg1,ireg2)) then
              equalrow=.false.
              return
           endif
        enddo
     enddo
  enddo

  equalrow=.true.
end function equalrow

! ######################################################################

function equal1(ii,jj,kk,ll,ireg1,ireg2)

  integer ii,jj,kk,ll,ireg1,ireg2
  logical equal1

  equal1=.false.

  if (iCC(ireg1,ii,jj,kk,ll)==0.and.iCC(ireg2,ii,jj,kk,ll)==0) then   ! fixed literal
     if (CC(ireg1,ii,jj,kk,ll) == CC(ireg2,ii,jj,kk,ll)) then
        equal1=.true.
     endif
  else if (iCC(ireg1,ii,jj,kk,ll)==iCC(ireg1,ii,jj,kk,ll)) then       ! xyz, NL or token
     if (sCC(ireg1,ii,jj,kk,ll) == sCC(ireg2,ii,jj,kk,ll)) then
        equal1=.true.
     endif
  endif

end function equal1

! ######################################################################

  subroutine fixQ
    use common
    integer ip

     do ip=1,lenQ
        if (iunk(iQ(ip)) == 0) then
           if (iQ(ip) == jQ(ip)) then
              Qval(ip)=1
           else
              Qval(ip)=0
           endif
        endif
     enddo
   end subroutine fixQ

! ######################################################################

   subroutine fixR
     use common
     integer i
     do i=1,ndof
        if (iunk(i) == 0) then
           RR(i)=vec(i)
        endif
     enddo
   end subroutine fixR

! ######################################################################

   subroutine Qinfo(closed)

! ----------------------------------------------------------------------
! Write information about Q matrix, R matrix and Jacobian
! and other basic info
! ----------------------------------------------------------------------

     use common
     integer, parameter :: nQentries=3000,irow=26

     logical closed
     integer i,ip,ndofx
     double precision Qmax,Rmax

     if (closed) then
        ndofx=ndofred
     else
        ndofx=ndof
     endif

     write (1,*)
     write (1,*) 'Gaussian quadrature order, ng=',ng_xi,ng_eta,ng_mu
     
     write (1,*) 'No. of degrees of freedom: ',ndofx
     write (1,*) 'No. of non-zeros in Q:',lenQ,' out of maximum',27*nvar**2*nnod
     
     write (1,*) 'Allocated ',27*nvar**2*nnod, 'floats in Qval',&
          ' giving ',27*nvar**2*nnod*8*2/1e6,&
          'Mbytes for Qval,iQ,jQ (8 byte double, 4 byte integer)'
     write (1,*) 'Total memory allocation will be somewhat higher than this'
     
     if (key_db == 1) then
        write (1,*) 'Sample of row archive:'

        do i=1,min(ndofx,nQentries)
           write (1,*) 'row:',i,irowst(i),irowed(i),' iunk=',iunk(i)
        enddo

        write (1,*) 'Start of Q matrix:'
        write (1,'(3a12,a14)') 'ip','iQ','jQ','Qval'
        do ip=1,min(lenQ,nQentries)
           write (1,'(3i12,e14.5)') ip,iQ(ip),jQ(ip),Qval(ip)
        enddo

        if (lenQ>nQentries) then
           write (1,*) 'End of Q matrix:'
           write (1,'(3a12,a14)') 'ip','iQ','jQ','Qval'
           do ip=max(lenQ-nQentries,1),lenQ
              write (1,'(3i12,e14.5)') ip,iQ(ip),jQ(ip),Qval(ip)
           enddo
        endif

        if (irow.le.ndofx) then
           write (1,*) 'Sample Row. Row number: ',irow
           write (1,'(3a12,2a14)') 'ip','iQ','jQ','Qval','vec(jQ)'
           do ip=irowst(irow),irowed(irow)
              write (1,'(3i12,2e14.5)') ip,iQ(ip),jQ(ip),Qval(ip),vec(jQ(ip))
           enddo
        endif

        Qmax=0
        do ip=1,lenQ
           Qmax=max(abs(Qval(ip)),Qmax)
        enddo
        write (1,'(a,e13.5)') 'Maximum |Q-value| found (whole matrix): ',Qmax

        write (1,*) 'Sample of R matrix:'
        do i=1,min(ndofx,nQentries)
           write (1,'(i8,e13.5)') i,RR(i)
        enddo

        Rmax=0
        do i=1,ndofx
           Rmax=max(abs(RR(i)),Rmax)
        enddo
        write (1,'(a,e13.5)') 'Maximum |R-value| found (whole vector): ',Rmax

     endif
   end subroutine Qinfo

end module matrices
