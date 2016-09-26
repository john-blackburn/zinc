module global

  use ifwinty
  save

!  integer, parameter :: HANDLE=4, MAX_PATH=100  for gfortran command prompt

  logical :: fastdraw=.false.
  
  integer imax,jmax,kmax,npoly,nreg,nregnode,nvar,ivar,nview,iview
  real, allocatable :: rnode(:,:,:,:),pp(:,:,:,:),plist(:,:,:),pxcent(:),u(:,:,:,:)
  integer, allocatable :: iregnd(:,:,:),iregup(:,:,:),iprtnd(:,:,:),iprtup(:,:,:)
  integer, allocatable :: showreg(:),nodeshowreg(:)

  integer, allocatable :: preg(:)

  integer(HANDLE) hwndMain, hwndDebug,hInst
  integer(HANDLE) :: hDlgHelp=0
  character(100) message

  integer(HANDLE), allocatable :: cols(:),nodecols(:)
  integer(HANDLE) :: btnon,btnoff,legend(0:256)

  integer(HANDLE) :: excols(100),extrabrush(100)
  integer nexcols,nextra
  logical :: extra=.false., showextra=.false., showgridlines=.true.
  real, allocatable :: rextra(:,:,:),rextraold(:,:,:)
  integer, allocatable :: extraind(:),extratype(:),extraview(:)

  integer icutmax,icutmin,jcutmax,jcutmin,kcutmax,kcutmin
  integer xscale,yscale,zscale

  integer slice   ! 0, full, 1, 2, 3 =x,y,z

  integer :: btnstate(7)=[0,0,0, 0,0,0,1]
  integer btntxtpos(7,2)

  real bb(3,4,2,3),ff(3,2,3)

  character(100) :: ffstr(3,2)

  real angx,angy,angz,scale,redge(3,2)
  integer xc,yc,width,height,cxChar,cyChar

  character(MAX_PATH) MinName,StemNAME

  logical :: loaded=.false.,zou=.false.,resid=.false.

  real, allocatable :: umax(:),umin(:),uspan(:),umaxUser(:),uminUser(:),uspanUser(:)
  character(20), allocatable :: label(:)

  integer ScrollEl,ScrollNd

contains

! ######################################################################

  subroutine free
    deallocate (rnode,pp,plist,pxcent,iregup,iregnd,iprtup,iprtnd,preg)
  end subroutine free
    
! ######################################################################

  subroutine alloc

    integer mxpoly,mxnode,allocerr(10),i

    allocate (rnode(0:imax,0:jmax,0:kmax,3),pp(0:imax,0:jmax,0:kmax,3),stat=allocerr(1))
    allocate (iregup(0:imax,0:jmax,0:kmax),stat=allocerr(2))
    allocate (iregnd(0:imax,0:jmax,0:kmax),stat=allocerr(3))
    allocate (iprtup(0:imax,0:jmax,0:kmax),stat=allocerr(4))
    allocate (iprtnd(0:imax,0:jmax,0:kmax),stat=allocerr(5))
    
    mxpoly=imax*jmax*(kmax+1)+jmax*kmax*(imax+1)+imax*kmax*(jmax+1)
    mxnode=(imax+1)*(jmax+1)*(kmax+1)
    
    allocate (plist(mxpoly+mxnode+nextra,4,3),pxcent(mxpoly+mxnode+nextra),preg(mxpoly+mxnode+nextra),stat=allocerr(6))
    
    do i=1,6
       if (allocerr(i).ne.0) call snaperror('(global) Failed to allocate rnode etc')
    enddo

  end subroutine alloc

end module global
