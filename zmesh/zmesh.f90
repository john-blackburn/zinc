! instead of rotating pp, just rotate plist. Resort plist when we rotate.
! When visibility is changed, regenerate pp and then get plist.
! But this would mean need to accumulate total rotation when we do need
! to regenerate pp

! We need two legends each of which can be scrolled (scroll icons on screen).

! strange bug where scene is not shown from the front when file is
! loaded. Each time we load a different angle is shown
! This is due to the GetOpenFileName dialog. If you select a file
! by double clicking, Windows sends a WM_MOVE message (with lbuttondown)
! to main window!
! ######################################################################

integer(4) function WinMain(hInstance,hPrevinstance,szCmdLine,iCmdShow)
!DEC$ ATTRIBUTES STDCALL, DECORATE, ALIAS : 'WinMain' :: WinMain

! ----------------------------------------------------------------------
! Visualise an mtf file produced by metamesh.
! Read in file. Plot interface polygons between switched on regions
! User selects element, surface and node regions to activate
! Eg show region 2 elements, region 4 nodes
! Colour code all interface element faces according to the element region
! (but not if surface will be drawn)
! For surfaces, draw element faces and colour according to region of face
! Draw nodes according to region of node. (also need to be sorted)
! Either 3D or 2D mode. In latter, just draw one slice through the system
! slice perp to i,j,k direction. In this case, draw all elements around
! the surface of the "cuboid" slice.
! ----------------------------------------------------------------------

  use user32
  use gdi32
  use kernel32

  use global
  use const
  use win

  implicit none

  integer(HANDLE) hInstance,hPrevInstance
  integer(LPVOID) szCmdLine
  integer(DWORD) iCmdShow

  character(100) :: szAppName= "3dint"C
  character(100) :: szDebug= "debug"C
  type(T_MSG) msg
  type(T_WNDCLASSEX) wndclass

  integer ret

  hInst=hInstance

  wndclass%cbSize = sizeof(wndclass)
  wndclass%style            = 0
!  wndclass%style            = ior(CS_HREDRAW,CS_VREDRAW)
  wndclass%lpfnWndProc      = LOC(WndProc)
  wndclass%cbClsExtra       = 0
  wndclass%cbWndExtra       = 0
  wndclass%hInstance        = hInstance
  wndclass%hIcon            = LoadIcon(GetModuleHandle(NULL), 999)
  wndclass%hCursor          = LoadCursor(NULL, IDC_ARROW)
  wndclass%hbrBackground    = GetStockObject(LTGRAY_BRUSH)
  wndclass%lpszMenuName     = int(100,LPVOID)
  wndclass%lpszClassName    = LOC(szAppName)
  wndclass%hIconSm = LoadImage(GetModuleHandle(NULL), 999, IMAGE_ICON, 16, 16, 0)

  if (RegisterClassEx(wndclass)==0) then
     WinMain=0
     return
  endif
  
  wndclass%lpfnWndProc=LOC(debug)
  wndclass%lpszClassName=LOC(szDebug)
  wndclass%lpszMenuName     = NULL

  if (RegisterClassEx(wndclass)==0) then
     WinMain=0
     return
  endif

  hwndMain = CreateWindowEx (0,szAppName, "ZMesh"C,  &
       WS_OVERLAPPEDWINDOW,          &
       0, 0, &
       stwidth, stheight, &
       NULL, NULL, hInstance, NULL) 
  
  ret=ShowWindow (hwndMain, iCmdShow)
  ret=UpdateWindow (hwndMain)

  hwndDebug = CreateWindowEx (0,szDebug, "Output"C,  &
       WS_OVERLAPPEDWINDOW,          &
       stwidth+1, 0, &
       stconwidth, stheight, &
       NULL, NULL, hInstance, NULL) 
  
  ret=ShowWindow (hwndDebug, iCmdShow)
  ret=UpdateWindow (hwndDebug)

  message='Welcome to ZMesh 3.6 (2016)'
  ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

  do while (GetMessage (msg, NULL, 0, 0)>0)
     if (hDlgHelp.ne.0) then
        if (IsDialogMessage(hDlgHelp,msg).eq.FALSE) then
           ret=TranslateMessage(msg) 
           ret=DispatchMessage(msg)
        endif
     else
        ret=TranslateMessage(msg) 
        ret=DispatchMessage(msg)
     endif
  enddo

  WinMain=msg%wParam

end function WinMain

! ######################################################################

subroutine swap(a,i,j)
  
  use global

  real a(0:*),temp,tempv(4,3)
  integer i,j,i1,j1,itemp
  
  temp=a(i)
  a(i)=a(j)
  a(j)=temp
  
  i1=i+1
  j1=j+1
  
  tempv=plist(i1,:,:)
  plist(i1,:,:)=plist(j1,:,:)
  plist(j1,:,:)=tempv

  itemp=preg(i1)
  preg(i1)=preg(j1)
  preg(j1)=itemp
end subroutine swap

! ######################################################################

subroutine snaperror(s)
  use kernel32
  use user32
  use global, only : message,hwndDebug
  use const

  integer ret
  character(*) s

  message=s
  ret=SendMessage(hwndDebug,WM_WDEBUG,1,0)
end subroutine snaperror

subroutine createmessage(s)

  use kernel32
  use user32
  use global, only : message,hwndDebug
  use const

  integer ret
  character(*) s

  message=s
  ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

end subroutine createmessage

subroutine createerror(s)

  use kernel32
  use user32
  use global, only : message,hwndDebug
  use const

  integer ret
  character(*) s

  message=s
  ret=SendMessage(hwndDebug,WM_WDEBUG,1,0)

end subroutine createerror
