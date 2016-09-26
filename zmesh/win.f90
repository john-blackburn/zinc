! WndProc: main window procedure
! debug: debug window
! dlgproc: cropping dialog
! dlgcol: colour selection (superseded by comdlg ChooseColor)
! getpoly
! hiword,loword,rgb: replace C macros
! getDlgItemInt. success should be pointer to an integer(BOOL). or NULL
! movetoEx is ok because returns a struct
! in general need to be careful with non-pure Windows functions

! also in case of fastdraw getpoly needs to be called when toggling extra

module win

  implicit none

contains

! ######################################################################

  integer(4) function WndProc(hwnd,iMsg,wParam,lParam)
    !DEC$ ATTRIBUTES STDCALL :: WndProc

    use user32
    use gdi32
    use kernel32
    use comdlg32
    use gdi32, Polygon => MSFWIN$Polygon

    use const
    use global
    use geo
    use createmesh

    type(T_POINT) pline(4)

    integer(HANDLE) hwnd,wParam,lParam,hdc,hMenu
    integer(UINT) iMsg

    type (T_TEXTMETRIC) tm
    type (T_PAINTSTRUCT) ps

    integer, save :: xmouse, ymouse

    integer i,j,k,ret,l,r,g,b,ind,nshow,iS
    integer xarrow1,xarrow2,yarrow1,yarrow2

    real x,y,z
    real(8) u1
    integer ic,jc,kc,mreg,mregup,ioss(3),allocerr,xview
    integer colour,xint,yint,xmouse2,ymouse2,ios,i1,j1,k1,ii,ii1
    character(20) cdum,zou_format

    type(T_OPENFILENAME) :: ofn
    character(MAX_PATH) filter,PathName,extraName

    integer(HANDLE), save :: hdcMem=NULL, hbitMem=NULL

    type(T_CHOOSECOLOR) cc
    integer(ULONG)   crCustColors(16)
    integer(HANDLE) hdcEMF,hemf
    character(MAX_PATH) metafile
    character(1000) line,left,right
    type(T_RECT), save :: rect
    integer(HANDLE) hRgnClip
    integer type,val,ip
    logical found

    WndProc=0

    if (iMsg.eq.WM_CREATE) then

! ----------------------------------------------------------------------
! WM_CREATE. Get text metrics. Set button text positions
! ----------------------------------------------------------------------

       hdc=GetDC(hwnd)

       ret=SelectObject(hdc,GetStockObject(SYSTEM_FIXED_FONT))

       ret=GetTextMetrics(hdc,tm)
       cxChar=tm%tmAveCharWidth
       cyChar=tm%tmHeight

       ret=ReleaseDC(hwnd,hdc)

       do i=1,7
          btntxtpos(i,1)=btnwidth/2-len_trim(btntext(i))*cxChar/2
          btntxtpos(i,2)=btnheight/2-cyChar/2
       enddo

       slice=0

       btnoff=CreateSolidBrush(rgb(0,255,0))
       btnon= CreateSolidBrush(rgb(255,255,0))

! ----------------------------------------------------------------------
! Set angles, scale, speed
! ----------------------------------------------------------------------

       angx=10*pi/180
       angy=10*pi/180
       angz=10*pi/180

       scale=100

!       do i=0,255
!          legend(i)=CreateSolidBrush(rgb(i,i,i))
!       enddo

       ind=0
       do i=0,63        ! Red to yellow  0..63
          r=255; g=i*4; b=0
          legend(ind)=CreateSolidBrush(rgb(r,g,b))
          ind=ind+1
       enddo

       do i=0,63        ! Yellow to green   64..127
          r=255-i*4; g=255; b=0
          legend(ind)=CreateSolidBrush(rgb(r,g,b))
          ind=ind+1
       enddo

       do i=0,63        ! Green to cyan    128..191
          r=0; g=255; b=i*4
          legend(ind)=CreateSolidBrush(rgb(r,g,b))
          ind=ind+1
       enddo

       do i=0,63        ! cyan to blue   192..255
          r=0; g=255-i*4; b=255
          legend(ind)=CreateSolidBrush(rgb(r,g,b))
          ind=ind+1
       enddo

       legend(0:255)=legend(255:0:-1)
       legend(256)=GetStockObject(WHITE_BRUSH)

       xscale=1
       yscale=1
       zscale=1

       return

    else if (iMsg.eq.WM_SIZE) then

! ----------------------------------------------------------------------
! WM_SIZE
! ----------------------------------------------------------------------

       width= loword(lParam)
       height=hiword(lParam)

       rect%top=0
       rect%left=0
       rect%right=loword(lParam)
       rect%bottom=hiword(lParam)

       hdc=GetDC(hwnd)

       if (hdcMem.ne.NULL) then
          ret=DeleteDC(hdcMem)
          ret=DeleteObject(hbitMem)
       endif

       hdcMem=CreateCompatibleDC(hdc)
       hbitMem=CreateCompatibleBitmap(hdc,width,height)
       ret=SelectObject(hdcMem,hbitMem)

       ret=ReleaseDC(hwnd,hdc)

       if (.not.loaded) return

       call rescale(width,height)

       ret=InvalidateRect(hwnd,NULL,FALSE)

       return

    else if (iMsg.eq.WM_RBUTTONDOWN) then

! ----------------------------------------------------------------------
! WM_RBUTTONDOWN. Region buttons colour selection
! ----------------------------------------------------------------------

       if (.not.loaded) return

       xmouse=loword(lParam)
       ymouse=hiword(lParam)

       nshow=min(nreg-ScrollEl+1,maxbutton)

       do i=1,nshow
          xint=width-keyxoff-keywidth
          yint=keyyoff+(i-1)*keypitch       

          if (xmouse.gt.xint.and.xmouse.lt.xint+keywidth.and. &
              ymouse.gt.yint.and.ymouse.lt.yint+keyheight) then

!             colour=DialogBoxParam(hInst,int(300,lpvoid),hwnd,loc(DlgCol),0)

             iS=i+ScrollEl-1

             cc%lStructSize    = sizeof (cc) ;
             cc%hwndOwner      = hwnd;
             cc%hInstance      = NULL ;
             cc%rgbResult      = 255
             cc%lpCustColors   = loc(crCustColors) ;
             cc%Flags          = ior(CC_RGBINIT,CC_FULLOPEN)
             cc%lCustData      = 0 ;
             cc%lpfnHook       = NULL ;
             cc%lpTemplateName = NULL ;

             ret=ChooseColor(cc)
             colour=cc%rgbResult

!             if (colour.ne.-1) then
             if (ret.ne.0) then
                ret=DeleteObject(cols(iS))
                cols(iS)=CreateSolidBrush(colour)
                ret=InvalidateRect(hwnd,NULL,FALSE)
             endif
             exit
          endif
       enddo

! ----------------------------------------------------------------------
! node regions
! ----------------------------------------------------------------------

       nshow=min(nregnode-ScrollNd+1,maxbutton)

       do i=1,nshow
          xint=width-2*keyxoff-2*keywidth
          yint=keyyoff+(i-1)*keypitch       

          if (xmouse.gt.xint.and.xmouse.lt.xint+keywidth.and. &
              ymouse.gt.yint.and.ymouse.lt.yint+keyheight) then

!             colour=DialogBoxParam(hInst,int(300,lpvoid),hwnd,loc(DlgCol),0)

             iS=i+ScrollNd-1

             cc%lStructSize    = sizeof (cc) ;
             cc%hwndOwner      = hwnd;
             cc%hInstance      = NULL ;
             cc%rgbResult      = 255
             cc%lpCustColors   = loc(crCustColors) ;
             cc%Flags          = ior(CC_RGBINIT,CC_FULLOPEN)
             cc%lCustData      = 0 ;
             cc%lpfnHook       = NULL ;
             cc%lpTemplateName = NULL ;

             ret=ChooseColor(cc)
             colour=cc%rgbResult

!             if (colour.ne.-1) then
             if (ret.ne.0) then
                ret=DeleteObject(nodecols(iS))
                nodecols(iS)=CreateSolidBrush(colour)
                ret=InvalidateRect(hwnd,NULL,FALSE)
             endif
             exit
          endif
       enddo

       return

    else if (iMsg.eq.WM_LBUTTONDOWN) then

! ----------------------------------------------------------------------
! WM_LBUTTONDOWN. Check region buttons
! account for ScrollEl, ScrollNd which points to the first button to be shown
! maxbutton is max number of colour buttons shown at once
! ----------------------------------------------------------------------

       if (.not.loaded) return

       xmouse=loword(lParam)
       ymouse=hiword(lParam)

       nshow=min(nreg-ScrollEl+1,maxbutton)

       do i=1,nshow
          xint=width-keyxoff-keywidth
          yint=keyyoff+(i-1)*keypitch       

          if (xmouse.gt.xint.and.xmouse.lt.xint+keywidth.and. &
              ymouse.gt.yint.and.ymouse.lt.yint+keyheight) then

             iS=i+ScrollEl-1

             if (zou) then
                if (showreg(iS).eq.0) then
                   showreg(iS)=1
                else if (showreg(iS).eq.1) then
                   showreg(iS)=2
                else
                   showreg(iS)=0
                endif
             else
                showreg(iS)=1-showreg(iS)
             endif

             if (fastdraw) then
                if (zou) then
                   call getpolyzou
                else
                   call getpoly
                endif
             endif

             ret=InvalidateRect(hwnd,NULL,FALSE)
             exit
          endif
       enddo

       nshow=min(nregnode-ScrollNd+1,maxbutton)

       do i=1,nshow
          xint=width-2*keyxoff-2*keywidth
          yint=keyyoff+(i-1)*keypitch       

          if (xmouse.gt.xint.and.xmouse.lt.xint+keywidth.and. &
              ymouse.gt.yint.and.ymouse.lt.yint+keyheight) then

             iS=i+ScrollNd-1

             if (zou) then
                if (nodeshowreg(iS).eq.0) then
                   nodeshowreg(iS)=1
                else if (nodeshowreg(iS).eq.1) then
                   nodeshowreg(iS)=2
                else
                   nodeshowreg(iS)=0
                endif
             else
                nodeshowreg(iS)=1-nodeshowreg(iS)
             endif

             ret=InvalidateRect(hwnd,NULL,FALSE)
             exit
          endif
       enddo

! ----------------------------------------------------------------------
! Scroll buttons
! ----------------------------------------------------------------------

       xarrow1=width-keyxoff-keywidth
       xarrow2=xarrow1+keywidth

       yarrow1=0
       yarrow2=keyyoff/2

       if (xmouse.gt.xarrow1.and.xmouse.lt.xarrow2.and.&
           ymouse.gt.yarrow1.and.ymouse.lt.yarrow2) then

          if (ScrollEl.gt.1) ScrollEl=ScrollEl-1
          ret=InvalidateRect(hwnd,NULL,FALSE)
       endif

       yarrow1=keyyoff*1.5+(maxbutton-1)*keypitch+keyheight
       yarrow2=yarrow1+keyyoff/2

       if (xmouse.gt.xarrow1.and.xmouse.lt.xarrow2.and.&
           ymouse.gt.yarrow1.and.ymouse.lt.yarrow2) then

          if (ScrollEl.lt.nreg-maxbutton+1) ScrollEl=ScrollEl+1
          ret=InvalidateRect(hwnd,NULL,FALSE)
       endif

       xarrow1=width-2*keyxoff-2*keywidth
       xarrow2=xarrow1+keywidth

       yarrow1=0
       yarrow2=keyyoff/2

       if (xmouse.gt.xarrow1.and.xmouse.lt.xarrow2.and.&
           ymouse.gt.yarrow1.and.ymouse.lt.yarrow2) then

          if (scrollNd.gt.1) ScrollNd=ScrollNd-1
          ret=InvalidateRect(hwnd,NULL,FALSE)
       endif

       yarrow1=keyyoff*1.5+(maxbutton-1)*keypitch+keyheight
       yarrow2=yarrow1+keyyoff/2

       if (xmouse.gt.xarrow1.and.xmouse.lt.xarrow2.and.&
           ymouse.gt.yarrow1.and.ymouse.lt.yarrow2) then

          if (ScrollNd.lt.nregnode-maxbutton+1) ScrollNd=ScrollNd+1
          ret=InvalidateRect(hwnd,NULL,FALSE)
       endif

! ----------------------------------------------------------------------
! Check bottom row of buttons
! ----------------------------------------------------------------------

       do i=1,7
          xint=btnxoff+(i-1)*btnpitch
          yint=height-btnyoff-btnheight

          if (xmouse.gt.xint.and.xmouse.lt.xint+btnwidth.and. &
              ymouse.gt.yint.and.ymouse.lt.yint+btnwidth) then

             if (i.eq.1) then
                ! pp=rnode
                pp(:,:,:,1)=rnode(:,:,:,1)*xscale
                pp(:,:,:,2)=rnode(:,:,:,2)*yscale
                pp(:,:,:,3)=rnode(:,:,:,3)*zscale
                   
                call setbb
                ! rextra=rextraold
                rextra(:,:,1)=rextraold(:,:,1)*xscale
                rextra(:,:,2)=rextraold(:,:,2)*yscale
                rextra(:,:,3)=rextraold(:,:,3)*zscale

                if (fastdraw) then
                   if (zou) then
                      call getpolyzou
                   else
                      call getpoly
                   endif
                endif

             else if (i.eq.2) then
                ! pp=rnode
                pp(:,:,:,1)=rnode(:,:,:,1)*xscale
                pp(:,:,:,2)=rnode(:,:,:,2)*yscale
                pp(:,:,:,3)=rnode(:,:,:,3)*zscale

                call setbb
                ! rextra=rextraold
                rextra(:,:,1)=rextraold(:,:,1)*xscale
                rextra(:,:,2)=rextraold(:,:,2)*yscale
                rextra(:,:,3)=rextraold(:,:,3)*zscale

                call rotz(pi/2)

                if (fastdraw) then
                   if (zou) then
                      call getpolyzou
                   else
                      call getpoly
                   endif
                endif
             else if (i.eq.3) then
                ! pp=rnode
                pp(:,:,:,1)=rnode(:,:,:,1)*xscale
                pp(:,:,:,2)=rnode(:,:,:,2)*yscale
                pp(:,:,:,3)=rnode(:,:,:,3)*zscale

                call setbb
!                rextra=rextraold
                rextra(:,:,1)=rextraold(:,:,1)*xscale
                rextra(:,:,2)=rextraold(:,:,2)*yscale
                rextra(:,:,3)=rextraold(:,:,3)*zscale

                call roty(pi/2)

                if (fastdraw) then
                   if (zou) then
                      call getpolyzou
                   else
                      call getpoly
                   endif
                endif

             else if (i.ge.4.and.i.le.7) then
                if (btnstate(i).eq.1) then
                   return
                else
                   btnstate(4:7)=0
                   btnstate(i)=1
                endif

                if (i.eq.4) then
                   icutmin=0
                   icutmax=1
                   
                   jcutmin=0
                   jcutmax=jmax
                   
                   kcutmin=0
                   kcutmax=kmax
                   
                   slice=1
                else if (i.eq.5) then
                   icutmin=0
                   icutmax=imax
                   
                   jcutmin=0
                   jcutmax=1
                   
                   kcutmin=0
                   kcutmax=kmax
                   
                   slice=2
                else if (i.eq.6) then
                   icutmin=0
                   icutmax=imax
                   
                   jcutmin=0
                   jcutmax=jmax
                   
                   kcutmin=0
                   kcutmax=1
                   
                   slice=3

                else if (i.eq.7) then
                   icutmin=0
                   icutmax=imax
                   
                   jcutmin=0
                   jcutmax=jmax
                   
                   kcutmin=0
                   kcutmax=kmax

                   slice=0
                endif

                if (fastdraw) then
                   if (zou) then
                      call getpolyzou
                   else
                      call getpoly
                   endif
                endif

             endif
             ret=InvalidateRect(hwnd,NULL,FALSE)
             exit
          endif
       enddo

! ----------------------------------------------------------------------
! Variable change buttons (ZOU File)
! Change legend limits
! ----------------------------------------------------------------------

       if (zou) then
          if (xmouse.gt.xvar.and.xmouse.lt.xvar+varsize.and. &
              ymouse.gt.yvar.and.ymouse.lt.yvar+varsize) then
             if (ivar.gt.1) then
                ivar=ivar-1
                ret=InvalidateRect(hwnd,NULL,FALSE)
                return
             endif
          endif
          
          if (xmouse.gt.xvar+varsize+varsep.and.xmouse.lt.xvar+varsize*2+varsep.and. &
              ymouse.gt.yvar.and.ymouse.lt.yvar+varsize) then
             if (ivar.lt.nvar) then
                ivar=ivar+1
                ret=InvalidateRect(hwnd,NULL,FALSE)
                return
             endif
          endif

          if (xmouse.gt.legendxoff+legendwidth+5.and.xmouse.lt.legendxoff+legendwidth+68.and. &
              ymouse.gt.legendyoff.and.ymouse.lt.legendyoff+20) then
             ret=DialogBoxParam(hInst,int(700,lpvoid),hwnd,loc(DlgSetMax),0)
             if (ret.eq.702.or.ret.eq.704) then       ! OK or Max
                ret=InvalidateRect(hwnd,NULL,FALSE)
                return
             endif
          endif

          if (xmouse.gt.legendxoff+legendwidth+5.and.xmouse.lt.legendxoff+legendwidth+68.and. &
              ymouse.gt.height-legendyoff-35.and.ymouse.lt.height-legendyoff-15) then
             ret=DialogBoxParam(hInst,int(800,lpvoid),hwnd,loc(DlgSetMin),0)
             if (ret.eq.802.or.ret.eq.804) then       ! OK or Min
                ret=InvalidateRect(hwnd,NULL,FALSE)
                return
             endif
          endif

          if (iand(wParam,MK_CONTROL).ne.0) then
             do ip=npoly,1,-1

                if (preg(ip)>0) then
                   do i=1,4
                      pline(i)%x= plist(ip,i,2)*scale+xc
                      pline(i)%y=-plist(ip,i,3)*scale+yc
                   enddo
  
                   if (pointInPolygon(pline,xmouse,ymouse)) then
                      hdc=GetDC(hwnd)
                      ret=SelectObject(hdc,GetStockObject(GRAY_BRUSH))
                      ret=Polygon(hdc,pline(1),4)
                      ret=ReleaseDC(hwnd,hdc)

                      if (preg(ip)>nreg.and.preg(ip).ne.256+nreg+1) then
                         write (cdum,'(a,g13.5)') 'value=',(preg(ip)-nreg-1)*uspanUser(ivar)/255.0+uminUser(ivar)
                         call winmess(cdum)
                      endif

                      exit
                   endif

                endif
             enddo
          endif

       endif

! ----------------------------------------------------------------------
! View change buttons in case of extra
! ----------------------------------------------------------------------

       if (extra.and.showextra) then
          if (xmouse.gt.xvar.and.xmouse.lt.xvar+varsize.and. &
              ymouse.gt.yvar.and.ymouse.lt.yvar+varsize) then
             if (iview.gt.1) then
                iview=iview-1
                ret=InvalidateRect(hwnd,NULL,FALSE)
             endif
          endif
          
          if (xmouse.gt.xvar+varsize+varsep.and.xmouse.lt.xvar+varsize*2+varsep.and. &
              ymouse.gt.yvar.and.ymouse.lt.yvar+varsize) then
             if (iview.lt.nview) then
                iview=iview+1
                ret=InvalidateRect(hwnd,NULL,FALSE)
             endif
          endif
       endif

       return

! ----------------------------------------------------------------------
! WM_COMMAND. First, quit and bounding box dialog
! ----------------------------------------------------------------------

    else if (iMsg.eq.WM_COMMAND) then

       if (loword(wParam).eq.112) then
          hMenu=GetMenu(hwnd)

          if (showextra) then
             ret=CheckMenuItem(hMenu,112,MF_UNCHECKED)
             showextra=.false.
             ret=InvalidateRect(hwnd,NULL,FALSE)
          else
             ret=CheckMenuItem(hMenu,112,MF_CHECKED)
             showextra=.true.
             ret=InvalidateRect(hwnd,NULL,FALSE)
          endif

       else if (loword(wParam).eq.113) then
          hMenu=GetMenu(hwnd)

          if (showgridlines) then
             ret=CheckMenuItem(hMenu,113,MF_UNCHECKED)
             showgridlines=.false.
             ret=InvalidateRect(hwnd,NULL,FALSE)
          else
             ret=CheckMenuItem(hMenu,113,MF_CHECKED)
             showgridlines=.true.
             ret=InvalidateRect(hwnd,NULL,FALSE)
          endif

       else if (loword(wParam).eq.101) then
          ret=SendMessage(hwnd, WM_CLOSE, 0, 0)

       else if (loword(wParam).eq.102) then           ! Bounding box
          ret=DialogBoxParam(hInst,int(200,lpvoid),hwnd,loc(DlgProc),0)
          ret=InvalidateRect(hwnd,NULL,FALSE)

       else if (loword(wParam).eq.114) then           ! distort
          if (DialogBoxParam(hInst,int(600,lpvoid),hwnd,loc(DlgDistort),0).eq.604) then  ! OK
             pp(:,:,:,1)=rnode(:,:,:,1)*xscale
             pp(:,:,:,2)=rnode(:,:,:,2)*yscale
             pp(:,:,:,3)=rnode(:,:,:,3)*zscale

             call setbb
             ! rextra=rextraold
             rextra(:,:,1)=rextraold(:,:,1)*xscale
             rextra(:,:,2)=rextraold(:,:,2)*yscale
             rextra(:,:,3)=rextraold(:,:,3)*zscale

             ret=InvalidateRect(hwnd,NULL,FALSE)
          endif
             
! ----------------------------------------------------------------------
! Help and Help about
! ----------------------------------------------------------------------

       else if (loword(wParam).eq.110) then
!          ret=DialogBoxParam(hInst,int(500,lpvoid),hwnd,loc(DlgHelp),0)
          hDlgHelp=CreateDialogParam(hInst,int(500,lpvoid),hwnd,loc(DlgHelp),0)
          ret=ShowWindow(hDlgHelp,SW_SHOW)
       else if (loword(wParam).eq.111) then
             ret=MessageBox(hwnd,'Zinc 3.6 (Zmesh). NPL, 2016'C,'About Zinc'C,MB_OK)

! ----------------------------------------------------------------------
! Export View
! ----------------------------------------------------------------------

       else if (loword(wParam).eq.103) then

          filter="Enhanced MetaFile (*.emf)\0*.emf\0"C
          metafile=""C

          ofn%lStructSize = sizeof(ofn)
          ofn%hwndOwner = GetForegroundWindow()
          ofn%hInstance = NULL

          ofn%lpstrFilter = LOC(filter)

          ofn%lpstrCustomFilter = NULL
          ofn%nMaxCustFilter = 0
          ofn%nFilterIndex=1

          ofn%lpstrFile = LOC(metafile)
          ofn%nMaxFile = sizeof(metafile)

          ofn%lpstrFileTitle = NULL
          ofn%nMaxFileTitle = 0
          ofn%lpstrInitialDir=NULL
          ofn%lpstrTitle=loc("Select file"C)

          ofn%Flags = NULL

          ofn%nFileOffset=0
          ofn%nFileExtension=0

          ofn%lpstrDefExt = LOC("txt"C)

          ofn%lCustData=0
          ofn%lpfnHook=NULL
          ofn%lpTemplateName=NULL

          if (GetSaveFileName(ofn).eq.TRUE) then
             hdcEMF=CreateEnhMetafile(NULL,metafile,NULL,'Zmesh\0View\0'C)

             hRgnClip=CreateRectRgn(0,0,width,height)
             ret=SelectClipRgn(hdcEMF,hRgnClip)

             call drawit(hdcEMF)
             hemf=CloseEnhMetaFile(hdcEMF)
             ret=DeleteEnhMetaFile(hemf)

             call winmess('Writing EMF file:')
             call winmess(metafile)
          endif
          
! ----------------------------------------------------------------------
! Copy view
! ----------------------------------------------------------------------

       else if (loword(wParam).eq.108) then

          hdcEMF=CreateEnhMetafile(NULL,NULL,NULL,NULL)
          
          hRgnClip=CreateRectRgn(0,0,width,height)
          ret=SelectClipRgn(hdcEMF,hRgnClip)

          call drawit(hdcEMF)
          hemf=CloseEnhMetaFile(hdcEMF)
          
          ret=OpenClipboard(hwnd)
          ret=EmptyClipboard()
          ret=SetClipboardData(CF_ENHMETAFILE,hemf)
          ret=CloseClipboard()

          call winmess('View copied to clipboard')

! do NOT delete the metafile. hemf now belongs to the OS

! ----------------------------------------------------------------------
! Process Mesh: Read in mesh and process it. We already know filename.
! ----------------------------------------------------------------------

       else if (loword(wParam).eq.107) then

          if (processMIN().ne.0) return   ! allocates rnode etc

! ----------------------------------------------------------------------
! Figure out how many regions there are. Set up the key
! ----------------------------------------------------------------------

          nreg=0; nregnode=0
          do i=0,imax
             do j=0,jmax
                do k=0,kmax
                   mregup=iregup(i,j,k)
                   mreg  =iregnd(i,j,k)
                   nreg=max(nreg,mregup)
                   nregnode=max(nregnode,mreg)
                enddo
             enddo
          enddo

!          write (cdum,'(i5,i5)') nreg,nregnode
!          call winmess('Found:'//cdum)

          if (allocated(cols)) call freekey
          call allockey

          zou=.false.
          extra=.false.

          call initview(hwnd)

! ----------------------------------------------------------------------
! Open MIN file. This only displays the MIN file and records the filename
! When we hit "process" this file is loaded again.
! ----------------------------------------------------------------------

       else if (loword(wParam).eq.105) then

          filter="Mesh Input File (*.min)\0*.min\0"C
          MinName=""C

          ofn%lStructSize = sizeof(ofn)
          ofn%hwndOwner = GetForegroundWindow()
          ofn%hInstance = NULL

          ofn%lpstrFilter = LOC(filter)

          ofn%lpstrCustomFilter = NULL
          ofn%nMaxCustFilter = 0
          ofn%nFilterIndex=1

          ofn%lpstrFile = LOC(MinName)
          ofn%nMaxFile = sizeof(MinName)

          ofn%lpstrFileTitle = NULL
          ofn%nMaxFileTitle = 0
          ofn%lpstrInitialDir=NULL
          ofn%lpstrTitle=loc("Select file"C)

          ofn%Flags = OFN_PATHMUSTEXIST

          ofn%nFileOffset=0
          ofn%nFileExtension=0

          ofn%lpstrDefExt = LOC("txt"C)

          ofn%lCustData=0
          ofn%lpfnHook=NULL
          ofn%lpTemplateName=NULL

!          message=filter
!          ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

! ----------------------------------------------------------------------
! Set StemName from MinName, light up View and Process options
! Write to Output window. Open MIN file in view dialog
! ----------------------------------------------------------------------

          if (GetOpenFileName(ofn).eq.TRUE) then
             l=index(MinName,'.min')
             StemName=MinName(:l-1)

             message='Opened file:'
             ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

             message=MinName
             ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

             hMenu=GetMenu(hwnd)
             ret=EnableMenuItem(hMenu,106,MF_ENABLED)

             hMenu=GetMenu(hwnd)
             ret=EnableMenuItem(hMenu,107,MF_ENABLED)
             !light up view, process

             ret=DialogBoxParam(hInst,int(400,lpvoid),hwnd,loc(DlgView),0)

!             loaded=.false.
          endif

! ----------------------------------------------------------------------
! View MIN file
! ----------------------------------------------------------------------

       else if (loword(wParam).eq.106) then
          
          ret=DialogBoxParam(hInst,int(400,lpvoid),hwnd,loc(DlgView),0)

! ----------------------------------------------------------------------
! Open MTF file or Open ZOU file
! In the case of MTF, is .extra is present read this file and set extra
! ----------------------------------------------------------------------

       else if (loword(wParam).eq.104.or.loword(wParam).eq.109.or.loword(wParam).eq.120) then

          if (loword(wParam).eq.109) then
             zou=.true.
             resid=.false.
             filter="Zinc Output File (*.zou)\0*.zou\0"C
          else if (loword(wParam).eq.104) then
             zou=.false.
             resid=.false.
             filter="Mesh Text File (*.mtf)\0*.mtf\0"C
          else
             zou=.true.
             resid=.true.
             filter="Residual File (*.resid)\0*.resid\0"C
          endif

          PathName=""C

          ofn%lStructSize = sizeof(ofn)
          ofn%hwndOwner = GetForegroundWindow()
          ofn%hInstance = NULL

          ofn%lpstrFilter = LOC(filter)

          ofn%lpstrCustomFilter = NULL
          ofn%nMaxCustFilter = 0
          ofn%nFilterIndex=1

          ofn%lpstrFile = LOC(PathName)
          ofn%nMaxFile = sizeof(PathName)

          ofn%lpstrFileTitle = NULL
          ofn%nMaxFileTitle = 0
          ofn%lpstrInitialDir=NULL
          ofn%lpstrTitle=loc("Select file"C)

          ofn%Flags = OFN_PATHMUSTEXIST

          ofn%nFileOffset=0
          ofn%nFileExtension=0

          ofn%lpstrDefExt = LOC("txt"C)

          ofn%lCustData=0
          ofn%lpfnHook=NULL
          ofn%lpTemplateName=NULL

! ----------------------------------------------------------------------
! remove the \0 from PathName. For ZOU replace .ZOU with .MTF for now
! since first thing is to read the MTF file.
! ----------------------------------------------------------------------

          if (GetOpenFileName(ofn).eq.TRUE) then

             l=index(PathName,char(0))
             PathName=PathName(:l-1)

             l=index(PathName,'.',back=.true.)
             extraName=PathName(:l-1)//'.extra'

             if (zou) then
                PathName=PathName(:l-1)//'.mtf'
             endif

             if (.not.zou) then
                message='Opening file:'
                ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)
                
                message=PathName
                ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)
             endif

! ----------------------------------------------------------------------
! Read .extra file if it exists. build up colour table
! ----------------------------------------------------------------------

             nextra=0           ! no of extra point sets
             nexcols=0          ! no of unique colours in extra list
             extra=.false.

             if (.not.zou) then
                open (1,file=extraName,status='old',iostat=ios)

                if (ios.eq.0) then

                   call winmess('Opening extra file:')
                   call winmess(trim(extraName))
                   extra=.true.
                   showextra=.true.

                   hMenu=GetMenu(hwnd)
                   ret=CheckMenuItem(hMenu,112,MF_CHECKED)

                   hMenu=GetMenu(hwnd)
                   ret=CheckMenuItem(hMenu,113,MF_CHECKED)

                   do
                      read (1,*,iostat=ios) type
                      if (ios.ne.0) exit

                      backspace(1)
                      if (type.eq.1) then
                         read (1,*,iostat=ios) type,cdum,cdum,cdum,r,g,b
                      else
                         read (1,*,iostat=ios) type,(cdum,i=1,12),r,g,b
                      endif

                      if (ios.ne.0) then
                         call winerror('Incomprehensible line in extra file')
                         return
                      endif

                      nextra=nextra+1
                      val=rgb(r,g,b)
                      
                      found=.false.
                      do i=1,nexcols
                         if (excols(i).eq.val) then
                            found=.true.
                            exit
                         endif
                      enddo
                      
                      if (.not.found) then
                         nexcols=nexcols+1
                         if (nexcols.gt.100) then
                            call winerror('Too many unique colours in .extra file')
                            return
                         endif
                         excols(nexcols)=val
                      endif
                   enddo
                   
                   rewind(1)
                   
                   do i=1,nexcols
                      extrabrush(i)=CreateSolidBrush(excols(i))
                   enddo
                                      
                   if (allocated(rextra)) deallocate(rextra,rextraold,extraind,extratype,extraview)
                   allocate (rextra(nextra,4,3),rextraold(nextra,4,3),extraind(nextra),extratype(nextra),extraview(nextra))

                   write (line,*) nextra
                   call winmess(trim(line))

                   nview=0
                   do i=1,nextra
                      read (1,'(a)') line
                      read (line,*) type

                      if (type.eq.1) then

                         read (line,*,iostat=ios) extratype(i),rextra(i,1,:),r,g,b,extraview(i)

                         if (ios.ne.0) then
                            read (line,*,iostat=ios) extratype(i),rextra(i,1,:),r,g,b
                            extraview(i)=1
                         endif

                      else
                         read (line,*,iostat=ios) extratype(i),((rextra(i,j,k),k=1,3),j=1,4),r,g,b,extraview(i)

                         if (ios.ne.0) then
                            read (line,*,iostat=ios) extratype(i),((rextra(i,j,k),k=1,3),j=1,4),r,g,b
                            extraview(i)=1
                         endif

                      endif

                      if (ios.ne.0) then
                         call winerror('Incomprehensible line in extra file')
                         return
                      endif

                      nview=max(nview,extraview(i))
                      
                      val=rgb(r,g,b)
                      do j=1,nexcols
                         if (excols(j).eq.val) then
                            extraind(i)=j
                            exit
                         endif
                      enddo

                   enddo

                   rextraold=rextra
                   rextra(:,:,1)=rextra(:,:,1)*xscale
                   rextra(:,:,2)=rextra(:,:,2)*yscale
                   rextra(:,:,3)=rextra(:,:,3)*zscale
                endif
             endif

! ----------------------------------------------------------------------
! Read from MTF file
! ----------------------------------------------------------------------
          
             open (1,file=PathName,status='old',iostat=ios)
       
             if (ios.ne.0) then
                call winerror('Could not open MTF file'//trim(PathName))
                return
             endif

             read (1,*)
             read (1,*,iostat=ioss(1)) cdum,imax
             read (1,*,iostat=ioss(2)) cdum,jmax
             read (1,*,iostat=ioss(3)) cdum,kmax
             
             if (ioss(1).ne.0) then
                call winerror('Incomprehensible imax')
                return
             endif

             if (ioss(2).ne.0) then
                call winerror('Incomprehensible jmax')
                return
             endif

             if (ioss(3).ne.0) then
                call winerror('Incomprehensible kmax')
                return
             endif

             read (1,*)
             read (1,*)
             read (1,*)

             if (allocated(rnode)) call free
             call alloc                               ! takes account of nextra above

             nreg=0
             nregnode=0

             do k=0,kmax
                do j=0,jmax
                   do i=0,imax
                      
                      read (1,*,iostat=ios) ic,jc,kc,x,y,z,mreg,mregup
                      
                      if (ios.ne.0) then
                         call winerror('Incomprehensible line in MTF file')
                         return
                      endif

                      if (.not.(ic.eq.i.and.jc.eq.j.and.kc.eq.k)) then
                         call winerror('Error in .mtf file')
                         return
                      endif
                
                      nreg=max(nreg,mregup)
                      nregnode=max(nregnode,mreg)

                      rnode(i,j,k,1)=x
                      rnode(i,j,k,2)=y
                      rnode(i,j,k,3)=z
                      iregnd(i,j,k)=mreg
                      iregup(i,j,k)=mregup
                      
                   enddo
                enddo
             enddo
       
             close (1)

! ----------------------------------------------------------------------
! For ZOU, ALSO load ZOU or RESID file
! Pathname is a local so we can play with it without danger
! ----------------------------------------------------------------------

             if (zou) then

! ----------------------------------------------------------------------
! If .ZIN file exists and there is a labels command, load the labels
! Also look for zou_format
! ----------------------------------------------------------------------

                l=index(PathName,'.',back=.true.)
                PathName=PathName(:l-1)//'.zin'

                zou_format='BINARY'
                nvar=0

                open (1,file=PathName,status='old',iostat=ios)

                if (ios.ne.0) then
                   call winerror('Could not open ZIN file')
                   return
                endif

                do
                   read (1,'(a)',iostat=ios) line
                   if (ios.ne.0) exit
                   
                   i=index(line,'=')
                   
                   if (i.ne.0) then
                      left= toupper(adjustl(line(:i-1)))
                      right=toupper(adjustl(line(i+1:)))
                      
                      if (left(1:6).eq.'LABELS') then
                         if (nvar==0) then
                            call winerror('NVAR must be set before LABELS')
                            return
                         endif
                         
                         read (right,*,iostat=ios) (label(ii),ii=1,nvar)
                         
                         if (ios.ne.0) then
                            call winerror('Incomprehensible LABELS line in ZIN file')
                            close (1)
                            return
                         endif
                      else if (left(1:10).eq.'ZOU_FORMAT') then
                         zou_format=right
                         if (zou_format.ne.'BINARY'.and.zou_format.ne.'ASCII') then
                            call winerror('zou_format in ZIN must be BINARY or ASCII')
                            close (1)
                            return
                         endif
                      else if (left(1:4).eq.'NVAR') then
                         read (right,*,iostat=ios) nvar
                         if (ios.ne.0) then
                            call winerror('Incomprehensible NVAR specification in ZIN file')
                            close (1)
                            return
                         endif
                         
                         if (allocated(label)) deallocate(label)
                         allocate (label(nvar))
                         
                         do ii=1,nvar
                            write (label(ii),'(a,i1)') 'u',ii
                         enddo
                         
                      endif
                   endif
                   
                enddo

                close (1)

                if (resid) nvar=nvar*3
                         
                if (nvar.eq.0) then
                   call winerror('Could not find NVAR in ZIN file')
                   return
                endif

! ----------------------------------------------------------------------
! Open ZOU or RESID file
! ----------------------------------------------------------------------

                l=index(PathName,'.',back=.true.)
                if (resid) then
                   PathName=PathName(:l-1)//'.resid'
                else
                   PathName=PathName(:l-1)//'.zou'
                endif

                call winmess('Opening file:')
                call winmess(PathName)

                if (zou_format=='BINARY') then
                   open (1,file=PathName,status='old',iostat=ios,form='unformatted')
                else
                   open (1,file=PathName,status='old',iostat=ios,form='formatted')
                endif

                if (ios.ne.0) then
                   call winerror('Could not open ZOU file'//trim(PathName))
                   return
                endif

! ----------------------------------------------------------------------
! Search for final
! ----------------------------------------------------------------------

                if (zou_format=='BINARY') then
                   read (1)
                   do
                      read (1,iostat=ios) i
                      if (ios /= 0) then
                         call winerror('Unexpected end of file (binary file)')
                         loaded=.false.
                         close (1)
                         return
                      endif

                      if (i==-2) exit
                   enddo
                else
                   do
                      read (1,'(a)',iostat=ios) line
                      if (ios.ne.0) then
                         call winerror('Unexpected end of file, searching for "final"')
                         close (1)
                         return
                      endif
                      
                      if (index(line,'final').ne.0) exit
                   enddo
                endif

! ----------------------------------------------------------------------
! Allocate memory
! ----------------------------------------------------------------------

                if (allocated(u)) deallocate(u,umax,umin,uspan,umaxUser,uminUser,uspanUser)
                allocate(u(0:imax,0:jmax,0:kmax,nvar),umax(nvar),umin(nvar),uspan(nvar),&
                     umaxUser(nvar),uminUser(nvar),uspanUser(nvar),stat=allocerr)

                if (allocerr.ne.0) then
                   call winerror('Failed to allocate arrays for ZOU or RESID display')
                   close (1)
                   return
                endif

! ----------------------------------------------------------------------
! Read in ZOU or RESID
! ----------------------------------------------------------------------

                do i=0,imax
                   do j=0,jmax
                      do k=0,kmax
                         do ii=1,nvar

                            if (zou_format=='BINARY') then
                               read (1,iostat=ios) i1,j1,k1,ii1,u1
                            else
                               read (1,*,iostat=ios) i1,j1,k1,ii1,u1
                            endif

                            if (ios.ne.0) then
                               call winerror('Incomprehensible number in .ZOU or .RESID file')
                               close (1)
                               return
                            endif

                            if (.not.(i.eq.i1.and.j.eq.j1.and.k.eq.k1.and.ii.eq.ii1)) then
                               call winerror('ZOU file out of sequence: vary k first, then j, then i')
                               close (1)
                               return
                            endif

                            u(i,j,k,ii)=u1
                            
                            if (i.eq.0.and.j.eq.0.and.k.eq.0) then
                               umax(ii)=u1
                               umin(ii)=u1
                            else if (u1.gt.umax(ii)) then
                               umax(ii)=u1
                            else if (u1.lt.umin(ii)) then
                               umin(ii)=u1
                            endif
                            
                         enddo
                      enddo
                   enddo
                enddo

                uspan=umax-umin
                do ii=1,nvar
                   if (uspan(ii)==0) uspan(ii)=1
                enddo

                umaxUser=umax
                uminUser=umin
                uspanUser=uspan
     
                close (1)

             endif

! ----------------------------------------------------------------------
! Set up cols and nodecols. Initialise view
! ----------------------------------------------------------------------

             if (allocated(cols)) call freekey
             call allockey

             call initview(hwnd)

             hMenu=GetMenu(hwnd)
             ret=EnableMenuItem(hMenu,106,MF_GRAYED)
             ret=EnableMenuItem(hMenu,107,MF_GRAYED)

             if (extra) then
                ret=EnableMenuItem(hMenu,112,MF_ENABLED)
             else
                ret=EnableMenuItem(hMenu,112,MF_GRAYED)
             endif
          endif

       endif           ! switch for command selected

       return

! ----------------------------------------------------------------------
! WM_MOUSEMOVE
! ----------------------------------------------------------------------

    else if (iMsg.eq.WM_MOUSEMOVE) then

       if (.not.loaded) return

       if (iand(wParam,MK_LBUTTON).ne.0) then
          if (iand(wParam,MK_SHIFT).eq.0) then
             xmouse2=loword(lParam)
             ymouse2=hiword(lParam)
             
             call rotz( (xmouse2-xmouse)*ascale)
             call roty(-(ymouse2-ymouse)*ascale)
             
             xmouse=xmouse2
             ymouse=ymouse2

             ret=InvalidateRect(hwnd,NULL,FALSE)
          else
             xmouse2=loword(lParam)
             ymouse2=hiword(lParam)

             call shift((xmouse2-xmouse)*sscale/scale,(ymouse-ymouse2)*sscale/scale)

             xmouse=xmouse2
             ymouse=ymouse2

             ret=InvalidateRect(hwnd,NULL,FALSE)
          endif
       endif

       return

    else if (iMsg.eq.WM_KEYDOWN) then

! ----------------------------------------------------------------------
! WM_KEYDOWN
! Need fastdraw switch for up/down when in slice>0
! ----------------------------------------------------------------------

       if (.not.loaded) return

       if (wParam.eq.VK_LEFT) then
          call rotz(-angz)

       else if (wParam.eq.VK_RIGHT) then
          call rotz(angz)

       else if (wParam.eq.VK_UP) then
          if (slice.eq.0) then
             call roty(angy)
          else if (slice.eq.1.and.icutmax.lt.imax) then
             icutmax=icutmax+1
             icutmin=icutmin+1

          else if (slice.eq.2.and.jcutmax.lt.jmax) then
             jcutmax=jcutmax+1
             jcutmin=jcutmin+1

          else if (slice.eq.3.and.kcutmax.lt.kmax) then
             kcutmax=kcutmax+1
             kcutmin=kcutmin+1
          else
             return
          endif

       else if (wParam.eq.VK_DOWN) then
          if (slice.eq.0) then
             call roty(-angy)
          else if (slice.eq.1.and.icutmin.gt.0) then
             icutmax=icutmax-1
             icutmin=icutmin-1

          else if (slice.eq.2.and.jcutmin.gt.0) then
             jcutmax=jcutmax-1
             jcutmin=jcutmin-1

          else if (slice.eq.3.and.kcutmin.gt.0) then
             kcutmax=kcutmax-1
             kcutmin=kcutmin-1
          else
             return
          endif

       else if (char(wParam).eq.'M') then
          call rotx(angx)

       else if (char(wParam).eq.'N') then
          call rotx(-angx)

       else if (char(wParam).eq.'A') then
!          scale=scale+speed
          scale=scale*1.1

       else if (char(wParam).eq.'Z') then
!          scale=scale-speed
          scale=scale*0.9

       else
          return
       endif

       ret=InvalidateRect(hwnd,NULL,FALSE)
       return

    else if (iMsg.eq.WM_PAINT) then

! ----------------------------------------------------------------------
! WM_PAINT
! ----------------------------------------------------------------------

       hdc=BeginPaint(hwnd,ps)

       if (loaded) then
          call drawit(hdcMem)
          ret=BitBlt(hdc,0,0,width,height,hdcMem,0,0,SRCCOPY)
       endif

       ret=EndPaint(hwnd,ps)

       return
       
    else if (iMsg.eq.WM_DESTROY) then

! ----------------------------------------------------------------------
! WM_DESTROY and DefWindowProc
! ----------------------------------------------------------------------

       if (allocated(rnode)) then
          call free
          call freekey
       endif

       ret=DeleteObject(hbitMem)
       ret=DeleteDC(hdcMem)
       ret=DeleteObject(btnoff)
       ret=DeleteObject(btnon)

       call PostQuitMessage(0)
       return
    else
       WndProc=DefWindowProc(hwnd, iMsg, wParam, lParam)
       return
    endif

  end function WndProc

! ######################################################################

  subroutine rescale(w,h)
    use global, only : scale,xc,yc,redge
    integer w,h
    
    scale=0.8*min(w,h)/max(redge(2,2)-redge(2,1),redge(3,2)-redge(3,1))
    
    xc=w/2
    yc=h/2
  end subroutine rescale
  
! ######################################################################

  subroutine drawit(hdc)

! ----------------------------------------------------------------------
! Draw the scene on hdc. First, the polygons
! In case of fastdraw, do not call getpoly but reorder existing lists
! ----------------------------------------------------------------------

    use user32
    use gdi32, Polygon => MSFWIN$Polygon, Rectangle => MSFWIN$Rectangle,&
         LineTo => MSFWIN$LineTo, Ellipse => MSFWIN$Ellipse

    use const
    use global
    use heap

    integer ip,ireg,i,oldreg,ret,xint,yint,ibrush,nshow,iS
    integer(HANDLE) hdc
    character(100) string
    integer iedge,x1int,y1int,x2int,y2int,ic,ivar2,ivar3
    real xav,xavmin,lowedge(3),lowface(3)
    type(T_POINT) pline(4),oldpoint

    real bht
    integer ytop,ybot

! ----------------------------------------------------------------------
! Clearscreen
! rear parts of bounding box and text
! ----------------------------------------------------------------------

    ret=SelectObject(hdc,GetStockObject(SYSTEM_FIXED_FONT))      
    ret=SetBkMode(hdc,TRANSPARENT)

    ret=SelectObject(hdc,GetStockObject(LTGRAY_BRUSH))
    ret=Rectangle(hdc,0,0,width,height)
    
    ret=SetTextAlign(hdc,TA_CENTER)

    do ic=1,3
       do iedge=1,4
          xav=(bb(ic,iedge,1,1)+bb(ic,iedge,2,1))/2
          
          if (iedge.eq.1) then
             xavmin=xav
             lowedge(ic)=iedge
          else if (xav.lt.xavmin) then
             xavmin=xav
             lowedge(ic)=iedge
          endif
       enddo

       x1int= bb(ic,lowedge(ic),1,2)*scale+xc
       y1int=-bb(ic,lowedge(ic),1,3)*scale+yc
       x2int= bb(ic,lowedge(ic),2,2)*scale+xc
       y2int=-bb(ic,lowedge(ic),2,3)*scale+yc
       
       ret=MoveToEx(hdc,x1int,y1int,oldpoint)
       ret=LineTo(hdc,x2int,y2int)

       if (ff(ic,1,1).lt.ff(ic,2,1)) then
          lowface(ic)=1
       else
          lowface(ic)=2
       endif
       
       x1int= ff(ic,lowface(ic),2)*scale+xc
       y1int=-ff(ic,lowface(ic),3)*scale+yc
       
       string=ffstr(ic,lowface(ic))

       ret=TextOut(hdc,x1int,y1int,string,len_trim(string))
    enddo

! ----------------------------------------------------------------------
! Draw the polygons or ellipses
! ----------------------------------------------------------------------

    if (showgridlines) then
       ret=SelectObject(hdc,GetStockObject(BLACK_PEN))
    else
       ret=SelectObject(hdc,GetStockObject(NULL_PEN))
    endif

    if ((.not.zou).and.(.not.extra)) then

       if (.not.fastdraw) call getpoly
       call heapsort(pxcent,npoly)
       
       do ip=1,npoly
          ireg=preg(ip)
          
          ibrush=abs(ireg)
          
          if (ip.eq.1) then
             if (ireg.lt.0) then
                ret=SelectObject(hdc,nodecols(ibrush))
             else
                ret=SelectObject(hdc,cols(ibrush))
             endif
          else if (ireg.ne.oldreg) then
             if (ireg.lt.0) then
                ret=SelectObject(hdc,nodecols(ibrush))
             else
                ret=SelectObject(hdc,cols(ibrush))
             endif
          endif
          
          oldreg=ireg

          if (ireg.lt.0) then
             xint= plist(ip,1,2)*scale+xc
             yint=-plist(ip,1,3)*scale+yc
             
             ret=Ellipse(hdc,xint-5,yint-5,xint+5,yint+5)
          else
             do i=1,4
                pline(i)%x= plist(ip,i,2)*scale+xc
                pline(i)%y=-plist(ip,i,3)*scale+yc
             enddo
             
             ret=Polygon(hdc,pline(1),4)
          endif
       enddo
       
! ----------------------------------------------------------------------
! In case of zou. If |preg(ip)|<= nreg, use key else colour bar legend
! Eg 3 regions, but preg(ip)=4, use legend(0) (4-3-1=0)
! ----------------------------------------------------------------------

    else if (zou) then

       if (.not.fastdraw) call getpolyzou
       call heapsort(pxcent,npoly)

       do ip=1,npoly
          ireg=preg(ip)
          
          ibrush=abs(ireg)
          
          if (ireg.lt.0) then
             if (ibrush.le.nregnode) then
                ret=SelectObject(hdc,nodecols(ibrush))
             else
                ret=SelectObject(hdc,legend(ibrush-nregnode-1))
             endif
             
             xint= plist(ip,1,2)*scale+xc
             yint=-plist(ip,1,3)*scale+yc
             
             ret=Ellipse(hdc,xint-5,yint-5,xint+5,yint+5)
          else
             if (ibrush.le.nreg) then
                ret=SelectObject(hdc,cols(ibrush))
             else
                ret=SelectObject(hdc,legend(ibrush-nreg-1))
             endif
             
             do i=1,4
                pline(i)%x= plist(ip,i,2)*scale+xc
                pline(i)%y=-plist(ip,i,3)*scale+yc
             enddo
             
             ret=Polygon(hdc,pline(1),4)
          endif
       enddo

! ----------------------------------------------------------------------
! case of extra (same as zou code but extrabrush instead of legend)
! ----------------------------------------------------------------------

    else if (extra) then

       if (.not.fastdraw) call getpoly                      ! knows about extra
       call heapsort(pxcent,npoly)

       do ip=1,npoly
          ireg=preg(ip)
          
          ibrush=abs(ireg)
          
          if (ireg.lt.0) then
             if (ibrush.le.nregnode) then
                ret=SelectObject(hdc,nodecols(ibrush))
             else
                ret=SelectObject(hdc,extrabrush(ibrush-nregnode))
             endif
             
             xint= plist(ip,1,2)*scale+xc
             yint=-plist(ip,1,3)*scale+yc
             
             ret=Ellipse(hdc,xint-5,yint-5,xint+5,yint+5)
          else
             if (ibrush.le.nreg) then
                ret=SelectObject(hdc,cols(ibrush))
             else
                ret=SelectObject(hdc,extrabrush(ibrush-nreg))
             endif
             
             do i=1,4
                pline(i)%x= plist(ip,i,2)*scale+xc
                pline(i)%y=-plist(ip,i,3)*scale+yc
             enddo
             
             ret=Polygon(hdc,pline(1),4)
          endif
       enddo

    endif

    ret=SelectObject(hdc,GetStockObject(BLACK_PEN))

! ----------------------------------------------------------------------
! Front part of bounding box and text
! ----------------------------------------------------------------------

    ret=SetTextAlign(hdc,TA_CENTER)
    
    do ic=1,3
       do iedge=1,4
          if (iedge.ne.lowedge(ic)) then
             
             x1int= bb(ic,iedge,1,2)*scale+xc
             y1int=-bb(ic,iedge,1,3)*scale+yc
             x2int= bb(ic,iedge,2,2)*scale+xc
             y2int=-bb(ic,iedge,2,3)*scale+yc
             
             ret=MoveToEx(hdc,x1int,y1int,oldpoint)
             ret=LineTo(hdc,x2int,y2int)
          endif
       enddo
       
       x1int= ff(ic,3-lowface(ic),2)*scale+xc
       y1int=-ff(ic,3-lowface(ic),3)*scale+yc
       
       string=ffstr(ic,3-lowface(ic))
       
       ret=TextOut(hdc,x1int,y1int,string,len_trim(string))
    enddo

! ----------------------------------------------------------------------
! draw the key
! ----------------------------------------------------------------------

    ret=SetTextAlign(hdc,ior(TA_LEFT,TA_TOP))

    nshow=min(nreg-ScrollEl+1,maxbutton)
    
    do i=1,nshow
       xint=width-keyxoff-keywidth
       yint=keyyoff+(i-1)*keypitch

       iS=i+ScrollEl-1
       
       ret=SelectObject(hdc,cols(iS))
       ret=Rectangle(hdc,xint,yint,xint+keywidth,yint+keyheight)
       
       if (showreg(iS).eq.0) then
          ret=MoveToEx(hdc,xint,yint,oldpoint)
          ret=LineTo(hdc,xint+keywidth,yint+keyheight)
          
          ret=MoveToEx(hdc,xint+keywidth,yint,oldpoint)
          ret=LineTo(hdc,xint,yint+keyheight)
       endif
       
       if (showreg(iS).eq.2) then
          write (string,'(a,i3)') 'DATA',iS
       else
          write (string,'(a,i3)') 'Elem',iS
       endif

       ret=TextOut(hdc,xint+2,yint,string,7)
    enddo

    if (nreg.gt.maxbutton) then
       ret=SelectObject(hdc,btnon)
       
       pline(1)%x=width-keyxoff-keywidth
       pline(1)%y=keyyoff/2
       
       pline(2)%x=width-keyxoff
       pline(2)%y=keyyoff/2
       
       pline(3)%x=width-keyxoff-keywidth/2
       pline(3)%y=0
       
       ret=Polygon(hdc,pline(1),3)
       
       pline(1)%y=keyyoff*1.5+(maxbutton-1)*keypitch+keyheight
       pline(2)%y=pline(1)%y
       pline(3)%y=pline(1)%y+keyyoff/2
       
       ret=Polygon(hdc,pline(1),3)
    endif

! ----------------------------------------------------------------------
! node key
! ----------------------------------------------------------------------

    nshow=min(nregnode-ScrollNd+1,maxbutton)

    do i=1,nshow
       xint=width-2*keyxoff-2*keywidth
       yint=keyyoff+(i-1)*keypitch
       
       iS=i+ScrollNd-1

       ret=SelectObject(hdc,nodecols(iS))
       ret=Rectangle(hdc,xint,yint,xint+keywidth,yint+keyheight)
       
       if (nodeshowreg(iS).eq.0) then
          ret=MoveToEx(hdc,xint,yint,oldpoint)
          ret=LineTo(hdc,xint+keywidth,yint+keyheight)
          
          ret=MoveToEx(hdc,xint+keywidth,yint,oldpoint)
          ret=LineTo(hdc,xint,yint+keyheight)
       endif
       
       if (nodeshowreg(iS).eq.2) then
          write (string,'(a,i3)') 'DATA',iS
       else
          write (string,'(a,i3)') 'Node',iS
       endif

       ret=TextOut(hdc,xint+2,yint,string,7)
    enddo

    if (nregnode.gt.maxbutton) then
       ret=SelectObject(hdc,btnon)
       
       pline(1)%x=width-2*keyxoff-2*keywidth
       pline(1)%y=keyyoff/2
       
       pline(2)%x=pline(1)%x+keywidth
       pline(2)%y=keyyoff/2
       
       pline(3)%x=pline(1)%x+keywidth/2
       pline(3)%y=0
       
       ret=Polygon(hdc,pline(1),3)
       
       pline(1)%y=keyyoff*1.5+(maxbutton-1)*keypitch+keyheight
       pline(2)%y=pline(1)%y
       pline(3)%y=pline(1)%y+keyyoff/2
       
       ret=Polygon(hdc,pline(1),3)
    endif

! ----------------------------------------------------------------------
! Draw the buttons
! ----------------------------------------------------------------------

    do i=1,nbtn
       xint=btnxoff+(i-1)*btnpitch
       yint=height-btnyoff-btnheight
         
       if (btnstate(i).eq.0) then
          ret=SelectObject(hdc,btnoff)
       else
          ret=SelectObject(hdc,btnon)
       endif
       
       ret=Rectangle(hdc,xint,yint,xint+btnwidth,yint+btnheight)
       
       ret=TextOut(hdc,xint+btntxtpos(i,1),yint+btntxtpos(i,2),&
            btntext(i),len_trim(btntext(i)))
    enddo
    
! ----------------------------------------------------------------------
! In case of ZOU, Draw the colour legend
! ----------------------------------------------------------------------

    if (zou) then
!       ret=SetBkMode(hdc,OPAQUE)
       ret=SetTextAlign(hdc,ior(TA_LEFT,TA_BOTTOM))
       ret=SelectObject(hdc,GetStockObject(NULL_PEN))
       bht=(height-legendyoff*2)/63.0
       do i=0,63
          ytop=height-legendyoff-i*bht-bht
          ybot=height-legendyoff-i*bht+1

          ret=SelectObject(hdc,legend(i*4))
          ret=Rectangle(hdc,legendxoff,ytop,legendxoff+legendwidth,ybot)
       enddo
       ret=SelectObject(hdc,GetStockObject(BLACK_PEN))
       
       write (string,'(g13.5)') umaxUser(ivar)
       string=adjustl(string)
       ret=TextOut(hdc,legendxoff+legendwidth+5,int(legendyoff-bht+6),string,len_trim(string))

       write (string,'(g13.5)') uminUser(ivar)
       string=adjustl(string)
       ret=TextOut(hdc,legendxoff+legendwidth+5,height-legendyoff+5,string,len_trim(string))

       ret=SelectObject(hdc,btnon)
       ret=Rectangle(hdc,legendxoff+legendwidth+5, legendyoff, &
                         legendxoff+legendwidth+68,legendyoff+20)

       ret=Rectangle(hdc,legendxoff+legendwidth+5, height-legendyoff-35, &
                         legendxoff+legendwidth+68,height-legendyoff-15)

       ret=TextOut(hdc,legendxoff+legendwidth+7,legendyoff+17,       'set max',7)
       ret=TextOut(hdc,legendxoff+legendwidth+7,height-legendyoff-17,'set min',7)

       ret=SetTextAlign(hdc,ior(TA_LEFT,TA_TOP))

       write (string,'(a,i2,a,i2)') 'Showing variable',ivar,' of',nvar
       ret=TextOut(hdc,10,10,string,len_trim(string))

       if (resid) then
          ivar2=(ivar-1)/3+1   ! 1,2,3 -> 1; 4,5,6 -> 2
          ivar3=ivar-(ivar2-1)*3
          if (ivar3.eq.1) then
             write (string,'(a)') 'Variable name: '//trim(label(ivar2))//'.lhs'
          else if (ivar3.eq.2) then
             write (string,'(a)') 'Variable name: '//trim(label(ivar2))//'.rhs'
          else
             write (string,'(a)') 'Variable name: '//trim(label(ivar2))//'.resid'
          endif

       else
          write (string,'(a)') 'Variable name: '//trim(label(ivar))
       endif

       ret=TextOut(hdc,10,30,string,len_trim(string))

       if (nvar.gt.1) then

          ret=SelectObject(hdc,btnon)
          ret=Rectangle(hdc,xvar,yvar,xvar+varsize,yvar+varsize)
          ret=Rectangle(hdc,xvar+varsize+varsep,yvar,xvar+varsize*2+varsep,yvar+varsize)
          
          ret=TextOut(hdc,xvar+vartxt+1,yvar+vartxt,'<',1)
          ret=TextOut(hdc,xvar+varsize+varsep+vartxt+1,yvar+vartxt,'>',1)
       endif
    endif

! ----------------------------------------------------------------------
! In case of EXTRA state which view is being seen and show < > buttons
! ----------------------------------------------------------------------

    if (extra.and.showextra) then
 
       ret=SetTextAlign(hdc,ior(TA_LEFT,TA_TOP))

       write (string,'(a,i2,a,i2)') 'Showing EXT view',iview,' of',nview
       ret=TextOut(hdc,10,10,string,len_trim(string))

       if (nview.gt.1) then

          ret=SelectObject(hdc,btnon)
          ret=Rectangle(hdc,xvar,yvar,xvar+varsize,yvar+varsize)
          ret=Rectangle(hdc,xvar+varsize+varsep,yvar,xvar+varsize*2+varsep,yvar+varsize)
          
          ret=TextOut(hdc,xvar+vartxt,yvar+vartxt,'<',1)
          ret=TextOut(hdc,xvar+varsize+varsep+vartxt,yvar+vartxt,'>',1)
       endif

    endif

  end subroutine drawit

! ######################################################################

  subroutine setbb

! ----------------------------------------------------------------------
! bounding box setup. Also sets text for box size
! ----------------------------------------------------------------------

    use global

    real x1,x2,y1,y2,z1,z2
    integer i

    x1=redge(1,1)*xscale
    x2=redge(1,2)*xscale
    
    y1=redge(2,1)*yscale
    y2=redge(2,2)*yscale
    
    z1=redge(3,1)*zscale
    z2=redge(3,2)*zscale
      
    bb(1,1,1,:)=[x1,y1,z1]
    bb(1,1,2,:)=[x2,y1,z1]
      
    bb(1,2,1,:)=[x1,y2,z1]
    bb(1,2,2,:)=[x2,y2,z1]
      
    bb(1,3,1,:)=[x1,y1,z2]
    bb(1,3,2,:)=[x2,y1,z2]
    
    bb(1,4,1,:)=[x1,y2,z2]
    bb(1,4,2,:)=[x2,y2,z2]
    
      
    bb(2,1,1,:)=[x1,y1,z1]
    bb(2,1,2,:)=[x1,y2,z1]
    
    bb(2,2,1,:)=[x2,y1,z1]
    bb(2,2,2,:)=[x2,y2,z1]
    
    bb(2,3,1,:)=[x1,y1,z2]
    bb(2,3,2,:)=[x1,y2,z2]
    
    bb(2,4,1,:)=[x2,y1,z2]
    bb(2,4,2,:)=[x2,y2,z2]
      
      
    bb(3,1,1,:)=[x1,y1,z1]
    bb(3,1,2,:)=[x1,y1,z2]
    
    bb(3,2,1,:)=[x2,y1,z1]
    bb(3,2,2,:)=[x2,y1,z2]
    
    bb(3,3,1,:)=[x1,y2,z1]
    bb(3,3,2,:)=[x1,y2,z2]
    
    bb(3,4,1,:)=[x2,y2,z1]
    bb(3,4,2,:)=[x2,y2,z2]
    
    ff(1,1,:)=[x1,(y1+y2)/2,(z1+z2)/2]
    ff(1,2,:)=[x2,(y1+y2)/2,(z1+z2)/2]
    
    ff(2,1,:)=[(x1+x2)/2,y1,(z1+z2)/2]
    ff(2,2,:)=[(x1+x2)/2,y2,(z1+z2)/2]
    
    ff(3,1,:)=[(x1+x2)/2,(y1+y2)/2,z1]
    ff(3,2,:)=[(x1+x2)/2,(y1+y2)/2,z2]
      
    do i=1,2
       write (ffstr(1,i),'(a,g13.5)') 'x=',ff(1,i,1)/xscale
       write (ffstr(2,i),'(a,g13.5)') 'y=',ff(2,i,2)/yscale
       write (ffstr(3,i),'(a,g13.5)') 'z=',ff(3,i,3)/zscale
    enddo
  end subroutine setbb

! ######################################################################

!  print *,ibits(i,16,16)   ! HIWORD
!  print *,ibits(i,0,16)    ! LOWORD

  function hiword(i)
    integer hiword,i
    hiword=ibits(i,16,16)
  end function hiword

  function loword(i)
    integer loword,i
    loword=ibits(i,0,16)
  end function loword

  function rgb(r,g,b)
    integer rgb,r,g,b
    rgb=r+256*g+65536*b
  end function rgb

! ######################################################################

  subroutine getpoly

! ----------------------------------------------------------------------
! Go through pp and generate a list of polygons in plist.
! also store pxcent (x coord of poly) and preg (region number of poly)
! Set npoly as the number of polygons
! ----------------------------------------------------------------------

    use user32
    use const
    use global

    integer i,j,k,ireg1,ireg2,ireg1fnd,ireg2fnd,icol,ireg,ret

    npoly=0
    
! ----------------------------------------------------------------------
! First, the nodes
! ----------------------------------------------------------------------

    if (any(nodeshowreg.eq.1)) then

       do i=icutmin,icutmax
          do j=jcutmin,jcutmax
             do k=kcutmin,kcutmax
                
                ireg1=iregnd(i,j,k)
                
                if (nodeshowreg(ireg1).eq.1) then
                   npoly=npoly+1
                   
                   preg(npoly)=-ireg1
                   plist(npoly,1,:)=pp(i,j,k,:)
                   pxcent(npoly)=pp(i,j,k,1)
                endif
             enddo
          enddo
       enddo
    endif

! ----------------------------------------------------------------------
! i-directed faces
! ----------------------------------------------------------------------

    do i=icutmin,icutmax
       do j=jcutmin,jcutmax-1
          do k=kcutmin,kcutmax-1

             if (i.eq.icutmax) then
                ireg1=-1
             else
                ireg1=iregup(i,j,k)
             endif

             if (i.eq.icutmin) then
                ireg2=-1
             else
                ireg2=iregup(i-1,j,k)
             endif

             ireg1fnd=0
             ireg2fnd=0

             do ireg=1,nreg
                if (showreg(ireg).eq.1.and.ireg1.eq.ireg) then
                   ireg1fnd=1
                   icol=ireg1
                endif
                
                if (showreg(ireg).eq.1.and.ireg2.eq.ireg) then
                   ireg2fnd=1
                   icol=ireg2
                endif
             enddo
             
             if (ireg1fnd+ireg2fnd.eq.1) then
                
                npoly=npoly+1
                
                preg(npoly)=icol
                
                plist(npoly,1,:)=pp(i,j,k,:)
                plist(npoly,2,:)=pp(i,j,k+1,:)
                plist(npoly,3,:)=pp(i,j+1,k+1,:)
                plist(npoly,4,:)=pp(i,j+1,k,:)
                
                pxcent(npoly)=sum(plist(npoly,:,1))/4
             endif
          enddo
       enddo
    enddo

! ----------------------------------------------------------------------
! j-directed faces
! ----------------------------------------------------------------------

    do i=icutmin,icutmax-1
       do j=jcutmin,jcutmax
          do k=kcutmin,kcutmax-1

             if (j.eq.jcutmax) then
                ireg1=-1
             else
                ireg1=iregup(i,j,k)
             endif

             if (j.eq.jcutmin) then
                ireg2=-1
             else
                ireg2=iregup(i,j-1,k)
             endif

             ireg1fnd=0
             ireg2fnd=0

             do ireg=1,nreg
                if (showreg(ireg).eq.1.and.ireg1.eq.ireg) then
                   ireg1fnd=1
                   icol=ireg1
                endif
                
                if (showreg(ireg).eq.1.and.ireg2.eq.ireg) then
                   ireg2fnd=1
                   icol=ireg2
                endif
             enddo
             
             if (ireg1fnd+ireg2fnd.eq.1) then
                
                npoly=npoly+1
                
                preg(npoly)=icol
                
                plist(npoly,1,:)=pp(i,j,k,:)
                plist(npoly,2,:)=pp(i,j,k+1,:)
                plist(npoly,3,:)=pp(i+1,j,k+1,:)
                plist(npoly,4,:)=pp(i+1,j,k,:)
                 
                pxcent(npoly)=sum(plist(npoly,:,1))/4
             endif
          enddo
       enddo
    enddo

! ----------------------------------------------------------------------
! k-directed faces
! ----------------------------------------------------------------------

    do i=icutmin,icutmax-1
       do j=jcutmin,jcutmax-1
          do k=kcutmin,kcutmax

             if (k.eq.kcutmax) then
                ireg1=-1
             else
                ireg1=iregup(i,j,k)
             endif

             if (k.eq.kcutmin) then
                ireg2=-1
             else
                ireg2=iregup(i,j,k-1)
             endif

             ireg1fnd=0
             ireg2fnd=0

             do ireg=1,nreg
                if (showreg(ireg).eq.1.and.ireg1.eq.ireg) then
                   ireg1fnd=1
                   icol=ireg1
                endif
                
                if (showreg(ireg).eq.1.and.ireg2.eq.ireg) then
                   ireg2fnd=1
                   icol=ireg2
                endif
             enddo
             
             if (ireg1fnd+ireg2fnd.eq.1) then
                
                npoly=npoly+1
                
                preg(npoly)=icol

                plist(npoly,1,:)=pp(i,j,k,:)
                plist(npoly,2,:)=pp(i+1,j,k,:)
                plist(npoly,3,:)=pp(i+1,j+1,k,:)
                plist(npoly,4,:)=pp(i,j+1,k,:)
                
                pxcent(npoly)=sum(plist(npoly,:,1))/4
             endif
          enddo
       enddo
    enddo

! ----------------------------------------------------------------------
! Add extra elements if necessary
! preg is index number (1,2,...) of colour specified. We add this to
! nreg, nregnode. Eg nregnode=4, extra col index=1, preg=5
! ----------------------------------------------------------------------

    if (extra.and.showextra) then

       do i=1,nextra

          if (extraview(i).eq.iview) then

             npoly=npoly+1

             if (extratype(i).eq.1) then
                preg(npoly)=-(extraind(i)+nregnode)
                plist(npoly,1,:)=rextra(i,1,:)
                pxcent(npoly)=rextra(i,1,1)
             else
                preg(npoly)=extraind(i)+nreg
                plist(npoly,:,:)=rextra(i,:,:)
                pxcent(npoly)=sum(rextra(i,:,1))/4
             endif
          endif

       enddo

    endif

  end subroutine getpoly

! ######################################################################

  subroutine getpolyzou

! ----------------------------------------------------------------------
! Go through pp and generate a list of polygons in plist.
! also store pxcent (x coord of poly) and preg (region number of poly)
! Set npoly as the number of polygons
! ----------------------------------------------------------------------

    use user32
    use const
    use global

    integer i,j,k,ireg1,ireg2,ireg1fnd,ireg2fnd,icol,ireg,ret
    real uav

    npoly=0
    
! ----------------------------------------------------------------------
! First, the nodes
! ----------------------------------------------------------------------

    if (any(nodeshowreg.ne.0)) then

       do i=icutmin,icutmax
          do j=jcutmin,jcutmax
             do k=kcutmin,kcutmax
                
                ireg1=iregnd(i,j,k)
                
                if (nodeshowreg(ireg1).eq.1) then
                   npoly=npoly+1
                   
                   preg(npoly)=-ireg1
                   plist(npoly,1,:)=pp(i,j,k,:)
                   pxcent(npoly)=pp(i,j,k,1)
                else if (nodeshowreg(ireg1).eq.2) then
                   npoly=npoly+1

                   uav=u(i,j,k,ivar)
                   if (uav.gt.umaxUser(ivar).or.uav.lt.uminUser(ivar)) then
                      preg(npoly)=-(256+nregnode+1)
                   else
                      preg(npoly)=-(int((uav-uminUser(ivar))/uspanUser(ivar)*255)+nregnode+1)
                   endif

                   plist(npoly,1,:)=pp(i,j,k,:)
                   pxcent(npoly)=pp(i,j,k,1)
                endif
             enddo
          enddo
       enddo
    endif

! ----------------------------------------------------------------------
! i-directed faces
! ----------------------------------------------------------------------

    do i=icutmin,icutmax
       do j=jcutmin,jcutmax-1
          do k=kcutmin,kcutmax-1

             if (i.eq.icutmax) then
                ireg1=-1
             else
                ireg1=iregup(i,j,k)
             endif

             if (i.eq.icutmin) then
                ireg2=-1
             else
                ireg2=iregup(i-1,j,k)
             endif

             ireg1fnd=0
             ireg2fnd=0

! if (showreg(ireg1).eq.1) icol=ireg1?? would need showreg(-1)=0

             do ireg=1,nreg
                if (ireg1.eq.ireg) then
                   if (showreg(ireg).eq.1) then
                      ireg1fnd=1
                      icol=ireg1
                   else if (showreg(ireg).eq.2) then
                      ireg1fnd=1
                      icol=-1
                   endif
                endif

                if (ireg2.eq.ireg) then
                   if (showreg(ireg).eq.1) then
                      ireg2fnd=1
                      icol=ireg2
                   else if (showreg(ireg).eq.2) then
                      ireg2fnd=1
                      icol=-1
                   endif
                endif
             enddo
             
             if (ireg1fnd+ireg2fnd.eq.1) then
                
                npoly=npoly+1
                
                if (icol.eq.-1) then
                   uav=(u(i,j,k,ivar)+u(i,j,k+1,ivar)+u(i,j+1,k+1,ivar)+u(i,j+1,k,ivar))/4

                   if (uav.gt.umaxUser(ivar).or.uav.lt.uminUser(ivar)) then
                      preg(npoly)=256+nreg+1
                   else
                      preg(npoly)=int((uav-uminUser(ivar))/uspanUser(ivar)*255)+nreg+1
                   endif
                else
                   preg(npoly)=icol
                endif
                
                plist(npoly,1,:)=pp(i,j,k,:)
                plist(npoly,2,:)=pp(i,j,k+1,:)
                plist(npoly,3,:)=pp(i,j+1,k+1,:)
                plist(npoly,4,:)=pp(i,j+1,k,:)
                
                pxcent(npoly)=sum(plist(npoly,:,1))/4
             endif
          enddo
       enddo
    enddo

! ----------------------------------------------------------------------
! j-directed faces
! ----------------------------------------------------------------------

    do i=icutmin,icutmax-1
       do j=jcutmin,jcutmax
          do k=kcutmin,kcutmax-1

             if (j.eq.jcutmax) then
                ireg1=-1
             else
                ireg1=iregup(i,j,k)
             endif

             if (j.eq.jcutmin) then
                ireg2=-1
             else
                ireg2=iregup(i,j-1,k)
             endif

             ireg1fnd=0
             ireg2fnd=0

             do ireg=1,nreg
                if (ireg1.eq.ireg) then
                   if (showreg(ireg).eq.1) then
                      ireg1fnd=1
                      icol=ireg1
                   else if (showreg(ireg).eq.2) then
                      ireg1fnd=1
                      icol=-1
                   endif
                endif

                if (ireg2.eq.ireg) then
                   if (showreg(ireg).eq.1) then
                      ireg2fnd=1
                      icol=ireg2
                   else if (showreg(ireg).eq.2) then
                      ireg2fnd=1
                      icol=-1
                   endif
                endif
             enddo

             if (ireg1fnd+ireg2fnd.eq.1) then
                
                npoly=npoly+1
                
                if (icol.eq.-1) then
                   uav=(u(i,j,k,ivar)+u(i,j,k+1,ivar)+u(i+1,j,k+1,ivar)+u(i+1,j,k,ivar))/4

                   if (uav.gt.umaxUser(ivar).or.uav.lt.uminUser(ivar)) then
                      preg(npoly)=256+nreg+1
                   else
                      preg(npoly)=int((uav-uminUser(ivar))/uspanUser(ivar)*255)+nreg+1
                   endif

                else
                   preg(npoly)=icol
                endif

                plist(npoly,1,:)=pp(i,j,k,:)
                plist(npoly,2,:)=pp(i,j,k+1,:)
                plist(npoly,3,:)=pp(i+1,j,k+1,:)
                plist(npoly,4,:)=pp(i+1,j,k,:)
                 
                pxcent(npoly)=sum(plist(npoly,:,1))/4
             endif
          enddo
       enddo
    enddo

! ----------------------------------------------------------------------
! k-directed faces
! ----------------------------------------------------------------------

    do i=icutmin,icutmax-1
       do j=jcutmin,jcutmax-1
          do k=kcutmin,kcutmax

             if (k.eq.kcutmax) then
                ireg1=-1
             else
                ireg1=iregup(i,j,k)
             endif

             if (k.eq.kcutmin) then
                ireg2=-1
             else
                ireg2=iregup(i,j,k-1)
             endif

             ireg1fnd=0
             ireg2fnd=0

             do ireg=1,nreg
                if (ireg1.eq.ireg) then
                   if (showreg(ireg).eq.1) then
                      ireg1fnd=1
                      icol=ireg1
                   else if (showreg(ireg).eq.2) then
                      ireg1fnd=1
                      icol=-1
                   endif
                endif

                if (ireg2.eq.ireg) then
                   if (showreg(ireg).eq.1) then
                      ireg2fnd=1
                      icol=ireg2
                   else if (showreg(ireg).eq.2) then
                      ireg2fnd=1
                      icol=-1
                   endif
                endif
             enddo
             
             if (ireg1fnd+ireg2fnd.eq.1) then
                
                npoly=npoly+1
                
                if (icol.eq.-1) then
                   uav=(u(i,j,k,ivar)+u(i+1,j,k,ivar)+u(i+1,j+1,k,ivar)+u(i,j+1,k,ivar))/4

                   if (uav.gt.umaxUser(ivar).or.uav.lt.uminUser(ivar)) then
                      preg(npoly)=256+nreg+1
                   else
                      preg(npoly)=int((uav-uminUser(ivar))/uspanUser(ivar)*255)+nreg+1
                   endif
                else
                   preg(npoly)=icol
                endif

                plist(npoly,1,:)=pp(i,j,k,:)
                plist(npoly,2,:)=pp(i+1,j,k,:)
                plist(npoly,3,:)=pp(i+1,j+1,k,:)
                plist(npoly,4,:)=pp(i,j+1,k,:)
                
                pxcent(npoly)=sum(plist(npoly,:,1))/4
             endif
          enddo
       enddo
    enddo

!    write (message,*) 'npoly=',npoly
!    ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

!    write (message,*) 'rnode,pp=',rnode(1,1,1,1),pp(1,1,1,1)
!    ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)

  end subroutine getpolyzou

! ######################################################################

  integer(4) function debug(hwnd,iMsg,wParam,lParam)
    !DEC$ ATTRIBUTES STDCALL :: debug

! ----------------------------------------------------------------------
! Accepts WM_WDEBUG messages in which wParam=1 means error message
! ----------------------------------------------------------------------

    use user32
    use gdi32, SetTextColor => MSFWIN$SetTextColor
    use kernel32

    use const
    use global, only : message

    integer(HANDLE) hwnd,wParam,lParam,hdc
    integer(UINT) iMsg

    integer, save :: cyChar,cxChar
    type(T_TEXTMETRIC) tm
    type(T_RECT), save :: rect
    type(T_PAINTSTRUCT) :: ps

    integer ret,i

    integer, parameter :: nbuf=100

    character(100), save :: buffer(nbuf)

    debug=0

    if (iMsg.eq.WM_CREATE) then
       hdc=GetDC(hwnd)

       ret=SelectObject(hdc,GetStockObject(SYSTEM_FIXED_FONT))

       ret=GetTextMetrics(hdc,tm)
       cxChar=tm.tmAveCharWidth
       cyChar=tm.tmHeight

       ret=ReleaseDC(hwnd,hdc)

       buffer=" "

       return

    else if (iMsg.eq.WM_SIZE) then
       rect%top=0
       rect%left=0
       rect%right=loword(lParam)
       rect%bottom=hiword(lParam)

       ret=UpdateWindow(hwnd)
       return

    else if (iMsg.eq.WM_WDEBUG) then

       do i=1,len(message)
          if (message(i:i).eq.char(0)) message(i:i)=' '
       enddo
       
       ret=ScrollWindow(hwnd,0,-cyChar,rect,rect)

       hdc=GetDC(hwnd)

       if (wParam.eq.1) then
          ret=SetTextColor(hdc,rgb(255,0,0))
          message='ERROR: '//trim(message)
       else
          ret=SetTextColor(hdc,0)
       endif

       ret=SelectObject(hdc,GetStockObject(SYSTEM_FIXED_FONT))
       ret=TextOut(hdc,0,rect%bottom-cyChar,message,100)

       do i=1,nbuf-1  ! 1 is off the top of the window
          buffer(i)=buffer(i+1)
       enddo

       buffer(nbuf)=message

       ret=ReleaseDC(hwnd,hdc)
       ret=ValidateRect(hwnd,NULL)

       return

    else if (iMsg.eq.WM_PAINT) then
       ret=InvalidateRect(hwnd,NULL,FALSE)

       hdc=BeginPaint(hwnd,ps)

       ret=SelectObject(hdc,GetStockObject(SYSTEM_FIXED_FONT))

       do i=1,nbuf
          if (buffer(i)(1:5).eq.'ERROR') then
             ret=SetTextColor(hdc,rgb(255,0,0))
          else
             ret=SetTextColor(hdc,0)
          endif

          ret=TextOut(hdc,0,rect%bottom-cyChar*(nbuf-i+1),buffer(i),100)
       enddo

       ret=EndPaint(hwnd,ps)

       return

    else if (iMsg.eq.WM_DESTROY) then
       call PostQuitMessage(0)
       return
    else
       debug=DefWindowProc(hwnd,iMsg,wParam,lParam)
       return
    endif
  end function debug

! ######################################################################

  integer(4) function DlgProc(hwnd,message,wParam,lParam)
    !DEC$ ATTRIBUTES STDCALL :: DlgProc

! ----------------------------------------------------------------------
! Bounding box dialog
! ----------------------------------------------------------------------
    
    use ifwinty
    use user32
    use global, only : icutmax,icutmin,jcutmax,jcutmin,kcutmax,kcutmin,&
         slice,imax,jmax,kmax,btnstate
    
    integer(HANDLE) hwnd,wParam,lParam
    integer(UINT) message

    integer ret,i1,i2,j1,j2,k1,k2
    integer(LPBOOL) success
    logical iok,jok,kok
    character(100) string

    DlgProc=FALSE

    if (message.eq.WM_INITDIALOG) then
       write (string,'(a,i5,a,i5,a)') 'i=[',0,',',imax,']'C
       ret=SetDlgItemText(hwnd,209,string)

       write (string,'(a,i5,a,i5,a)') 'j=[',0,',',jmax,']'C
       ret=SetDlgItemText(hwnd,210,string)

       write (string,'(a,i5,a,i5,a)') 'k=[',0,',',kmax,']'C
       ret=SetDlgItemText(hwnd,211,string)

       ret=SetDlgItemInt(hwnd,203,icutmin,FALSE)
       ret=SetDlgItemInt(hwnd,204,icutmax,FALSE)
       
       ret=SetDlgItemInt(hwnd,205,jcutmin,FALSE)
       ret=SetDlgItemInt(hwnd,206,jcutmax,FALSE)
       
       ret=SetDlgItemInt(hwnd,207,kcutmin,FALSE)
       ret=SetDlgItemInt(hwnd,208,kcutmax,FALSE)

       DlgProc=TRUE
    else if (message.eq.WM_COMMAND) then
       if (loword(wParam).eq.201) then

          i1=GetDlgItemInt(hwnd,203,success,FALSE)
          i2=GetDlgItemInt(hwnd,204,success,FALSE)

          j1=GetDlgItemInt(hwnd,205,success,FALSE)
          j2=GetDlgItemInt(hwnd,206,success,FALSE)

          k1=GetDlgItemInt(hwnd,207,success,FALSE)
          k2=GetDlgItemInt(hwnd,208,success,FALSE)

          iok= i1.ge.0.and.i1.le.imax.and. &
               i2.ge.0.and.i2.le.imax.and.i1.lt.i2

          jok= j1.ge.0.and.j1.le.jmax.and. &
               j2.ge.0.and.j2.le.jmax.and.j1.lt.j2

          kok= k1.ge.0.and.k1.le.kmax.and. &
               k2.ge.0.and.k2.le.kmax.and.k1.lt.k2

          if (.not.(iok.and.jok.and.kok)) then
             ret=MessageBox(hwnd,&
                  'Selection out of range.\ni, j, k must be in range shown\nAnd imax>imin; jmax>jmin; kmax>kmin'C,NULL,MB_OK)
          else
             icutmin=i1; icutmax=i2
             jcutmin=j1; jcutmax=j2
             kcutmin=k1; kcutmax=k2

             btnstate(4:7)=[0,0,0,1]
             slice=0

             ret=EndDialog(hwnd,201)
          endif
          DlgProc=TRUE

       else if (loword(wParam).eq.202) then
          ret=EndDialog(hwnd,202)
          DlgProc=TRUE

       else if (loword(wParam).eq.212) then
          ret=SetDlgItemInt(hwnd,203,0,FALSE)
          ret=SetDlgItemInt(hwnd,204,imax,FALSE)

          ret=SetDlgItemInt(hwnd,205,0,FALSE)
          ret=SetDlgItemInt(hwnd,206,jmax,FALSE)

          ret=SetDlgItemInt(hwnd,207,0,FALSE)
          ret=SetDlgItemInt(hwnd,208,kmax,FALSE)
          
          DlgProc=TRUE
       else if (loword(wParam).eq.213) then
          i1 =GetDlgItemInt(hwnd,203,success,FALSE)
          ret=SetDlgItemInt(hwnd,204,i1+1,FALSE)

          j1 =GetDlgItemInt(hwnd,205,success,FALSE)
          ret=SetDlgItemInt(hwnd,206,j1+1,FALSE)

          k1 =GetDlgItemInt(hwnd,207,success,FALSE)
          ret=SetDlgItemInt(hwnd,208,k1+1,FALSE)

          DlgProc=TRUE

       endif
    endif
    
  end function DlgProc

! ######################################################################

  integer(4) function DlgDistort(hwnd,message,wParam,lParam)
    !DEC$ ATTRIBUTES STDCALL :: DlgDistort

! ----------------------------------------------------------------------
! Distort dialog
! ----------------------------------------------------------------------
    
    use ifwinty
    use user32
    use global, only : xscale,yscale,zscale
    
    integer(HANDLE) hwnd,wParam,lParam
    integer(UINT) message

    integer ret
    integer(LPBOOL) success
    character(100) string

    DlgDistort=FALSE

    if (message.eq.WM_INITDIALOG) then

       ret=SetDlgItemInt(hwnd,601,xscale,FALSE)
       ret=SetDlgItemInt(hwnd,602,yscale,FALSE)
       ret=SetDlgItemInt(hwnd,603,zscale,FALSE)

       DlgDistort=TRUE
    else if (message.eq.WM_COMMAND) then
       if (loword(wParam).eq.604) then         ! OK

          xscale=GetDlgItemInt(hwnd,601,success,FALSE)
          yscale=GetDlgItemInt(hwnd,602,success,FALSE)
          zscale=GetDlgItemInt(hwnd,603,success,FALSE)

          ret=EndDialog(hwnd,604)
          DlgDistort=TRUE

       else if (loword(wParam).eq.605) then    ! Cancel
          ret=EndDialog(hwnd,605)
          DlgDistort=TRUE

       else if (loword(wParam).eq.606) then    ! Cubic
          
          DlgDistort=TRUE
       endif
    endif
    
  end function DlgDistort

! ######################################################################

  integer(4) function DlgSetMax(hwnd,message,wParam,lParam)
    !DEC$ ATTRIBUTES STDCALL :: DlgSetMax

    use ifwinty
    use user32
    use global, only : umax,umaxUser,umin,uminUser,uspanUser,ivar

    integer(HANDLE) hwnd,wParam,lParam
    integer(UINT) message
    integer ret,ios,i
    character(100) string
    double precision val

    DlgSetMax=FALSE

    if (message.eq.WM_INITDIALOG) then
       DlgSetMax=TRUE
    else if (message.eq.WM_COMMAND) then
       if (loword(wParam).eq.703) then          ! Cancel
          ret=EndDialog(hwnd,703)
          DlgSetMax=TRUE
       else if (loword(wParam).eq.702) then     ! OK

          DlgSetMax=TRUE

          string=''
          ret=GetDlgItemText(hwnd,701,string,100)
          i=index(string,char(0))
          if (i.ne.0) string(i:i)=' '

          read (string,*,iostat=ios) val

          if (ios.ne.0) then
             call winerror('Incomprehensible number when setting upper limit: '//trim(string)) 
             return
          endif

          if (val.gt.umax(ivar).or.val.lt.umin(ivar)) then
             write (string,'(a,g13.5,a,g13.5,a)') 'Limit must be in range [',umin(ivar),',',umax(ivar),']'
             call winerror(string)
             return
          endif

          umaxUser(ivar)=val
          uspanUser(ivar)=umaxUser(ivar)-uminUser(ivar)
          ret=EndDialog(hwnd,702)

       else if (loword(wParam).eq.704) then            ! Max
          umaxUser(ivar)=umax(ivar)
          uspanUser(ivar)=umaxUser(ivar)-uminUser(ivar)
          ret=EndDialog(hwnd,704)
          
          DlgSetMax=TRUE
       endif
    endif

  end function DlgSetMax

! ######################################################################

  integer(4) function DlgSetMin(hwnd,message,wParam,lParam)
    !DEC$ ATTRIBUTES STDCALL :: DlgSetMin

    use ifwinty
    use user32
    use global, only : umax,umaxUser,umin,uminUser,uspanUser,ivar

    integer(HANDLE) hwnd,wParam,lParam
    integer(UINT) message
    integer ret,ios,i
    character(100) string
    double precision val

    DlgSetMin=FALSE

    if (message.eq.WM_INITDIALOG) then
       DlgSetMin=TRUE
    else if (message.eq.WM_COMMAND) then
       if (loword(wParam).eq.803) then          ! Cancel
          ret=EndDialog(hwnd,803)
          DlgSetMin=TRUE
       else if (loword(wParam).eq.802) then     ! OK

          DlgSetMin=TRUE

          string=''
          ret=GetDlgItemText(hwnd,801,string,100)
          i=index(string,char(0))
          if (i.ne.0) string(i:i)=' '

          read (string,*,iostat=ios) val

          if (ios.ne.0) then
             call winerror('Incomprehensible number when setting upper limit: '//trim(string)) 
             return
          endif

          if (val.gt.umax(ivar).or.val.lt.umin(ivar)) then
             write (string,'(a,g13.5,a,g13.5,a)') 'Limit must be in range [',umin(ivar),',',umax(ivar),']'
             call winerror(string)
             return
          endif

          uminUser(ivar)=val
          uspanUser(ivar)=umaxUser(ivar)-uminUser(ivar)
          ret=EndDialog(hwnd,802)

       else if (loword(wParam).eq.804) then            ! Max
          uminUser(ivar)=umin(ivar)
          uspanUser(ivar)=umaxUser(ivar)-uminUser(ivar)
          ret=EndDialog(hwnd,804)
          
          DlgSetMin=TRUE
       endif
    endif

  end function DlgSetMin

! ######################################################################

  integer(4) function DlgCol(hwnd,message,wParam,lParam)
    !DEC$ ATTRIBUTES STDCALL :: DlgCol

! ----------------------------------------------------------------------
! Colour dialog
! ----------------------------------------------------------------------
    
    use ifwinty
    use user32
    
    integer(HANDLE) hwnd,wParam,lParam
    integer(UINT) message

    integer ret,red,green,blue
    integer(LPBOOL) success
    
    DlgCol=FALSE

    if (message.eq.WM_INITDIALOG) then
       DlgCol=TRUE
    else if (message.eq.WM_COMMAND) then
       if (loword(wParam).eq.304) then

          red=GetDlgItemInt(hwnd,301,success,FALSE)
          green=GetDlgItemInt(hwnd,302,success,FALSE)
          blue=GetDlgItemInt(hwnd,303,success,FALSE)

          if (red.gt.255.or.green.gt.255.or.blue.gt.255) then
             ret=MessageBox(hwnd,'Colours must be in range 0..255'C,NULL,MB_OK)
             DlgCol=TRUE
             return
          endif

          ret=EndDialog(hwnd,rgb(red,green,blue))
          DlgCol=TRUE
       else if (loword(wParam).eq.305) then
          ret=EndDialog(hwnd,-1)
          DlgCol=TRUE
       endif
    endif
    
  end function DlgCol

! ######################################################################

  integer(4) function DlgView(hwnd,message,wParam,lParam)
    !DEC$ ATTRIBUTES STDCALL :: DlgView

! ----------------------------------------------------------------------
! View MIN file. MinName is global and has been found during FILE > open min
! note that pszFileText must be passed as a pointer to ReadFile
! (see interface in kernel32)
! on the other hand SetDlgItemText expect a character(*)
! pzFileText is an array of character(1) since we cannot easily set
! character(len=whatever get file size returns)
! ----------------------------------------------------------------------

    use kernel32 ! CreateFile,GetFileSize,ReadFile,CloseHandle
    use user32   ! SetDlgItemText,EndDialog
    use global, only : MinName

    integer(HANDLE) hwnd,wParam,lParam
    integer(UINT) message
    
    integer(HANDLE) hFile
    integer(DWORD) dwFileSize
    character(1), allocatable :: pszFileText(:)
    integer(DWORD) dwRead
    integer ret

    DlgView=FALSE

    if (message.eq.WM_INITDIALOG) then

       ret=SetDlgItemText(hwnd,403,MinName)

       hFile = CreateFile(MinName, GENERIC_READ, FILE_SHARE_READ, NULL, &
            OPEN_EXISTING, 0, NULL)

       if (hFile .ne. INVALID_HANDLE_VALUE) then
      
          dwFileSize = GetFileSize(hFile, NULL)

          if (dwFileSize.gt.0) then

             allocate(pszFileText(dwFileSize+1))

             if (ReadFile(hFile, loc(pszFileText), dwFileSize, loc(dwRead), NULL).eq.TRUE) then
                pszFileText(dwFileSize+1) = char(0) ! Add null terminator
                ret=SetDlgItemText(hwnd,401,pszFileText(1))
             endif

             deallocate(pszFileText)
          endif
          ret=CloseHandle(hFile)
       endif
       DlgView=TRUE
    else if (message.eq.WM_COMMAND) then
       if (loword(wParam).eq.402) then
          ret=EndDialog(hwnd,-1)
          DlgView=TRUE
       endif
    endif

  end function DlgView

! ######################################################################

  integer(4) function DlgHelp(hwnd,message,wParam,lParam)
    !DEC$ ATTRIBUTES STDCALL :: DlgHelp

! ----------------------------------------------------------------------
! View MIN file. MinName is global and has been found during FILE > open min
! note that pszFileText must be passed as a pointer to ReadFile
! (see interface in kernel32)
! on the other hand SetDlgItemText expect a character(*)
! pzFileText is an array of character(1) since we cannot easily set
! character(len=whatever get file size returns)
! ----------------------------------------------------------------------

    use user32   ! SetDlgItemText,EndDialog
    use global, only : hDlgHelp

    integer(HANDLE) hwnd,wParam,lParam
    integer(UINT) message
    
    integer ret

    character(2), parameter :: NL=char(13)//char(10)

    character(*), parameter :: helptxt=&
         'General:'//NL// &
         '* Use the mouse to drag and rotate the view'//NL// &
         '* Hold <shift> and drag to shift the view'//NL// &
         '* Use A/Z to zoom in and out'//NL// &
         '* Use N/M to rotate view in the plane'//NL// &
         '* You can also rotate the view using the arrow keys instead of the mouse'//NL// &
         ''//NL// &
         'Colour Key (on right of screen):'//NL// &
         '* Click with left mouse button to toggle region on or off'//NL// &
         ' (you can also show DATA if you are viewing a .ZOU file. The sequence is on -> DATA -> off)'//NL// &
         '* Click the Key with the right mouse button to change the colour scheme'//NL// &
         ''//NL// &
         'Slices:'//NL// &
         '* Select x/y/z-slice using buttons at the bottom of the screen'//NL// &
         '* Only a slice through the system is shown letting you get a look at what is going on inside'//NL// &
         '* To change which slice is being viewed, use the up and down arrow keys'//NL// &
         '* You can continue to rotate/shift in the normal way and can still toggle regions on or off'//NL// &
         '* Combining these features gives a variety of views'//NL// &
         '* Click "all" to switch off slices and show the full view again'//NL// &
         ''//NL// &
         'Bounding Box:'//NL// &
         '* Allows you to view a cuboid section from the total mesh'//NL// &
         '* Use View > Bounding Box and set imin, imax etc to numbers in the range given'//NL// &
         '* The "all" button resets the selection to the whole system'//NL// &
         '* The "one" button selects just one element, given by [imin,jmin,kmin]'//NL// &
         ''//NL// &
         'x/y/z buttons (bottom): sets the view along the indicated axis'//NL// &
         ''//NL// &
         'Mesh viewing'//NL// &
         '* Load an existing mesh file from the File > Open MTF option'//NL// &
         ''//NL// &
         'Mesh generation'//NL// &
         '* Open a .MIN file (contains geometry description) using FILE> Open MIN'//NL// &
         '* Then use File > Process mesh'//NL// &
         '* The resulting mesh is displayed on screen and automatically saved as <filename>.mtf'//NL// &
         ''//NL// &
         'Viewing simulation output'//NL// &
         '* Zmesh can be used to view the "raw" output of Zinc simulations'//NL// &
         '* For more advanced views and post processing, line/surface graphs etc. use zpp instead'//NL// &
         '* Select a .zou file using File > Open .ZOU'//NL// &
         '* The colour key on the right now has three states: on -> DATA -> mesh'//NL// &
         '* The DATA setting shows the value of one of the variables, colour coded with a legend on the left'//NL// &
         '* Change the limits of the legend by clicking "set max" and "set min"'//NL// &
         '* You can do slices as usual, both elements and nodes can show colour-coded data'//NL// &
         '* To change which variable you want to view, click the < and > buttons at the top-left of the screen'//NL// &
         '* Use CTRL+Click to print the solution value of one element'//NL// &
         '* Note that zmesh will also open <filename>.mtf and <filename>.zin.'//NL// &
         ''//NL// &
         'For further help, see the Zinc User Manual and Tutorial Manual'C

    DlgHelp=FALSE

    if (message.eq.WM_INITDIALOG) then

       ret=SetDlgItemText(hwnd,501,helptxt)
       DlgHelp=TRUE
    else if (message.eq.WM_COMMAND) then
       if (loword(wParam).eq.502) then
!          ret=EndDialog(hwnd,-1)
          ret=DestroyWindow(hwnd)
          DlgHelp=TRUE
          hDlgHelp=0
       endif
    else if (message.eq.WM_CLOSE) then
!       ret=EndDialog(hwnd,-1)
       ret=DestroyWindow(hwnd)
       DlgHelp=TRUE
       hDlgHelp=0
    endif

  end function DlgHelp

! ######################################################################

  subroutine freekey

    use global
    use gdi32

    integer i,ret

    do i=1,nreg
       ret=DeleteObject(cols(i))
    enddo

    do i=1,nregnode
       ret=DeleteObject(nodecols(i))
    enddo

    deallocate (cols,nodecols)
    deallocate (showreg,nodeshowreg)
  end subroutine freekey

! ######################################################################

  subroutine allockey

    use global
    use gdi32

    integer ind,intens,ir,ig,ib,colour,allocerr

    allocate (cols(nreg),nodecols(nregnode),stat=allocerr)
    if (allocerr.ne.0) then
       call winerror('Failed to allocate colour arrays')
       return
    endif

    allocate (showreg(nreg),nodeshowreg(nregnode),stat=allocerr)
    if (allocerr.ne.0) then
       call winerror('Failed to allocate node display arrays')
       return
    endif

    ind=0
    outer: do
       do intens=1,2
          do ir=0,1
             do ig=0,1
                do ib=0,1
                   if (ir+ig+ib.ne.0) then
                      ind=ind+1
                      colour=rgb(ir*255/intens,ig*255/intens,ib*255/intens)
                      cols(ind)=CreateSolidBrush(colour)
                      
                      if (ind.eq.nreg) exit outer
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo outer
    
    ind=0
    outernode: do
       do intens=1,2
          do ir=0,1
             do ig=0,1
                do ib=0,1
                   if (ir+ig+ib.ne.0) then
                      ind=ind+1
                      colour=rgb(ir*255/intens,ig*255/intens,ib*255/intens)
                      nodecols(ind)=CreateSolidBrush(colour)
                      
                      if (ind.eq.nregnode) exit outernode
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo outernode
  end subroutine allockey

! ######################################################################

  subroutine winerror(string)

    use global
    use user32
    use kernel32
    use const

    character(*) string
    integer ret

!    message='ERROR: '//trim(string)
    message=string
    ret=SendMessage(hwndDebug,WM_WDEBUG,1,0)
  end subroutine winerror

  subroutine winmess(string)

    use global
    use user32
    use kernel32
    use const

    character(*) string
    integer ret

    message=string
    ret=SendMessage(hwndDebug,WM_WDEBUG,0,0)
  end subroutine winmess

! ######################################################################

  subroutine initview(hwnd)

! ----------------------------------------------------------------------
! set pp. Find number of regions and set initial showreg
! ----------------------------------------------------------------------

    use global
    use user32
    use kernel32
    use comdlg32

    integer ic,ret
    integer(HANDLE) hwnd,hMenu

    !    pp=rnode     ! pp will get rotated
    pp(:,:,:,1)=rnode(:,:,:,1)*xscale
    pp(:,:,:,2)=rnode(:,:,:,2)*yscale
    pp(:,:,:,3)=rnode(:,:,:,3)*zscale

    if (ZOU) then
       showreg=2
    else
       showreg=1
    endif

    nodeshowreg=0
    ivar=1       ! default for ZOU
    iview=1      ! default for EXTRA view

! ----------------------------------------------------------------------
! find the bounding box in redge(3,2). Find number of regions
! for elements, faces and nodes, nreg(1-3). 
! note: imode =1,2,3 for elements, faces and node display
! rescale so that object is shown nicely
! set slice to zero (full view). Set bounding box to full
! ----------------------------------------------------------------------

    do ic=1,3
       redge(ic,1)=minval(rnode(:,:,:,ic))
       redge(ic,2)=maxval(rnode(:,:,:,ic))
    enddo
    
    call rescale(width,height)
    
    slice=0
    
    icutmax=imax
    icutmin=0
             
    jcutmax=jmax
    jcutmin=0
    
    kcutmax=kmax
    kcutmin=0
    
! ----------------------------------------------------------------------
! Set bb and ff for bounding box display
! InvalidateRect to draw new picture. Set loaded.
! Activate "View" menu items
! ----------------------------------------------------------------------
    
    call setbb
    
    hMenu=GetMenu(hwnd)
    ret=EnableMenuItem(hMenu,102,MF_ENABLED)
    ret=EnableMenuItem(hMenu,103,MF_ENABLED)
    ret=EnableMenuItem(hMenu,108,MF_ENABLED)

    ScrollEl=1
    ScrollNd=1

    loaded=.true.
    
    ret=InvalidateRect(hwnd,NULL,FALSE)
  end subroutine initview

! ######################################################################

  function toupper(string)

    character(*) string
    character(len(string)) toupper
    character(26), parameter :: lower='abcdefghijklmnopqrstuvwxyz'
    character(26), parameter :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    character(1) c
    integer i,j

    toupper=string

    do i=1,len_trim(string)
       c=string(i:i)
       j=index(lower,c)
       if (j.ne.0) toupper(i:i)=upper(j:j)
    enddo
  end function toupper

! ######################################################################

  function pointInPolygon(vertices,x,y)

    use gdi32, Polygon => MSFWIN$Polygon

    logical pointInPolygon
    type(T_POINT) vertices(4)
    integer x,y
    integer intersectionCount,x0,y0,x1,y1,i

    intersectionCount = 0
  
    x0 = vertices(4)%x - x
    y0 = vertices(4)%y - y
    
    do i=1,4
       x1 = vertices(i)%x - x
       y1 = vertices(i)%y - y
 
       if (y0 > 0 .and. y1 <= 0 .and. x1 * y0 > y1 * x0) then
          intersectionCount = intersectionCount + 1
       endif
    
       if (y1 > 0 .and. y0 <= 0 .and. x0 * y1 > y0 * x1) then
          intersectionCount = intersectionCount + 1
       endif
 
       x0 = x1
       y0 = y1	
    enddo
 
    pointInPolygon=mod(intersectionCount,2) == 1
  end function pointInPolygon

end module win
