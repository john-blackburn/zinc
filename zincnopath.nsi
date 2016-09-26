!include "EnvVarUpdate.nsh"

name "Zinc 3.31"
Icon "zmesh\z.ico"

# define name of installer
outFile "zinc331nopathsetup.exe"
 
# define installation directory
installDir "c:\program files\Zinc"

LicenseData "zinclicence.rtf"

Page license
Page directory
Page instfiles
#UninstPage uninstConfirm
 
# start default section
section
 
    # set the installation directory as the destination for the following actions
    setOutPath $INSTDIR

    file zinc.exe
    file zpp.exe
    file zincwin.exe
    file zppwin.exe
    file zmesh\zmesh.exe
    file zmesh\zmeshwin.exe
    file zincman.pdf
    file zinctheo.pdf
    file tutorial.pdf
    file nltemplate.f90
    file nltemplate.c

    file zincrun.bat
    file zpprun.bat

    file c:\mingw\bin\libgcc_s_dw2-1.dll

    SetOutPath "$INSTDIR\examples"
    file /r /x CVS examples\*.*

    SetOutPath "$INSTDIR\gnuplot"
    file /r c:\gnuplot42\*.*

    SetOutPath "$INSTDIR\mesh_examples"
    file /r zmesh\*.min

    # create the uninstaller
    writeUninstaller "$INSTDIR\uninstall.exe"
 
    # create a shortcut named "new shortcut" in the start menu programs directory
    # point the new shortcut at the program uninstaller

    CreateDirectory "$SMPROGRAMS\Zinc 3.31"

    createShortCut "$SMPROGRAMS\Zinc 3.31\gnuplot.lnk" "$INSTDIR\gnuplot\bin\gnuplot.exe"
    createShortCut "$SMPROGRAMS\Zinc 3.31\zmesh.lnk" "$INSTDIR\zmeshwin.exe"
    createShortCut "$SMPROGRAMS\Zinc 3.31\zinc.lnk" "$INSTDIR\zincwin.exe"
    createShortCut "$SMPROGRAMS\Zinc 3.31\zpp.lnk" "$INSTDIR\zppwin.exe"
    createShortCut "$SMPROGRAMS\Zinc 3.31\Manual.lnk" "$INSTDIR\zincman.pdf"
    createShortCut "$SMPROGRAMS\Zinc 3.31\Theory.lnk" "$INSTDIR\zinctheo.pdf"
    createShortCut "$SMPROGRAMS\Zinc 3.31\Tutorial.lnk" "$INSTDIR\Tutorial.pdf"

    createShortCut "$SMPROGRAMS\Zinc 3.31\Uninstall.lnk" "$INSTDIR\uninstall.exe"

#    ${EnvVarUpdate} $0 "PATH" "A" "HKLM" "$INSTDIR"
#    ${EnvVarUpdate} $0 "PATH" "A" "HKLM" "$INSTDIR\gnuplot\bin"

#    WriteRegExpandStr ${hklm_all_users} ZINCDIR $INSTDIR
#    SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000

sectionEnd
 
# uninstaller section start
section "uninstall"
 
    # first, delete the uninstaller
    delete "$INSTDIR\uninstall.exe"

    delete "$INSTDIR\zinc.exe"
    delete "$INSTDIR\zpp.exe"
    delete "$INSTDIR\zincwin.exe"
    delete "$INSTDIR\zppwin.exe"
    delete "$INSTDIR\zmesh.exe"
    delete "$INSTDIR\zmeshwin.exe"
    delete "$INSTDIR\zincman.pdf"
    delete "$INSTDIR\zinctheo.pdf"
    delete "$INSTDIR\tutorial.pdf"
    delete "$INSTDIR\nltemplate.f90"
    delete "$INSTDIR\nltemplate.c"

    delete "$INSTDIR\zincrun.bat"
    delete "$INSTDIR\zpprun.bat"

    delete "$INSTDIR\libgcc_s_dw2-1.dll"

    RMDir /r "$INSTDIR\examples"
    RMDir /r "$INSTDIR\gnuplot"
    RMDir /r "$INSTDIR\mesh_examples"

    RMDir $INSTDIR

    # second, remove the link from the start menu
    delete "$SMPROGRAMS\Zinc 3.31\Uninstall.lnk"

    delete "$SMPROGRAMS\Zinc 3.31\gnuplot.lnk"
    delete "$SMPROGRAMS\Zinc 3.31\zmesh.lnk"
    delete "$SMPROGRAMS\Zinc 3.31\zinc.lnk"
    delete "$SMPROGRAMS\Zinc 3.31\zpp.lnk"
    delete "$SMPROGRAMS\Zinc 3.31\Manual.lnk"
    delete "$SMPROGRAMS\Zinc 3.31\Tutorial.lnk"
    delete "$SMPROGRAMS\Zinc 3.31\Theory.lnk"

    RMDir "$SMPROGRAMS\Zinc 3.31"

    # third, remove path extras and the ZINCDIR environment variable
 
#    ${un.EnvVarUpdate} $0 "PATH" "R" "HKLM" "$INSTDIR"
#    ${un.EnvVarUpdate} $0 "PATH" "R" "HKLM" "$INSTDIR\gnuplot\bin"

#    DeleteRegValue ${hklm_all_users} ZINCDIR
#    SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000

# uninstaller section end
sectionEnd
