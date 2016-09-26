module const

  use ifwinty
  save

  integer, parameter :: WM_WDEBUG=WM_USER    ! WM_USER .. 0x7FFF
  real, parameter :: pi=3.141592654

  integer, parameter :: keywidth=65,keyheight=25,keyxoff=20,keyyoff=20,keypitch=40

  real, parameter :: ascale=0.01,sscale=1.0

  integer, parameter :: btnxoff=10,btnyoff=10,btnwidth=70,btnheight=20,btnpitch=90

  integer, parameter :: legendxoff=20,legendyoff=70,legendwidth=20
  integer, parameter :: xvar=200,yvar=8,varsize=20,varsep=5,vartxt=4

  integer, parameter :: nbtn=7
  character(20), parameter :: btntext(nbtn)=&
       ['x','y','z','xslice','yslice','zslice','all']

  integer, parameter :: stwidth=700, stheight=700, stconwidth=300

  integer, parameter :: maxbutton=10

end module const
