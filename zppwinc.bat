rem windres -I c:\progra~1\gfortran\include zppwin.rc zppwin_res.o

windres zppwin.rc zppwin_res.o
gcc -o zppwin.exe zppwin.c zppwin_res.o -lcomdlg32 -lgdi32
