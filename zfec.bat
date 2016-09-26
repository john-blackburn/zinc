windres zincfe.rc zincfe_res.o
gcc  -o zincwin.exe zincfe.c zincfe_res.o -lcomdlg32 -lgdi32
