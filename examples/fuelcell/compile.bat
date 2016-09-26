gfortran -c -fno-underscoring -fbounds-check cathodetoken.f90
gfortran -s -shared -mrtd  -o cathodetoken.dll cathodetoken.o
