gfortran -c -O2 -Wall -fno-underscoring fullfuelmod.f90
gfortran -c -O2 -Wall -fno-underscoring fullfuel.f90

gfortran -s -shared -mrtd -static-libgfortran -o fullfuel.dll fullfuelmod.o fullfuel.o
