set term emf mono dashed

set output 'vscan.emf'
set xlabel 'Scan distance (m)'
set ylabel 'V (Volts)'
unset key

plot 'dielectric01.out' u 1:5 w l

set output 'Ezscan.emf'
set ylabel 'Ez (V/m)'
set key

plot 'dielectric03.out' u 1:5 t 'x=0.1 m' w l,\
'dielectric04.out' u 1:5 t 'x=0.0 m' w l

set output 'Dzmap.emf'
set title 'Dz (C/m^2)'
set xlabel 'x (m)'
set ylabel 'z (m)'
set size square
set pm3d map
splot 'dielectric07.out' u 1:3:4
