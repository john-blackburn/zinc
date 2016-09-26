set term emf

set output 'uzcomp.emf'
set xlabel 'z (cm)'
set ylabel 'uz (pm)'
plot 'piezo09.out' u ($1*100):($5*1e12) t 'Zinc','uz.txt' i 0 u 1:($2*1e12) t 'Comsol' w l

set output 'uxcomp.emf'
set xlabel 'z (cm)'
set ylabel 'ux (pm)'
plot 'piezo08.out' u ($1*100):($5*1e12) t 'Zinc','ux.txt' i 0 u 1:($2*1e12) t 'Comsol' w l
