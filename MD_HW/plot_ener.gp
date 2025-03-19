set terminal png size 800,450 enhanced font "TeX Gyre Termes,20"
set output 'LaTEX/PBC_E.png'
set border linewidth 3.5

set xlabel 'Time (fs)'
set ylabel 'Energy (au)'

p 'ENERGIES' u 1:3 w l t 'Potenital E' lw 3.5, \
  'ENERGIES' u 1:4 w l t 'Kinetic E' lw 3.5, \
  'ENERGIES' u 1:($3+$4) w l t 'Total E' lw 3.5