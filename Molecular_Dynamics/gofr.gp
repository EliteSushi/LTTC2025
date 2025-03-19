set terminal png size 800,450 enhanced font "TeX Gyre Termes,20"
set output 'LaTEX/Gr.png'
set border linewidth 3.5

set xlabel 'r (Å)'
set ylabel 'g(r)'
set xrange [0:12.5]


p 'gofrArPBC.dat' u 1:2 w lp lw 3.5 t 'g(r) Ar-Ar', \