datafile = 'boltz_bin' 
binwidth = 0.2

set terminal png size 1200,850 enhanced font "TeX Gyre Termes,20"
set output 'report/Boltz.png'
set border linewidth 3.5

bin(x) = binwidth * floor(x / binwidth)
set style data histogram
set style fill pattern 5 border lt -1
set xlabel 'x'
set ylabel 'Frequency'
set xrange [0:6.2]

# boltzian function
mean = 0
stddev = 1
kbT_m = 1.8
boltz(x) = 4 * pi * sqrt(kbT_m / (2 * pi))**(3/2) * x**2 * exp(-x**2 / (2 * kbT_m))

stats datafile nooutput 

# Plot
plot datafile using 1:2 with boxes title '10^6 initial\_velocities()', \
     boltz(x)*0.01 title 'f(v)' lw 8
#/ boltz(0) * STATS_max_y 
