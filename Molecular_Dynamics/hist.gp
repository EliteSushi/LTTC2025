datafile = 'gauss_bin' 
binwidth = 0.2

set terminal png size 1200,850 enhanced font "TeX Gyre Termes,20"
set output 'report/Gauss.png'
set border linewidth 3.5

bin(x) = binwidth * floor(x / binwidth)
set style data histogram
set style fill pattern 5 border lt -1
set xlabel 'x'
set ylabel 'Frequency'
set xrange [-4:4]

# Gaussian function
mean = 0
stddev = 1
gauss(x) = (1 / (stddev * sqrt(2 * pi))) * exp(-(x - mean)**2 / (2 * stddev**2))

stats datafile nooutput 

# Plot
plot datafile using 1:2 with boxes title '10^6 gaussian\_distr()', \
     gauss(x) / gauss(0) * STATS_max_y title 'Gaussian (μ=0, σ=1)' lw 8
