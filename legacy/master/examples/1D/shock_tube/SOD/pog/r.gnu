set style func linespoints
set style line 1 lt 1 lc rgb "red" lw 1 ps 1 pt 6
set style line 2 lt 1 lc rgb "blue" lw 1 ps 1 pt 6
plot 'output-ref/sod.g01.b001.gnu.dat' u 1:4 ls 2 w lp notitle, 'output-ref/sod.g01.b002.gnu.dat' u 1:4 ls 2 w lp title "OFF Reference",\
         'output/sod.g01.b001.gnu.dat' u 1:4 ls 1 w lp notitle,     'output/sod.g01.b002.gnu.dat' u 1:4 ls 1 w lp title "OFF output"
pause -1
