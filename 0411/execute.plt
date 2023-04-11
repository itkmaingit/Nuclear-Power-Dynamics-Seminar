# グラフの設定
set title "Comparison of Exact and Numerical Graphs"
set xlabel "t"
set ylabel "X"
set key left top
set yrange[:5]
set xrange[:20]

# Exactグラフの描画
plot "leapflog.dat" using 1:3 with lines lw 2 lt rgb "blue" title "Exact"

# Numericalグラフの描画
replot "leapflog.dat" using 1:4 with lines lw 2 lt rgb "red" title "Numerical"

set term pngcairo
set output "leapflog.png"
replot
set output