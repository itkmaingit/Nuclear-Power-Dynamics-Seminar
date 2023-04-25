set encoding utf8
set size square
set key box font ",24"
set xran [0:10.000000]
set yran [0:10.000000]
max_iter = 400
progress_width = 50
set terminal gif animate delay 5 size 900,900 optimize
set output "cell_10.000000_temperature_1.300000.gif"
do for [i=0:max_iter] {
progress_percent = (100*i) / (max_iter-1)
progress_bar = "["
progress_num = int(progress_percent / (100.0 / progress_width))
do for [j=0:progress_width-1] {
    if (j <= progress_num) {
        progress_bar = progress_bar . "="
    } else {
        progress_bar = progress_bar . " "
    }
}
progress_bar = progress_bar . "]"
file_num = sprintf("%08d", i*100)
file_name = "lj" . file_num . ".dat"
title_str = file_num
print sprintf("Progress: %s %5.2f%%", progress_bar, progress_percent)
plot file_name using 1:2 title title_str with points pt 6 ps 5 lw 1
}
