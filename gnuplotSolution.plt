set terminal pdf
set output 'img/solution.pdf'
set title 'Solution'
set xlabel 'x'
set ylabel 'y'
splot "data/solution.txt" with pm3d
#pause -1
