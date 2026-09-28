set term gif animate

set output 'Prob.gif'
set yrange [0:10]



do for [a=0:100]{plot "fi_prob.dat" i a u 1:2}
