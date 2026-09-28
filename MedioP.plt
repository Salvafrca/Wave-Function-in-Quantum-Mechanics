set term gif animate delay 10 loop 0
set output 'MedioP.gif'

# Nuevos ejes ajustados a mediop.dat
set xrange [0:0.2]
set yrange [-50:50] 

# Mismo bucle con estela, pero leyendo el nuevo archivo
do for [a=0:100] {
    plot "mediop.dat" index 0:a using 1:2 with lines linewidth 2 linecolor "blue" , \
         "mediop.dat" index a using 1:2 with points pointtype 7 pointsize 2 linecolor "red" title "<p>"
}