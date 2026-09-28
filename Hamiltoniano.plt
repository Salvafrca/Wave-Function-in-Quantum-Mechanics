set term gif animate delay 10 loop 0
set output 'EnergiaMedia.gif'

# Nuevos ejes ajustados a mediop.dat
set yrange[80:90]
set xrange[0:0.2]

# Mismo bucle con estela, pero leyendo el nuevo archivo
do for [a=0:100] {
    plot "Hamiltoniano.dat" index 0:a using 1:2 with lines linewidth 2 linecolor "blue" , \
         "Hamiltoniano.dat" index a using 1:2 with points pointtype 7 pointsize 2 linecolor "red" title "<Energia Media>"
}