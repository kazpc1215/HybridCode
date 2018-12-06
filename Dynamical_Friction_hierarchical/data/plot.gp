unset multiplot
reset

set term aqua dashed font "Times-Roman,30" enhanced

PAUSE = -1


set bmargin 2
set lmargin 7
# set lmargin 6
set xtics offset 0,0.3
set ytics offset 0.5,0


set log

set format "10^{%L}"
# unset key

set grid

set bar 0
LW = 2
PS = 0.7
# DT = 2

###################################


set xl "time [yr]" offset 0,1
set yl "relative energy error" offset 2,0
set xr [1E-1:1E3+0.01]


set key left top box width -6.5 spacing 1.0 font "Times-Roman,20"
p "test/rand01/Energy.dat" u 1:(($3) > 0 ? ($3) : NaN) w lp lt 1 lw LW ps PS dt 1 t "w/ temporary file, > 0",\
"" u 1:(($3) < 0 ? -($3) : NaN) w lp lt 1 lw LW ps PS dt 2 t "< 0",\
"test_notemp/rand01/Energy.dat" u 1:(($3) > 0 ? ($3) : NaN) w lp lt 2 lw LW ps PS dt 1 t "w/o temporary file, > 0",\
"" u 1:(($3) < 0 ? -($3) : NaN) w lp lw LW ps PS lt 2 dt 2 t "< 0"

pause -1


set key left top box width -9 spacing 1.0 font "Times-Roman,20"
p "test/rand01/Energy.dat" u 1:(($3) > 0 ? ($3) : NaN) w lp lt 1 lw LW ps PS dt 1 t "relative energy error, > 0",\
"" u 1:(($3) < 0 ? -($3) : NaN) w lp lt 1 lw LW ps PS dt 2 t "< 0",\
"test/rand01/Energy.dat" u 1:(($4) > 0 ? ($4) : NaN) w lp lt 2 lw LW ps PS dt 1 t "collisional correction, > 0",\
"" u 1:(($4) < 0 ? -($4) : NaN) w lp lw LW ps PS lt 2 dt 2 t "< 0"


pause -1




set xl "step" offset 0,1
set xr [1E1:1E7]

set key left top box width -6.5 spacing 1.0 font "Times-Roman,20"
p "test/rand01/Energy.dat" u 5:(($3) > 0 ? ($3) : NaN) w lp lt 1 lw LW ps PS dt 1 t "w/ temporary file, > 0",\
"" u 5:(($3) < 0 ? -($3) : NaN) w lp lt 1 lw LW ps PS dt 2 t "< 0",\
"test_notemp/rand01/Energy.dat" u 5:(($3) > 0 ? ($3) : NaN) w lp lt 2 lw LW ps PS dt 1 t "w/o temporary file, > 0",\
"" u 5:(($3) < 0 ? -($3) : NaN) w lp lw LW ps PS lt 2 dt 2 t "< 0"

pause -1


set key left top box width -9 spacing 1.0 font "Times-Roman,20"
p "test/rand01/Energy.dat" u 5:(($3) > 0 ? ($3) : NaN) w lp lt 1 lw LW ps PS dt 1 t "relative energy error, > 0",\
"" u 5:(($3) < 0 ? -($3) : NaN) w lp lt 1 lw LW ps PS dt 2 t "< 0",\
"test/rand01/Energy.dat" u 5:(($4) > 0 ? ($4) : NaN) w lp lt 2 lw LW ps PS dt 1 t "collisional correction, > 0",\
"" u 5:(($4) < 0 ? -($4) : NaN) w lp lw LW ps PS lt 2 dt 2 t "< 0"