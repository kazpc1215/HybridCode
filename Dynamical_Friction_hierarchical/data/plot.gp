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
# PS = 0
# DT = 2

###################################

set xl "time [yr]" offset 0,1

OHTSUKI_NOFRAG_PLANET = "Meach3E-9_Mtot3E-7_Mmax5E-15_t1E9_dtlog_ecc1E-2_nofrag_dt/Planet.dat"
OHTSUKI_NOFRAG_PLANETESIMAL = "Meach3E-9_Mtot3E-7_Mmax5E-15_t1E9_dtlog_ecc1E-2_nofrag_dt/Planetesimal.dat"
RUN1 = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/RMS_randall.dat"

set yl "ecc" offset 2,0
set xr [1:1.E8+1]
 set yr [1E-3:1]
# set yr [0.003:1]
# set yr [0.005:0.5]


### nofrag_ecc  ###
set multiplot
set key left Left bottom box width -10 spacing 1.0 reverse font "Times-Roman,20"
plot OHTSUKI_NOFRAG_PLANET u 1:2 w l lt -1 lw 2 dt 1 t "Planet, semi-analytic, no frag.",\
RUN1 u 1:4:5 w yerrorlines lw 2 dt 1 lt 1 ps 0.5 t "Planet, r.m.s., no frag."
set key left Left top box width -12 spacing 1.0 reverse font "Times-Roman,20"
plot OHTSUKI_NOFRAG_PLANETESIMAL u 1:2 w l lt -1 lw 2 dt 2 t "Planetesimal, semi-analytic, no frag.",\
RUN1 u 1:6:7 w yerrorlines lw 2 dt 2 lt 1 ps 0.5 t "Planetesimal, r.m.s., no frag."
unset multiplot

pause -1


set yl "inc [rad]" offset 2,0
set xr [1:1.E8+1]
 set yr [1E-3:1]
# set yr [0.003:1]
# set yr [0.005:0.5]

### nofrag_inc  ###
set multiplot
set key left Left bottom box width -10 spacing 1.0 reverse font "Times-Roman,20"
plot OHTSUKI_NOFRAG_PLANET u 1:3 w l lt -1 lw 2 dt 1 t "Planet, semi-analytic, no frag.",\
RUN1 u 1:10:11 w yerrorlines lw 2 dt 1 lt 1 ps 0.5 t "Planet, r.m.s., no frag."
set key left Left top box width -12 spacing 1.0 reverse font "Times-Roman,20"
plot OHTSUKI_NOFRAG_PLANETESIMAL u 1:3 w l lt -1 lw 2 dt 2 t "Planetesimal, semi-analytic, no frag.",\
RUN1 u 1:12:13 w yerrorlines lw 2 dt 2 lt 1 ps 0.5 t "Planetesimal, r.m.s., no frag."
unset multiplot

pause -1


#########################################


set auto

set xl "time [yr]" offset 0,1
set yl "relative error of energy" offset 2,0
set xr [1E-1:1E8]



set key right bottom box width 0 spacing 1.0 font "Times-Roman,20"
p "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/rand01/Energy.dat" u 1:(($3) > 0 ? ($3) : NaN) w l lt 1 lw LW dt 1 t "> 0",\
"" u 1:(($3) < 0 ? -($3) : NaN) w l lt 1 lw LW dt 2 t "< 0",\

pause -1


set key right bottom box width -9.5 spacing 1.0 font "Times-Roman,20"
p "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/rand01/Energy.dat" u 1:(($3) > 0 ? ($3) : NaN) w l lt 1 lw LW dt 1 t "relative error of energy, > 0",\
"" u 1:(($3) < 0 ? -($3) : NaN) w l lt 1 lw LW dt 2 t "< 0",\
"Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/rand01/Energy.dat" u 1:(($4) > 0 ? ($4) : NaN) w l lt 2 lw LW dt 1 t "collisional correction, > 0",\
"" u 1:(($4) < 0 ? -($4) : NaN) w l lw LW lt 2 dt 2 t "< 0"


pause -1

########

set auto

set xl "time [yr]" offset 0,1
set yl "N" offset 2,0
set xr [1E2:1E8]

COLLISION = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/rand01/Collision_time.dat"


set key right top box width 1 spacing 1.0 font "Times-Roman,20"
p COLLISION u 1:2 w l t "N"

pause -1


########

set auto

set xl "time [yr]" offset 0,1
set yl "relative error of energy of planet" offset 2,0
set xr [1E2:1E8]

COLLISION = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/rand01/Collision_time.dat"


set key left top box width 0 spacing 1.0 font "Times-Roman,20"
p COLLISION u 1:(($8) > 0 ? ($8) : NaN) w l lt 1 lw LW dt 1 t "> 0",\
"" u 1:(($8) < 0 ? -($8) : NaN) w l lt 1 lw LW dt 2 t "< 0",\

pause -1


########

set auto

set xl "time [yr]" offset 0,1
set yl "relative error of total L" offset 2,0
set xr [1E2:1E8]

COLLISION = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/rand01/Collision_time.dat"


set key left top box width 0 spacing 1.0 font "Times-Roman,20"
p COLLISION u 1:(($10) > 0 ? ($10) : NaN) w l lt 1 lw LW dt 1 t "> 0",\
"" u 1:(($10) < 0 ? -($10) : NaN) w l lt 1 lw LW dt 2 t "< 0",\

pause -1


########

set auto

set xl "time [yr]" offset 0,1
set yl "relative error of L of planet" offset 2,0
set xr [1E2:1E8]

COLLISION = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/rand01/Collision_time.dat"


set key left top box width 0 spacing 1.0 font "Times-Roman,20"
p COLLISION u 1:(($12) > 0 ? ($12) : NaN) w l lt 1 lw LW dt 1 t "> 0",\
"" u 1:(($12) < 0 ? -($12) : NaN) w l lt 1 lw LW dt 2 t "< 0",\

pause -1



########

set auto

set xl "time [yr]" offset 0,1
set yl "relative collisional correction of L" offset 2,0
set xr [1E2:1E8]

COLLISION = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/rand01/Collision_time.dat"


set key left top box width 0 spacing 1.0 font "Times-Roman,20"
p COLLISION u 1:(($16) > 0 ? ($16) : NaN) w l lt 1 lw LW dt 1 t "x, > 0",\
"" u 1:(($16) < 0 ? -($16) : NaN) w l lt 1 lw LW dt 2 t "x, < 0",\
"" u 1:(($20) > 0 ? ($20) : NaN) w l lt 2 lw LW dt 1 t "y, > 0",\
"" u 1:(($20) < 0 ? -($20) : NaN) w l lt 2 lw LW dt 2 t "y, < 0",\
"" u 1:(($24) > 0 ? ($24) : NaN) w l lt 3 lw LW dt 1 t "z, > 0",\
"" u 1:(($24) < 0 ? -($24) : NaN) w l lt 3 lw LW dt 2 t "z, < 0"

pause -1