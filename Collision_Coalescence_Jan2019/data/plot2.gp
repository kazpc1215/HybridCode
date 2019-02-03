reset

set term aqua font "Times-Roman,30" enhanced
#set term pngcairo size 800,600 enhanced color font "Helvetica,30"

set auto

set xr [1E6:1E8]
YEAR=1.0
set xl "time [yr]" offset 0,1

#set xr [1E-2:1E2]
#YEAR=1.0E6
#set xl "time [Myr]" offset 0,1


set log x
set format x "10^{%L}"
#set format x "%g"

set bmargin 2
set lmargin 7
unset key

set xtics offset 0,0.3
set ytics offset 0,0

#set xl "time [kyr]" offset 0,1
#set yl "semi-major axis [AU]" offset 2.5,0
set yl "pericenter/apocenter distance [AU]" offset 2.5,0


# dirname = "N18_t1E8_dtlog_10RHM_2MMSN_Miso_ecc1E-2/rand01"
# i1 = 1
# i2 = 14
# i3 = 16
# i4 = 17


dirname = "N25_t1E8_dtlog_10RHM_1MMSN_Miso_ecc1E-2/rand01"
 i1 = 1
 i2 = 2
 i3 = 21
 i4 = 23
 i5 = 25

N=25

#YEAR=1.0E3
EVERY=1
LW=2


# plot for [i=1:N] sprintf("%s/Planet%02d.dat",dirname,i) every EVERY u (($1)/YEAR):(($3)*(1.0-($2))) w l lw 0.5

# pause
set yr [0:3]

plot for [i=1:N] sprintf("%s/Planet%02d.dat",dirname,i) every EVERY u (($1)/YEAR):(($3)*(1.0-($2))) w l lw 0.5 lc rgb "gray50" notitle,\
for [i=1:N] sprintf("%s/Planet%02d.dat",dirname,i) every EVERY u (($1)/YEAR):(($3)*(1.0+($2))) w l lw 0.5 lc rgb "gray50" notitle,\
sprintf("%s/Planet%02d.dat",dirname,i1) every EVERY u (($1)/YEAR):(($3)*(1.0-($2))) w l lw LW lt 1 notitle,\
sprintf("%s/Planet%02d.dat",dirname,i1) every EVERY u (($1)/YEAR):(($3)*(1.0+($2))) w l lw LW lt 1 notitle,\
sprintf("%s/Planet%02d.dat",dirname,i2) every EVERY u (($1)/YEAR):(($3)*(1.0-($2))) w l lw LW lt 2 notitle,\
sprintf("%s/Planet%02d.dat",dirname,i2) every EVERY u (($1)/YEAR):(($3)*(1.0+($2))) w l lw LW lt 2 notitle,\
sprintf("%s/Planet%02d.dat",dirname,i3) every EVERY u (($1)/YEAR):(($3)*(1.0-($2))) w l lw LW lt 3 notitle,\
sprintf("%s/Planet%02d.dat",dirname,i3) every EVERY u (($1)/YEAR):(($3)*(1.0+($2))) w l lw LW lt 3 notitle,\
sprintf("%s/Planet%02d.dat",dirname,i4) every EVERY u (($1)/YEAR):(($3)*(1.0-($2))) w l lw LW lt 4 notitle,\
sprintf("%s/Planet%02d.dat",dirname,i4) every EVERY u (($1)/YEAR):(($3)*(1.0+($2))) w l lw LW lt 4 notitle

# sprintf("%s/Planet%02d.dat",dirname,i5) every EVERY u (($1)/YEAR):(($3)*(1.0-($2))) w l lw LW lt 5 notitle,\
# sprintf("%s/Planet%02d.dat",dirname,i5) every EVERY u (($1)/YEAR):(($3)*(1.0+($2))) w l lw LW lt 5 notitle


# pause

set xr [1E5:1E8]

# set yr [0:20]
set yr [0:25]
set yl "number of protoplanets" offset 2.5,0
plot sprintf("%s/WeightedAverage.dat",dirname) every EVERY u (($1)/YEAR):5 w l lw LW notitle

pause -1


set log
set format "10^{%L}"

set xr [1E5:1E8]
set yr [1E-3:1]



# set yr [0:0.2]
# set ytics 0,0.05
# set mytics 5

set yl "mass weighted mean ecc" offset 2,0
plot sprintf("%s/WeightedAverage.dat",dirname) every EVERY u (($1)/YEAR):2 w l lw LW notitle

pause -1

set yl "mass weighted mean inc [rad]" offset 2,0
plot sprintf("%s/WeightedAverage.dat",dirname) every EVERY u (($1)/YEAR):3 w l lw LW notitle

pause -1

set yl "mass weighted mean (e^2 + i^2)^{1/2}" offset 2,0
plot sprintf("%s/WeightedAverage.dat",dirname) every EVERY u (($1)/YEAR):4 w l lw LW notitle