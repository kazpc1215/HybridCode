unset multiplot
reset

set term aqua dashed font "Times-Roman,30" enhanced

PAUSE = -1


set bmargin 2
# set lmargin 6
set lmargin 8
set xtics offset 0,0.3
set ytics offset 0.5,0
set xl "time [yr]" offset 0,1

# set log

# set format "10^{%L}"
# unset key

set bar 0
LW_DATA = 3
LW_ANALYTIC = 2
PS = 0.5
DT = 2

###################################

SMALL = "S8E2_t1E9_dtlog_ecc1E-4_Ohtsuki2002/Planetesimal.dat"
LARGE = "S8E2_t1E9_dtlog_ecc1E-4_Ohtsuki2002/Planet.dat"

S8E2 = "S8E2_t1E3_dtlog_M1E24g_ecc1E-4/RMS_randall.dat"

# set yl "ecc, inc [rad]" offset 2,0
set yl "ecc, inc [rad]" offset 3,0
set xr [0:1E3]
set yr [0:0.006]

set key left top box width -1 spacing 1.1 font "Times-Roman,20"
p SMALL u 1:2 w l lw LW_ANALYTIC lt 1 t "<e_S^2>^{1/2}",\
LARGE u 1:2 w l lw LW_ANALYTIC lt 2 t "<e_L^2>^{1/2}",\
SMALL u 1:3 w l lw LW_ANALYTIC lt 3 t "<i_S^2>^{1/2}",\
LARGE u 1:3 w l lw LW_ANALYTIC lt 4 t "<i_L^2>^{1/2}"

pause PAUSE

p S8E2 u 1:2:3 w yerrorlines lw LW_DATA ps PS lt 1 t "<e_S^2>^{1/2}",\
S8E2 u 1:4:5 w yerrorlines lw LW_DATA ps PS lt 2 t "<e_L^2>^{1/2}",\
S8E2 u 1:6:7 w yerrorlines lw LW_DATA ps PS lt 3 t "<i_S^2>^{1/2}",\
S8E2 u 1:8:9 w yerrorlines lw LW_DATA ps PS lt 4 t "<i_L^2>^{1/2}"


pause PAUSE




SMALL = "S1E3_t1E9_dtlog_ecc1E-4_Ohtsuki2002/Planetesimal.dat"

S1E3 = "S1E3_t1E3_dtlog_M1E24g_ecc1E-4/RMS_randall.dat"

# set yl "ecc, inc [rad]" offset 2,0
set yl "ecc, inc [rad]" offset 3,0
set xr [0:1E3]
set yr [0:0.004]
set ytics 0.001

set key left top box width 0 spacing 1.1 font "Times-Roman,20"
p SMALL u 1:2 w l lw LW_ANALYTIC lt 1 t "<e^2>^{1/2}",\
SMALL u 1:3 w l lw LW_ANALYTIC lt 2 t "<i^2>^{1/2}"

pause PAUSE

p S1E3 u 1:2:3 w yerrorlines lw LW_DATA ps PS lt 1 t "<e^2>^{1/2}",\
S1E3 u 1:6:7 w yerrorlines lw LW_DATA ps PS lt 2 t "<i^2>^{1/2}"
