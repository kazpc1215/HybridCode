unset multiplot
reset

set term aqua dashed font "Times-Roman,30" enhanced

PAUSE = -1

set bmargin 2
set lmargin 6
set xtics offset 0,0.3
set ytics offset 0.5,0
set xl "time [yr]" offset 0,1
set format "10^{%L}"

set bar 0.3

set log

# unset key


###########
ECC = 0.05
###########


#########################################



if (ECC == 0.01) {
 OHTSUKI_NOFRAG_PLANET = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc1E-2_nofrag_dt/Planet.dat"
 OHTSUKI_NOFRAG_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc1E-2_nofrag_dt/Planetesimal.dat"
 OHTSUKI_FRAG_1E19_PLANET = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc1E-2_frag_dt/Planet.dat"
 OHTSUKI_FRAG_1E19_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc1E-2_frag_dt/Planetesimal.dat"
 OHTSUKI_FRAG_1E16_PLANET = "Meach3E-8_Mtot3E-5_Mmax5E-18_t1E9_dtlog_ecc1E-2_frag_dt/Planet.dat"
 OHTSUKI_FRAG_1E16_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-18_t1E9_dtlog_ecc1E-2_frag_dt/Planetesimal.dat"
 NOFRAG = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc1E-2_nofrag_acc/RMS_randall.dat"
 FRAG_1E19 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc1E-2_frag_acc/RMS_randall.dat"
 FRAG_1E16 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc1E-2_frag_acc/RMS_randall.dat"
}

if (ECC == 0.03) {
 OHTSUKI_NOFRAG_PLANET = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc3E-2_nofrag_dt/Planet.dat"
 OHTSUKI_NOFRAG_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc3E-2_nofrag_dt/Planetesimal.dat"
 OHTSUKI_FRAG_1E19_PLANET = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc3E-2_frag_dt/Planet.dat"
 OHTSUKI_FRAG_1E19_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc3E-2_frag_dt/Planetesimal.dat"
 OHTSUKI_FRAG_1E16_PLANET = "Meach3E-8_Mtot3E-5_Mmax5E-18_t1E9_dtlog_ecc3E-2_frag_dt/Planet.dat"
 OHTSUKI_FRAG_1E16_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-18_t1E9_dtlog_ecc3E-2_frag_dt/Planetesimal.dat"
 NOFRAG = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc3E-2_nofrag_acc/RMS_randall.dat"
 FRAG_1E19 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc3E-2_frag_acc/RMS_randall.dat"
 FRAG_1E16 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc3E-2_frag_acc/RMS_randall.dat"
}

if (ECC == 0.05) {
 OHTSUKI_NOFRAG_PLANET = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc5E-2_nofrag_dt/Planet.dat"
 OHTSUKI_NOFRAG_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc5E-2_nofrag_dt/Planetesimal.dat"
 OHTSUKI_FRAG_1E19_PLANET = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc5E-2_frag_dt/Planet.dat"
 OHTSUKI_FRAG_1E19_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc5E-2_frag_dt/Planetesimal.dat"
 OHTSUKI_FRAG_1E16_PLANET = "Meach3E-8_Mtot3E-5_Mmax5E-18_t1E9_dtlog_ecc5E-2_frag_dt/Planet.dat"
 OHTSUKI_FRAG_1E16_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-18_t1E9_dtlog_ecc5E-2_frag_dt/Planetesimal.dat"
 NOFRAG = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_nofrag_acc/RMS_randall.dat"
 FRAG_1E19 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_frag_acc/RMS_randall.dat"
 FRAG_1E16 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc5E-2_frag_acc/RMS_randall.dat"
}

set yl "ecc" offset 2,0
set xr [1:1000]
if (ECC == 0.01) {
 set yr [0.001:0.1]
}
if (ECC == 0.03) {
 set yr [0.003:0.3]
}
if (ECC == 0.05) {
 set yr [0.005:0.5]
}


### nofrag_ecc  ###
set multiplot
set key left Left bottom box width -10 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_NOFRAG_PLANET u 1:2 w l lt -1 lw 2 dt 1 t "Planet, semi-analytic, no frag.",\
NOFRAG u 1:4:5 w yerrorlines lw 2 dt 1 lt 1 ps 0.5 t "Planet, r.m.s., no frag."
set key left Left top box width -12 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_NOFRAG_PLANETESIMAL u 1:2 w l lt -1 lw 2 dt 2 t "Planetesimal, semi-analytic, no frag.",\
NOFRAG u 1:10:11 w yerrorlines lw 2 dt 2 lt 1 ps 0.5 t "Planetesimal, r.m.s., no frag."
unset multiplot

pause PAUSE



### frag_1E19_ecc  ###
set multiplot
set key left Left bottom box width -10 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_FRAG_1E19_PLANET u 1:2 w l lt -1 lw 2 dt 1 t "Planet, semi-analytic, frag., 10^{19}g",\
FRAG_1E19 u 1:4:5 w yerrorlines lw 2 dt 1 lt 2 ps 0.5 t "Planet, r.m.s., frag., 10^{19}g"
set key left Left top box width -12 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_FRAG_1E19_PLANETESIMAL u 1:2 w l lt -1 lw 2 dt 2 t "Planetesimal, semi-analytic, frag., 10^{19}g",\
FRAG_1E19 u 1:10:11 w yerrorlines lw 2 dt 2 lt 2 ps 0.5 t "Planetesimal, r.m.s., frag., 10^{19}g"
unset multiplot

pause PAUSE



### frag_1E16_ecc  ###
set multiplot
set key left Left bottom box width -10 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_FRAG_1E16_PLANET u 1:2 w l lt -1 lw 2 dt 1 t "Planet, semi-analytic, frag., 10^{16}g",\
FRAG_1E16 u 1:4:5 w yerrorlines lw 2 dt 1 lt 3 ps 0.5 t "Planet, r.m.s., frag., 10^{16}g"
set key left Left top box width -12 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_FRAG_1E16_PLANETESIMAL u 1:2 w l lt -1 lw 2 dt 2 t "Planetesimal, semi-analytic, frag., 10^{16}g",\
FRAG_1E16 u 1:10:11 w yerrorlines lw 2 dt 2 lt 3 ps 0.5 t "Planetesimal, r.m.s., frag., 10^{16}g"
unset multiplot

pause PAUSE



### numerical_ecc  ###
set multiplot
set key left Left bottom box width -12 spacing 1.1 reverse font "Times-Roman,16"
plot NOFRAG u 1:4:5 w yerrorlines lw 2 dt 1 lt 1 ps 0.5 t "Planet, r.m.s., no frag.",\
FRAG_1E19 u 1:4:5 w yerrorlines lw 2 dt 1 lt 2 ps 0.5 t "Planet, r.m.s.,      frag., 10^{19}g",\
FRAG_1E16 u 1:4:5 w yerrorlines lw 2 dt 1 lt 3 ps 0.5 t "Planet, r.m.s.,      frag., 10^{16}g"
set key left Left top box width -13.5 spacing 1.1 reverse font "Times-Roman,16"
plot NOFRAG u 1:10:11 w yerrorlines lw 2 dt 2 lt 1 ps 0.5 t "Planetesimal, r.m.s., no frag.",\
FRAG_1E19 u 1:10:11 w yerrorlines lw 2 dt 2 lt 2 ps 0.5 t "Planetesimal, r.m.s.,      frag., 10^{19}g",\
FRAG_1E16 u 1:10:11 w yerrorlines lw 2 dt 2 lt 3 ps 0.5 t "Planetesimal, r.m.s.,      frag., 10^{16}g"
unset multiplot

pause PAUSE



### analytic_ecc  ###
set multiplot
set key left Left bottom box width -13.5 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_NOFRAG_PLANET u 1:2 w l lt 1 lw 2 dt 1 t "Planet, semi-analytic, no frag.",\
OHTSUKI_FRAG_1E19_PLANET u 1:2 w l lt 2 lw 2 dt 1 t "Planet, semi-analytic,      frag., 10^{19}g",\
OHTSUKI_FRAG_1E16_PLANET u 1:2 w l lt 3 lw 2 dt 1 t "Planet, semi-analytic,      frag., 10^{16}g"
set key left Left top box width -15.5 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_NOFRAG_PLANETESIMAL u 1:2 w l lt 1 lw 2 dt 2 t "Planetesimal, semi-analytic, no frag.",\
OHTSUKI_FRAG_1E19_PLANETESIMAL u 1:2 w l lt 2 lw 2 dt 2 t "Planetesimal, semi-analytic,      frag., 10^{19}g",\
OHTSUKI_FRAG_1E16_PLANETESIMAL u 1:2 w l lt 3 lw 2 dt 2 t "Planetesimal, semi-analytic,      frag., 10^{16}g"
unset multiplot

pause PAUSE

#####################################

set yl "inc [rad]" offset 2,0
set xr [1:1000]
if (ECC == 0.01) {
 set yr [0.0005:0.05]
}
if (ECC == 0.03) {
 set yr [0.0015:0.15]
}
if (ECC == 0.05) {
 set yr [0.0025:0.25]
}

### nofrag_inc  ###
set multiplot
set key left Left bottom box width -10 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_NOFRAG_PLANET u 1:3 w l lt -1 lw 2 dt 1 t "Planet, semi-analytic, no frag.",\
NOFRAG u 1:14:15 w yerrorlines lw 2 dt 1 lt 1 ps 0.5 t "Planet, r.m.s., no frag."
set key left Left top box width -12 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_NOFRAG_PLANETESIMAL u 1:3 w l lt -1 lw 2 dt 2 t "Planetesimal, semi-analytic, no frag.",\
NOFRAG u 1:20:21 w yerrorlines lw 2 dt 2 lt 1 ps 0.5 t "Planetesimal, r.m.s., no frag."
unset multiplot

pause PAUSE



### frag_1E19_inc  ###
set multiplot
set key left Left bottom box width -10 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_FRAG_1E19_PLANET u 1:3 w l lt -1 lw 2 dt 1 t "Planet, semi-analytic, frag., 10^{19}g",\
FRAG_1E19 u 1:14:15 w yerrorlines lw 2 dt 1 lt 2 ps 0.5 t "Planet, r.m.s., frag., 10^{19}g"
set key left Left top box width -12 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_FRAG_1E19_PLANETESIMAL u 1:3 w l lt -1 lw 2 dt 2 t "Planetesimal, semi-analytic, frag., 10^{19}g",\
FRAG_1E19 u 1:20:21 w yerrorlines lw 2 dt 2 lt 2 ps 0.5 t "Planetesimal, r.m.s., frag., 10^{19}g"
unset multiplot

pause PAUSE



### frag_1E16_inc  ###
set multiplot
set key left Left bottom box width -10 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_FRAG_1E16_PLANET u 1:3 w l lt -1 lw 2 dt 1 t "Planet, semi-analytic, frag., 10^{16}g",\
FRAG_1E16 u 1:14:15 w yerrorlines lw 2 dt 1 lt 3 ps 0.5 t "Planet, r.m.s., frag., 10^{16}g"
set key left Left top box width -12 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_FRAG_1E16_PLANETESIMAL u 1:3 w l lt -1 lw 2 dt 2 t "Planetesimal, semi-analytic, frag., 10^{16}g",\
FRAG_1E16 u 1:20:21 w yerrorlines lw 2 dt 2 lt 3 ps 0.5 t "Planetesimal, r.m.s., frag., 10^{16}g"
unset multiplot

pause PAUSE



### numerical_inc  ###
set multiplot
set key left Left bottom box width -12 spacing 1.1 reverse font "Times-Roman,16"
plot NOFRAG u 1:14:15 w yerrorlines lw 2 dt 1 lt 1 ps 0.5 t "Planet, r.m.s., no frag.",\
FRAG_1E19 u 1:14:15 w yerrorlines lw 2 dt 1 lt 2 ps 0.5 t "Planet, r.m.s.,      frag., 10^{19}g",\
FRAG_1E16 u 1:14:15 w yerrorlines lw 2 dt 1 lt 3 ps 0.5 t "Planet, r.m.s.,      frag., 10^{16}g"
set key left Left top box width -13.5 spacing 1.1 reverse font "Times-Roman,16"
plot NOFRAG u 1:20:21 w yerrorlines lw 2 dt 2 lt 1 ps 0.5 t "Planetesimal, r.m.s., no frag.",\
FRAG_1E19 u 1:20:21 w yerrorlines lw 2 dt 2 lt 2 ps 0.5 t "Planetesimal, r.m.s.,      frag., 10^{19}g",\
FRAG_1E16 u 1:20:21 w yerrorlines lw 2 dt 2 lt 3 ps 0.5 t "Planetesimal, r.m.s.,      frag., 10^{16}g"
unset multiplot

pause PAUSE



### analytic_inc  ###
set multiplot
set key left Left bottom box width -13.5 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_NOFRAG_PLANET u 1:3 w l lt 1 lw 2 dt 1 t "Planet, semi-analytic, no frag.",\
OHTSUKI_FRAG_1E19_PLANET u 1:3 w l lt 2 lw 2 dt 1 t "Planet, semi-analytic,      frag., 10^{19}g",\
OHTSUKI_FRAG_1E16_PLANET u 1:3 w l lt 3 lw 2 dt 1 t "Planet, semi-analytic,      frag., 10^{16}g"
set key left Left top box width -15.5 spacing 1.1 reverse font "Times-Roman,16"
plot OHTSUKI_NOFRAG_PLANETESIMAL u 1:3 w l lt 1 lw 2 dt 2 t "Planetesimal, semi-analytic, no frag.",\
OHTSUKI_FRAG_1E19_PLANETESIMAL u 1:3 w l lt 2 lw 2 dt 2 t "Planetesimal, semi-analytic,      frag., 10^{19}g",\
OHTSUKI_FRAG_1E16_PLANETESIMAL u 1:3 w l lt 3 lw 2 dt 2 t "Planetesimal, semi-analytic,      frag., 10^{16}g"
unset multiplot

pause PAUSE



###############################

set auto
unset log
set log x
set format y "%g"

if (ECC == 0.01) {
 NOFRAG = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc1E-2_nofrag_acc/axis_evo.dat"
 FRAG_1E19 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc1E-2_frag_acc/axis_evo.dat"
 FRAG_1E16 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc1E-2_frag_acc/axis_evo.dat"
}
if (ECC == 0.03) {
 NOFRAG = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc3E-2_nofrag_acc/axis_evo.dat"
 FRAG_1E19 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc3E-2_frag_acc/axis_evo.dat"
 FRAG_1E16 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc3E-2_frag_acc/axis_evo.dat"
}
if (ECC == 0.05) {
 NOFRAG = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_nofrag_acc/axis_evo.dat"
 FRAG_1E19 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_frag_acc/axis_evo.dat"
 FRAG_1E16 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc5E-2_frag_acc/axis_evo.dat"
}

set xr [1:1000]
set yr [-5:5]

RHILL_1 = 8.814745316374441e-03
RHILL_2 = 0.01
RHILL_3 = 1.134462725930811e-02


set yl "normalized migration distance" offset 2,0
set key left top box width -5 spacing 1.1 font "Times-Roman,16"
plot NOFRAG u 1:($8/RHILL_1):($9/RHILL_1) w yerrorlines lw 2 lt 1 ps 0.5 t "Planet1, no frag.",\
NOFRAG u 1:($10/RHILL_2):($11/RHILL_2) w yerrorlines lw 2 lt 2 ps 0.5 t "Planet2, no frag.",\
NOFRAG u 1:($12/RHILL_3):($13/RHILL_3) w yerrorlines lw 2 lt 3 ps 0.5 t "Planet3, no frag."

pause PAUSE

set yl "normalized migration distance" offset 2,0
set key left top box width -5 spacing 1.1 font "Times-Roman,16"
plot FRAG_1E19 u 1:($8/RHILL_1):($9/RHILL_1) w yerrorlines lw 2 lt 1 ps 0.5 t "Planet1, frag., 10^{19}g",\
FRAG_1E19 u 1:($10/RHILL_2):($11/RHILL_2) w yerrorlines lw 2 lt 2 ps 0.5 t "Planet2, frag., 10^{19}g",\
FRAG_1E19 u 1:($12/RHILL_3):($13/RHILL_3) w yerrorlines lw 2 lt 3 ps 0.5 t "Planet3, frag., 10^{19}g"

pause PAUSE

set yl "normalized migration distance" offset 2,0
set key left top box width -5 spacing 1.1 font "Times-Roman,16"
plot FRAG_1E16 u 1:($8/RHILL_1):($9/RHILL_1) w yerrorlines lw 2 lt 1 ps 0.5 t "Planet1, frag., 10^{16}g",\
FRAG_1E16 u 1:($10/RHILL_2):($11/RHILL_2) w yerrorlines lw 2 lt 2 ps 0.5 t "Planet2, frag., 10^{16}g",\
FRAG_1E16 u 1:($12/RHILL_3):($13/RHILL_3) w yerrorlines lw 2 lt 3 ps 0.5 t "Planet3, frag., 10^{16}g"

pause PAUSE





set yl "semi-major axis [AU]" offset 3,0
set yr [0.7:1.3]
set key left top box width -5 spacing 1.1 font "Times-Roman,16"
plot NOFRAG u 1:2:3 w yerrorlines lw 2 lt 1 ps 0.5 t "Planet1, no frag.",\
NOFRAG u 1:4:5 w yerrorlines lw 2 lt 2 ps 0.5 t "Planet2, no frag.",\
NOFRAG u 1:6:7 w yerrorlines lw 2 lt 3 ps 0.5 t "Planet3, no frag."

pause PAUSE

set yl "semi-major axis [AU]" offset 3,0
set yr [0.7:1.3]
set key left top box width -5 spacing 1.1 font "Times-Roman,16"
plot FRAG_1E19 u 1:2:3 w yerrorlines lw 2 lt 1 ps 0.5 t "Planet1, frag., 10^{19}g",\
FRAG_1E19 u 1:4:5 w yerrorlines lw 2 lt 2 ps 0.5 t "Planet2, frag., 10^{19}g",\
FRAG_1E19 u 1:6:7 w yerrorlines lw 2 lt 3 ps 0.5 t "Planet3, frag., 10^{19}g"

pause PAUSE

set yl "semi-major axis [AU]" offset 3,0
set yr [0.7:1.3]
set key left top box width -5 spacing 1.1 font "Times-Roman,16"
plot FRAG_1E16 u 1:2:3 w yerrorlines lw 2 lt 1 ps 0.5 t "Planet1, frag., 10^{16}g",\
FRAG_1E16 u 1:4:5 w yerrorlines lw 2 lt 2 ps 0.5 t "Planet2, frag., 10^{16}g",\
FRAG_1E16 u 1:6:7 w yerrorlines lw 2 lt 3 ps 0.5 t "Planet3, frag., 10^{16}g"

pause PAUSE



###############################

### sigma ###

set log
set format "10^{%L}"
set auto
set xr [0.1:1000]
#set yr [:1]
set yr [:1.01]
set yl "{/Symbol S}/{/Symbol S}_0" offset 2,0

#####

if (ECC == 0.01) {
 OHTSUKI_FRAG_1E19_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc1E-2_frag_dt/Planetesimal.dat"
 OHTSUKI_FRAG_1E16_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-18_t1E9_dtlog_ecc1E-2_frag_dt/Planetesimal.dat"
 FRAG_1E19 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc1E-2_frag_acc/Sigma_randall.dat"
 FRAG_1E16 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc1E-2_frag_acc/Sigma_randall.dat"
}
if (ECC == 0.03) {
 OHTSUKI_FRAG_1E19_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc3E-2_frag_dt/Planetesimal.dat"
 OHTSUKI_FRAG_1E16_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-18_t1E9_dtlog_ecc3E-2_frag_dt/Planetesimal.dat"
 FRAG_1E19 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc3E-2_frag_acc/Sigma_randall.dat"
 FRAG_1E16 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc3E-2_frag_acc/Sigma_randall.dat"
}
if (ECC == 0.05) {
 OHTSUKI_FRAG_1E19_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-15_t1E9_dtlog_ecc5E-2_frag_dt/Planetesimal.dat"
 OHTSUKI_FRAG_1E16_PLANETESIMAL = "Meach3E-8_Mtot3E-5_Mmax5E-18_t1E9_dtlog_ecc5E-2_frag_dt/Planetesimal.dat"
 FRAG_1E19 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_frag_acc/Sigma_randall.dat"
 FRAG_1E16 = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc5E-2_frag_acc/Sigma_randall.dat"
}


set key left Left bottom box width -5 spacing 1.2 reverse font "Times-Roman,16"
plot FRAG_1E19 u 1:6:7 w yerrorlines lw 2 lt 2 ps 0.5 t "numerical, m_{max}=10^{19}g",\
FRAG_1E16 u 1:6:7 w yerrorlines lw 2 lt 3 ps 0.5 t "numerical, m_{max}=10^{16}g",\
OHTSUKI_FRAG_1E19_PLANETESIMAL u 1:($5/3.0E-5) w l lt 2 lw 2 dt 2 t "semi-analytic, m_{max}=10^{19}g",\
OHTSUKI_FRAG_1E16_PLANETESIMAL u 1:($5/3.0E-5) w l lt 3 lw 2 dt 2 t "semi-analytic, m_{max}=10^{16}g"