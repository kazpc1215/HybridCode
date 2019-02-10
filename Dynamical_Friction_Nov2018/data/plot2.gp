reset

set term aqua dashed font "Times-Roman,30" enhanced


# set bmargin 2
set lmargin 6
set xtics offset 0,0.3
set ytics offset 0.5,0
set xl "time [yr]" offset 0,1


set bar 0.3

PS = 1.0

set log

unset key


##### ecc = 0.05 #####

NOFRAG_ecc = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_nofrag_acc/RMS_randall.dat"
FRAG19_ecc = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_frag_acc/RMS_randall.dat"

NOFRAG_sigma = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_nofrag_acc/tracerlistnumber_all.dat"
FRAG19_sigma = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_frag_acc/Sigma_randall.dat"



###
set xl "Time [yr]" offset 0,1
set xr [10:1000]
set format "10^{%L}"

set multiplot layout 2,1 upwards

set bmargin 2
set tmargin 0
set yl "ecc" offset 2,0
set yr [0.005:0.1]

p NOFRAG_ecc u 1:4:5 w yerrorlines lw 2 dt 1 lt 2 lc rgb "dark-blue" ps PS notitle,\
NOFRAG_ecc u 1:10:11 w yerrorlines lw 2 dt 2 lt 2 lc rgb "dark-blue" ps PS

###

set bmargin 0.5
set tmargin 1
unset xl
set format x ""
set yl "M_{center} [M_{earth}]" offset 1.0,0
set auto y

p NOFRAG_sigma u 1:(($4)*3.0E-8/3.0E-6):(($5)*3.0E-8/3.0E-6) w yerrorlines lw 2 dt 1 lt 2 lc rgb "dark-blue" ps PS notitle


unset multiplot
###



pause -1



###
set xl "Time [yr]" offset 0,1
set xr [10:1000]
set format "10^{%L}"

set multiplot layout 2,1 upwards

set bmargin 2
set tmargin 0
set yl "ecc" offset 2,0
# set yr [0.005:0.1]

p FRAG19_ecc u 1:4:5 w yerrorlines lw 2 dt 1 lt 5 lc rgb "red" ps PS notitle,\
FRAG19_ecc u 1:10:11 w yerrorlines lw 2 dt 2 lt 5 lc rgb "red" ps PS

###

set bmargin 0.5
set tmargin 1
unset xl
set format x ""
set yl "M_{center} [M_{earth}]" offset 2.0,0
# set yr [0.1:20]


p FRAG19_sigma u 1:(($4)*3.0E-8/3.0E-6*1.079105263157895e+03):(($5*3.0E-8/3.0E-6*1.079105263157895e+03)) w yerrorlines lw 2 dt 1 lt 5 lc rgb "red" ps PS notitle

unset multiplot
###




pause -1



###
set xl "Time [yr]" offset 0,1
set xr [10:1000]
set format "10^{%L}"

set multiplot layout 2,1 upwards

set bmargin 2
set tmargin 0
set yl "ecc" offset 2,0
set yr [0.005:0.1]

p NOFRAG_ecc u 1:4:5 w yerrorlines lw 2 dt 1 lt 2 lc rgb "dark-blue" ps PS notitle,\
FRAG19_ecc u 1:4:5 w yerrorlines lw 2 dt 1 lt 5 lc rgb "red" ps PS notitle

###

set bmargin 0.5
set tmargin 1
unset xl
set format x ""
set yl "M_{center} [M_{earth}]" offset 2.0,0
set auto y


p NOFRAG_sigma u 1:(($4)*3.0E-8/3.0E-6):(($5)*3.0E-8/3.0E-6) w yerrorlines lw 2 dt 1 lt 2 lc rgb "dark-blue" ps PS notitle,\
FRAG19_sigma u 1:(($4)*3.0E-8/3.0E-6*1.079105263157895e+03):(($5*3.0E-8/3.0E-6*1.079105263157895e+03)) w yerrorlines lw 2 dt 1 lt 5 lc rgb "red" ps PS notitle

unset multiplot
###