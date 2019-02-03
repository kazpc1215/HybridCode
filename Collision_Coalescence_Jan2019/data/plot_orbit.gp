reset
#set term postscript eps size 800,600 enhanced color font "Helvetica,30"
set term pngcairo size 1200,1200 enhanced color font "Helvetica,35"

set bmargin 4
# set auto

set size ratio -1
set xr [-2.0:2.0]
set yr [-2.0:2.0]


set xl "x [AU]" offset 0,0 tc rgb "white"
set yl "y [AU]" offset 0,0 tc rgb "white"

set xtics offset 0,0 tc rgb "white"
set ytics offset 0,0 tc rgb "white"

set border lw 2 lc rgb "white"

unset key

set object rect behind from screen 0,0 to screen 1,1 fc rgb "black" fill solid 1.0
set object circle fill solid fc rgb "red" at 0, 0

n=2
j=0


# dirname = "N18_t1E8_dtlog_10RHM_2MMSN_Miso_ecc1E-2/rand01"
dirname = "N25_t1E8_dtlog_10RHM_1MMSN_Miso_ecc1E-2/rand01"

timefile = sprintf("%s/Planet01_position.dat",dirname)

N=25

while(n<=75){

TIME = system("cat " . timefile . " | awk \'NR==" . n . "{printf(\"%e\", $1" . ")}\'")
TIME = TIME*1.0



###################################################################################################
set title sprintf("time : %.3e yr",TIME) offset 0,0 tc rgb "white"
set out sprintf("../image/%s/orbit_t%03d.png",dirname,j)


p sprintf("%s/orbit/Planet_orbit.dat.%03d",dirname,n-2) u 2:3 w l lw 3 lc rgb "white",\
for [i=1:N] sprintf("%s/Planet%02d_position.dat",dirname,i) every ::n-2::n-2 u 2:3:(($6)*1.5E3) w circles lc rgb "white" fill solid noborder
###################################################################################################

set output

print n

if(j==0){
n=n+49
}else{
n=n+1
}

j=j+1

}

