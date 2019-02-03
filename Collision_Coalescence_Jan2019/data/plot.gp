reset
#set term postscript eps size 800,600 enhanced color font "Helvetica,30"
set term pngcairo size 800,600 enhanced color font "Helvetica,20"
# set term pngcairo size 800,600 enhanced color font "Times-Roman,20"

#set key left reverse box spacing 1.5 width -2
unset key
set bmargin 4
set auto

# set xr [0.3:2.0]
set xr [0.3:2.5]


set log xy
set xl "semi-major axis [AU]"
set format y "10^{%L}"

# set xtics add (0.5,1.5,2.0)
set xtics add (0.5,1.5,2.0,2.5)

n=2
j=0


dirname = "N18_t1E8_dtlog_10RHM_2MMSN_Miso_ecc1E-2/rand01"
# dirname = "N25_t1E8_dtlog_10RHM_1MMSN_Miso_ecc1E-2/rand01"

timefile = sprintf("%s/Planet01.dat",dirname)


N=18

while(n<=316){

TIME = system("cat " . timefile . " | awk \'NR==" . n . "{printf(\"%e\", $1" . ")}\'")
TIME = TIME*1.0


###################################################################################################

#ecc
set title sprintf("time : %.3e yr",TIME) offset 0,0
set out sprintf("../image/%s/ecc_axis_t%02d.png",dirname,j)

set yl "ecc" offset 0,0
set yr [1E-3:1]


plot for [i=1:N] sprintf("%s/Planet%02d.dat",dirname,i) every ::n-2::n-2 u 3:2:(($9)*8.0E2*(($3))) w circles lc rgb "gray" fill solid border lc rgb "black" notitle

set output


#mass log
set title sprintf("time : %.3e yr",TIME) offset 0,0
set out sprintf("../image/%s/mass_axis_t%02d.png",dirname,j)

set yl "mass [M_E]" offset 0,0
#set yr [1E-2:1]
#set yr [1E-1:1E1]
set yr [1E-2:1E1]

plot for [i=1:N] sprintf("./%s/Planet%02d.dat",dirname,i) every ::n-2::n-2 u 3:(($10)/3.0E-6):(($9)*8.0E2*(($3))) w circles lc rgb "gray" fill solid border lc rgb "black" notitle

set output
###################################################################################################

print n

if(j==0){
  n=n+49
}else{
  if(n<59){
    n=n+1
  }else{
    n=n+16
  }
}

j=j+1

}


