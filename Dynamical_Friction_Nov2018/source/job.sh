#PBS -N test
#PBS -l mppwidth=24
#PBS -j oe
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
export CRAY_OMP_CHECK_AFFINITY=TRUE
aprun -n 1 -d $OMP_NUM_THREADS ./test.exe >> /work/isoyakz/Dynamical_Friction/log/test.txt
