#PBS -N test
#PBS -l nodes=1
#PBS -q test-md
#PBS -j oe
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
aprun -n 1 -d $OMP_NUM_THREADS ./N25_t1E8_dtlog_10RHM_1MMSN_Miso_ecc1E-2.exe
