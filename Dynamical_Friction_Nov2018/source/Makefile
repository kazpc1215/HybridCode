# Makefile

# SINGLE, DEBUG, OMP
MODE = SINGLE

# -O0,-O1,-O2,-O3
LEVEL = -O3

CC = gcc
OBJS = hybrid_main.o sub.o acc.o energy.o orbital_elements.o timestep.o neighbor.o massflux.o collision.o

JOBFILE = job.sh

ifeq ($(MODE),OMP)
# OpenMP option
OMPOP	= -fopenmp
#OMPOP	= -h omp
endif


ifeq ($(MODE),SINGLE)
	CFLAGS = $(LEVEL) -I/usr/local/include/ -L/usr/local/lib/openmpi/ -lmpi
ifeq ($(CC),cc)
	CFLAGS += -h noomp
endif
endif

ifeq ($(MODE),DEBUG)
	CFLAGS = -g -O0 -I/usr/local/include/ -L/usr/local/lib/openmpi/ -lmpi
ifeq ($(CC),cc)
	CFLAGS += -h noomp
endif
endif

ifeq ($(MODE),OMP)
	CFLAGS = $(OMPOP) $(LEVEL) -I/usr/local/include/ -L/usr/local/lib/openmpi/ -lmpi
endif

ifeq ($(CC),gcc)
CFLAGS += -lm -Wall -Wextra -MD
endif


define generate_shell_script_1
	echo ''
	echo '-----Macros in Makefile-----'
	echo "MODE = $(MODE)"
	echo "LEVEL = $(LEVEL)"
	echo "CC = $(CC)"
	echo "JOBFILE = $(JOBFILE)"
	echo "OMPOP = $(OMPOP)"
	echo "CFLAGS = $(CFLAGS)"
	echo '----------------------------'

	echo ''
	echo '-----The directory name in hybrid.h-----'
	cat hybrid.h | grep DIRECTORY
	echo '----------------------------------------'

	echo ''
	while [ -z "$$JOBNAME" ]; do \
	read -r -p "Type anything to set the 'job name' in a shell script file. (ex. job01) : " JOBNAME; \
	done; \
	echo "#PBS -N $$JOBNAME" > $(JOBFILE)

	echo '#PBS -l mppwidth=24' >> $(JOBFILE)
	echo '#PBS -j oe' >> $(JOBFILE)
	echo 'cd $$PBS_O_WORKDIR' >> $(JOBFILE)

	if [ $(MODE) = SINGLE ] ; then \
	echo 'export OMP_NUM_THREADS=1' >> $(JOBFILE); \
	elif [ $(MODE) = DEBUG ] ; then \
	echo 'export OMP_NUM_THREADS=1' >> $(JOBFILE); \
	elif [ $(MODE) = OMP ] ; then \
	while [ -z "$$THREADS" ]; do \
	read -r -p "Type any numbers to set OMP_NUM_THREADS. (2,...,24) : " THREADS; \
	done; \
	echo "export OMP_NUM_THREADS=$$THREADS" >> $(JOBFILE); \
	fi

#echo 'export OMP_WAIT_POLICY=PASSIVE' >> $(JOBFILE)
	echo 'export CRAY_OMP_CHECK_AFFINITY=TRUE' >> $(JOBFILE)
endef


define generate_shell_script_2
	echo 'aprun -n 1 -d $$OMP_NUM_THREADS' "./$$EXENAME.exe" '>>' "/work/isoyakz/Dynamical_Friction/log/$$EXENAME.txt" >> $(JOBFILE)
	echo ''
	echo 'The compilation has ended.'

	echo ''
	echo "-----$(JOBFILE)-----"
	cat $(JOBFILE)
	echo '------------------'

	echo ''
	echo 'The shell script file has generated.'
	echo ''
endef


all: $(OBJS)
	@$(call generate_shell_script_1)
	@while [ -z "$$EXENAME" ]; do \
	read -r -p "Type anything to set the 'execution file name'. (*****.exe) : " EXENAME; \
	done; \
	$(CC) $(CFLAGS) -o $$EXENAME.exe $(OBJS); \
	$(call generate_shell_script_2)


.c.o:
	$(CC) $(CFLAGS) -c $<


clean:
	rm -f $(OBJS) *~ *.dat *.d
	rm -f -r *.dtps/ *.dSYM/


ifeq ($(CC),gcc)
-include *.d
else
hybrid_main.o: hybrid_main.c hybrid.h func.h
acc.o: acc.c hybrid.h func.h
collision.o: collision.c hybrid.h func.h
energy.o: energy.c hybrid.h func.h
massflux.o: massflux.c hybrid.h func.h
neighbor.o: neighbor.c hybrid.h func.h
orbital_elements.o: orbital_elements.c hybrid.h func.h
sub.o: sub.c hybrid.h func.h
timestep.o: timestep.c hybrid.h func.h
endif

.PHONY: all clean
