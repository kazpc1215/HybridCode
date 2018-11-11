# Makefile

# SINGLE, DEBUG, OMP
MODE = SINGLE

# -O0,-O1,-O2,-O3
LEVEL = -O3

CC = gcc
OBJS = hybrid_main.o sub.o acc.o energy.o orbital_elements.o timestep.o neighbor.o massflux.o collision.o SFMT.o


ifeq ($(MODE),OMP)
# OpenMP option
# OMPOP = -fopenmp
OMPOP = -h omp
endif


ifeq ($(MODE),SINGLE)
	CFLAGS = $(LEVEL)
ifeq ($(CC),cc)
	CFLAGS += -h noomp
endif
endif

ifeq ($(MODE),DEBUG)
	CFLAGS = -g -O0
ifeq ($(CC),cc)
	CFLAGS += -h noomp
endif
endif

ifeq ($(MODE),OMP)
	CFLAGS = $(OMPOP) $(LEVEL)
endif

CFLAGS += -DSFMT_MEXP=216091

ifeq ($(CC),gcc)
CFLAGS += -I/usr/local/include/ -L/usr/local/lib/openmpi/ -lmpi -lm -Wall -Wextra -MD
endif


define show_macro
	echo ''; \
	echo '-----Macros in Makefile-----'; \
	echo "MODE = $(MODE)"; \
	echo "LEVEL = $(LEVEL)"; \
	echo "CC = $(CC)"; \
	echo "OMPOP = $(OMPOP)"; \
	echo "CFLAGS = $(CFLAGS)"; \
	echo '----------------------------'; \
	echo ''; \
	echo '-----The directory name in hybrid.h-----'; \
	cat hybrid.h | grep DIRECTORY; \
	echo '----------------------------------------'; \
	echo ''
endef

define generate_shell_script_1
	echo "#PBS -N $$JOBFILE" > "$$JOBFILE.sh"; \
	while [ -z "$$NUMNODES" ]; do \
        read -e -r -p "Type any numbers to set number of nodes. : " NUMNODES; \
        echo "$$NUMNODES"; \
	hitsory -s "$$NUMNODES"; \
        done; \
	echo "#PBS -l nodes=$$NUMNODES" >> "$$JOBFILE.sh"; \
	while [ -z "$$QUEUENAME" ]; do \
	read -e -r -p "Select the queue name. ('large', 'bulk', 'test', 'debug') (*****-md) : " QUEUENAME; \
	echo "$$QUEUENAME"; \
	history -s "$$QUEUENAME"; \
	done; \
	echo "#PBS -q $$QUEUENAME-md" >> "$$JOBFILE.sh"; \
	echo '#PBS -j oe' >> "$$JOBFILE.sh"; \
	echo 'cd $$PBS_O_WORKDIR' >> "$$JOBFILE.sh"; \
	if [ $(MODE) = SINGLE ] ; then \
	echo 'export OMP_NUM_THREADS=1' >> "$$JOBFILE.sh"; \
	elif [ $(MODE) = DEBUG ] ; then \
	echo 'export OMP_NUM_THREADS=1' >> "$$JOBFILE.sh"; \
	elif [ $(MODE) = OMP ] ; then \
	while [ -z "$$THREADS" ]; do \
	read -e -r -p "Type any numbers to set OMP_NUM_THREADS. (1,...,40) : " THREADS; \
	echo "$$THREADS"; \
	history -s "$$THREADS"; \
	done; \
	echo "export OMP_NUM_THREADS=$$THREADS" >> "$$JOBFILE.sh"; \
	fi
endef


define generate_shell_script_2
	echo 'aprun -n' "$$NUMPROCS" '-d $$OMP_NUM_THREADS' "./$$EXENAME.exe" >> "$$JOBFILE.sh"; \
	echo ''; \
	echo 'The compilation has ended.'; \
	echo ''; \
	echo "-----$$JOBFILE.sh-----"; \
	cat "$$JOBFILE.sh"; \
	echo '------------------'; \
	echo ''; \
	echo 'The shell script file has generated.'; \
	echo ''; \
	history -w .input_history
endef


all: $(OBJS)
	@$(call show_macro)
	@history -r .input_history; \
	while [ -z "$$JOBFILE" ]; do \
	read -e -r -p "Type anything to set the 'shell script file name'. (*****.sh) : " JOBFILE; \
	echo "$$JOBFILE"; \
	history -s "$$JOBFILE"; \
	done; \
	$(call generate_shell_script_1); \
	while [ -z "$$NUMPROCS" ]; do \
        read -e -r -p "Type any numbers to set MPI process number. (1,...,40 * numnodes) : " NUMPROCS; \
        echo "$$NUMPROCS"; \
	history -s "$$NUMPROCS"; \
        done; \
	while [ -z "$$EXENAME" ]; do \
	read -e -r -p "Type anything to set the 'execution file name'. (*****.exe) : " EXENAME; \
	echo "$$EXENAME"; \
	history -s "$$EXENAME"; \
	done; \
	$(CC) $(CFLAGS) -o $$EXENAME.exe $(OBJS); \
	$(call generate_shell_script_2)


.c.o:
	$(CC) $(CFLAGS) -c $<


clean:
	rm -f $(OBJS) *~ *.dat *.d



ifeq ($(CC),gcc)
-include *.d
else
hybrid_main.o: hybrid_main.c hybrid.h func.h SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
acc.o: acc.c hybrid.h func.h SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
collision.o: collision.c hybrid.h func.h SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
energy.o: energy.c hybrid.h func.h SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
massflux.o: massflux.c hybrid.h func.h SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
neighbor.o: neighbor.c hybrid.h func.h SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
orbital_elements.o: orbital_elements.c hybrid.h func.h SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
sub.o: sub.c hybrid.h func.h SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
timestep.o: timestep.c hybrid.h func.h SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
SFMT.o: SFMT.c SFMTdir/SFMT.h SFMTdir/SFMT-params.h SFMTdir/SFMT-params216091.h
endif

.PHONY: all clean
