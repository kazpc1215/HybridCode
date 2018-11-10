#!/bin/sh

history -r .input_history

while [ -z $SCRIPT ]
do
    read -e -r -p "Type the shell script name you want to run. : " SCRIPT
done
echo "$SCRIPT"
history -s "$SCRIPT"

while [ -z $YESNO ]
do
    read -e -r -p "Is there a dependent job? (y / n) : " YESNO
done
echo "$YESNO"
history -s "$YESNO"

if [ $YESNO = 'y' ] ; then
    
    while [ -z $JOBID ]
    do
        read -e -r -p "Type a number to set the dependent job ID. (ex. *****.sdb) : " JOBID
    done
    echo "$JOBID"
    history -s "$JOBID"
    
    while [ -z $NRUN ]
    do
	read -e -r -p "Type a number to decide how many times to execute the job.: " NRUN
    done
    echo "$NRUN"
    history -s "$NRUN"

    for i in `seq 1 $NRUN` ; do
        qsub -W depend=afterok:$JOBID $SCRIPT | tee tmp.txt
        JOBID=$(cat tmp.txt)
        echo ""
        echo "-----dependency of $JOBID-----"
        qstat -f $JOBID | grep depend
        echo -n "------------------------"
	for n in `seq 1 ${#JOBID}` ; do
            echo -n "-"
        done
        echo ""
	echo ""
    done
    
elif [ $YESNO = 'n' ] ; then

    while [ -z $NRUN ]
    do
        read -e -r -p "Type a number to decide how many times to execute the job.: " NRUN
    done
    echo "$NRUN"
    history -s "$NRUN"
    
    for i in `seq 1 $NRUN` ; do
	if [ $i = 1 ] ; then
	    qsub $SCRIPT | tee tmp.txt
	else
	    qsub -W depend=afterok:$JOBID $SCRIPT | tee tmp.txt
	fi
        JOBID=$(cat tmp.txt)
	echo ""
        echo "-----dependency of $JOBID-----"
        qstat -f $JOBID | grep depend
        echo -n "------------------------"
	for n in `seq 1 ${#JOBID}` ; do
	    echo -n "-"
	done
	echo ""
	echo ""
    done

else
	    
    echo "Input 'y' or 'n'."
    
fi

history -w .input_history
rm -f tmp.txt
