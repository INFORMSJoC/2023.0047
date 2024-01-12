#!/bin/sh

# RUN NAME
runname="run_name"

# PATHS
HOMEPATH="/home/user"
IRELPATH="realtive/path/to/instance_files"
IPATH=$HOMEPATH/$IRELPATH
RESULTS=$HOMEPATH/results/$runname/
LOGS=$HOMEPATH/logs/$runname
mkdir -p $RESULTS
mkdir -p $LOGS


# SET STANDARD PARAMETERS
ITYPE=1
DEBUG=1
RANDOMSEED=1
NSCENARIOS=10
NSCENARIOSEVAL=100000
KF=5
KL=0
CONSOLE="false"
PRECOMPILE="true"
WRITEOUTPUT="true"
HCALLBACK="true"
BENF="true"
PREPRO="false"
PROBF=0.1
PROBV=0.5
SSLTYPE=3
CLEVEL=0.1
HTYPE=3
OATYPE=0
LPMETHOD=0
RELCUTTOL=1.005
RELCUTTOLU=1.001
SOLVEROOT="false"
ROOTCUTS="false"
ROOTCUTSALPHA=0.5
ROOTCUTSLAMBDA=0.5
FRACTIONOFRESISTANTNODES=0.5
COMPUTESL="false"
PTYPE=1
NSSAITR=10
TIMELIMIT=7200
LEADERUTILITY=1
MEMLIMIT=15500
HVMEM=15.8G


# START RUNS
for files in $IPATH/*
do		
	IFILE=$(basename "$files")
	if [ -f $IPATH/$IFILE ]; then  
		echo "$IFILE is a file (not a directory)" 
	else
		continue #for
        fi

	[ "$IFILE" = "D-n25-k4-b0.3-beta2.0-18.0-i1" ] && continue  

for KF in 5 10 15 #5 10 15
do
for KL in 0 # #0 5 10 15 
do
for NSCENARIOS in 100 #100 250 500 750
do
for FRACTIONOFRESISTANTNODES in 0 #0.25 0.5 0.75
do
for MTYPE in 5 #1 2 3 4 5
do
for LEADERUTILITY in 0 #0 1 8 20 55 150
do
for PREPRO in "true" #"true" "false"
do
for HTYPE in 0	#0 3
do
for SSLTYPE in 4 #1 2 3 4 5
do
for RELCUTTOL in 1.005 #1.005 1.01 1.02 
do
for RELCUTTOLU in 1.01 #1.005 1.01 1.02 
do
for PROBV in 0.05 
do
for CUTTYPE in 5 #1 2 3 4 5 6 7 8 9
do
for OATYPE in 0 
do
for LPMETHOD in 0
do
for ROOTCUTSALPHA in 0.3
do
for ROOTCUTSLAMBDA in 0.5
do
for ROOTCUTS in "false"
do


# SET PARAMETERS
PARAMS="--ipath $IPATH --ifile $IFILE --itype $ITYPE --opath $RESULTS --nscenarios $NSCENARIOS --Nscenarios $NSCENARIOSEVAL --debug $DEBUG --console $CONSOLE --precompile $PRECOMPILE --writeoutput $WRITEOUTPUT --kF $KF --kL $KL --probF $PROBF --mtype $MTYPE --timelimit $TIMELIMIT --memlimit $MEMLIMIT --hcallback $HCALLBACK --benf $BENF --rseed $RANDOMSEED --ssltype $SSLTYPE --clevel $CLEVEL --prepro $PREPRO --htype $HTYPE --solveRoot $SOLVEROOT --ptype $PTYPE --nSSAItr $NSSAITR --computeSL $COMPUTESL --relcuttol $RELCUTTOL --relcuttolU $RELCUTTOLU --probV $PROBV --cutType $CUTTYPE --OAtype $OATYPE --LPMethod $LPMETHOD --rootCuts $ROOTCUTS --rootCutsAlpha $ROOTCUTSALPHA --rootCutsLambda $ROOTCUTSLAMBDA --fractionOfResistantNodes $FRACTIONOFRESISTANTNODES --leaderUtility $LEADERUTILITY"


# SUBMIT
OUTFILE=$LOGS/$IFILE-nsc$NSCENARIOS-kF$KF-kL$KL-ht$HTYPE-sr$SOLVEROOT-rc$ROOTCUTS-pV$PROBV-fr$FRACTIONOFRESISTANTNODES-rt$RELCUTTOL-rtU$RELCUTTOLU-lu$LEADERUTILITY-pp$PREPRO-m$MTYPE.log
ERRFILE=$LOGS/$IFILE-nsc$NSCENARIOS-kF$KF-kL$KL-ht$HTYPE-sr$SOLVEROOT-rc$ROOTCUTS-pV$PROBV-fr$FRACTIONOFRESISTANTNODESrt$RELCUTTOL-rtU$RELCUTTOLU-lu$LEADERUTILITY-pp$PREPRO-m$MTYPE.err
JOBNAME=$IFILE-nsc$NSCENARIOS-kF$KF-kL$KL-ht$HTYPE-sr$SOLVEROOT-rc$ROOTCUTS--pV$PROBV-fr$FRACTIONOFRESISTANTNODES-rt$RELCUTTOL-rtU$RELCUTTOLU-lu$LEADERUTILITY-pp$PREPRO-m$MTYPE		
echo $JOBNAME 
# submit to your favorite queue, e.g.:
#qsub -N $JOBNAME $emailflag -r y -j y -l h_vmem=$HVMEM -l s_core=0 -l longrun=0 -o $OUTFILE -e $ERRFILE ./runprogram_logit.sh "$PARAMS"

done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
