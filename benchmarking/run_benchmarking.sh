#!/usr/bin/bash
# nohup ./run_benchmarking > ../benchmarking/log.out &; echo $! > pid
if [ -z $2 ]; then
	#methods=("no_decontamination" "soupx:autoEstCont" "soupx:background_genes" "soupx:top_background_genes" "decontx:no_cell_types" "decontx:with_cell_types")
	methods=("fastcar")
	if [ -z $1 ]; then
		sample_ids=("BG1_BG20C" "BG3_BG21C" "BG5_BG22C" "BG52_BG20C_MeOH" "BG54_BG21C_MeOH" "BG56_BG22C_MeOH")
	else
		sample_ids=$1
	fi
else
	methods=$2
fi

USE_GPU='true'

dir=/data/Perkins/benchmarking/output_gpu

# START PIDSTAT for RAM AND CPU
nohup pidstat -h -r -u -C R 1 --human >> $dir/r_usage.txt &
pid_r=$!
nohup pidstat -h -r -u -C python3 1 --human >> $dir/python.txt &
pid_c=$!

echo $pid_r >> pid
echo $pid_c >> pid

echo "PID METHOD SAMPLE_ID" > $dir/pid.txt
echo "START" > $dir/output.log

# looping through all methods and samples
for i in {0..0}; do
	m=${methods[$i]}
	echo $m

	for id in ${sample_ids[@]}; do
		# record GPU usage > for cellbender run w/ cuda		
		if [ $USE_GPU = 'true' ]; then
			echo $m $id >> $dir/gpu_usage.csv
			nohup nvidia-smi --query-gpu=utilization.gpu,utilization.memory,memory.total,memory.free,memory.used --format=csv -l 1 >> $dir/gpu_usage.csv & 
			nvidia_pid=$!
		fi	
		
		# run r script
		nohup Rscript benchmarking/benchmarking.R $id $m $i >> $dir/output.log &

		pid=$!
		# record pid
		echo $pid $m $id >> $dir/pid.txt
		
		# wait for it to finish
		tail --pid=$pid -f /dev/null
		
		if [ $USE_GPU = 'true' ]; then
			kill -9 $nvidia_pid
		fi
		# wait 5 seconds
		sleep 5
	done
done

# killing pidstat processes
kill -9 $pid_r
kill -9 $pid_c
