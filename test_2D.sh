ompth=8
ndim=2
ntimes_eval=1
t0_eval=8
tdelta=0
nrk4=1
sched=2
chunk_size=1
print2file=1
export OMP_NUM_THREADS=$ompth
#./bin/gcc_O3_compute_flowmap $ndim $ntimes_eval $t0_eval $tdelta source/doublegire_input/coords.txt source/doublegire_input/faces.txt source/doublegire_input/times.txt source/doublegire_input/velocity.txt $nrk4 $sched $chunk_size $print2file
./bin/gcc_O3_compute_flowmap $ndim $t0_eval source/doublegire_input/coords.txt source/doublegire_input/faces.txt source/doublegire_input/times.txt source/doublegire_input/velocity.txt $nrk4 $sched $chunk_size $print2file
