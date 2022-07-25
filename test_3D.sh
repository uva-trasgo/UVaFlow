ompth=8
ndim=3
ntimes_eval=1
t0_eval=8
tdelta=0
nrk4=1
sched=2
chunk_size=1
print2file=0
export OMP_NUM_THREADS=$ompth
./bin/gcc_O3_compute_flowmap $ndim $ntimes_eval $t0_eval $tdelta source/abc_input/coords.txt source/abc_input/faces.txt source/abc_input/times.txt source/abc_input/velocity.txt $nrk4 $sched $chunk_size $print2file
