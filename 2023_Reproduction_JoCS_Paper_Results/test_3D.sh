# Create and enter the test folder
mkdir Test_3D
cd Test_3D

# Mesh dimension
nDim=3

# Time when we want to compute the flowmap
tEval=8

# Number of RK4 steps
nRK4=1

# Print result to file (disabled to complete the tests faster)
toFile=0

# Global mesh dimension (thousands) -  In this case, 200K
N=200

# Folder containing the mesh data files
source_folder="../../source/3D_200K"

# Mesh axis steps
nx=58
ny=58
nz=58

# Number of known time instants (and last one)
nt=500
tlim=10

# Execute the tests

# A) Sequential tests
th=1
export OMP_NUM_THREADS=$th
sched=1
# Combination of compiler and optimization flag to test
for comp in "gcc_O3" "gcc_Of" "icc_O3" "icc_Of" "aocc_O3" "aocc_Of"
do
	../bin/${comp}_compute_flowmap $nDim $tEval ${source_folder}/faces.txt ${source_folder}/times.txt $nRK4 $sched $toFile $nx $ny $nz $nt $tlim > ${nDim}D_${N}K_${comp}_${th}th_sched${sched}.txt
done

# B) Parallel tets
# Number of OpenMP threads to use
for th in 12 24 36 48 72 96
do
    export OMP_NUM_THREADS=$th
    # Sequential (1) or OpenMP parallel versions (static - 2, dynamic - 3, guided - 4)
    for sched in 2 3 4
    do
        # Combination of compiler and optimization flag to test
        for comp in "gcc_O3" "gcc_Of" "icc_O3" "icc_Of" "aocc_O3" "aocc_Of"
        do
            ../bin/${comp}_compute_flowmap $nDim $tEval ${source_folder}/faces.txt ${source_folder}/times.txt $nRK4 $sched $toFile $nx $ny $nz $nt $tlim > ${nDim}D_${N}K_${comp}_${th}th_sched${sched}.txt 
        done
    done
done

# Exit the test folder
cd ..
