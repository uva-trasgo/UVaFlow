<p align="center">
  <img src="https://github.com/uva-trasgo/UVaFlow/blob/master/UVaFlow_Logo.png"
 </p>

A parallel software (based on OpenMP) to compute flowmaps in the extraction of Lagrangian Coherent Structures (LCS) in 2D/3D scenarios.

Currently supported meshes must be composed of simplex faces, this is, triangles (in 2D) or tetrahedrons (in 3D).

## Compilation from scratch

All the make-derived files will be stored in the ```bin``` folder. 
__Make sure you have created this folder (```mkdir bin```) before compiling any version of the code using the Makefile provided.__

The __Makefile__ is currently ready to compile the code using the following compilers:
- GCC Compiler (type ```make gcc``` to compile both):
  - equipped with -O3 optimization flag: ```make gcc_O3_compute_flowmap```
  - equipped with -Ofast optimization flag: ```make gcc_Of_compute_flowmap```
- ICC Compiler (type ```make icc``` to compile both):
  - equipped with -O3 optimization flag: ```make icc_O3_compute_flowmap```
  - equipped with -Ofast optimization flag: ```make icc_Of_compute_flowmap```
- CLANG Compiler (type ```make aocc``` to compile both):
  - equipped with -O3 optimization flag: ```make aocc_O3_compute_flowmap```
  - equipped with -Ofast optimization flag: ```make aocc_Of_compute_flowmap```
  
NOTE: By default, ```make``` will compile the GCC and AOCC versions.
 
Besides, to see the __vectorization reports__ for each compiler/flag combination, type:
- GCC Compiler
  - Equipped with -O3 optimization flag: ```make gcc_O3_compute_flowmap_vect```
  - Equipped with -Ofast optimization flag: ```make gcc_Of_compute_flowmap_vect```
- ICC Compiler
  - Equipped with -O3 optimization flag: ```make icc_O3_compute_flowmap_vect```
  - Equipped with -Ofast optimization flag: ```make icc_Of_compute_flowmap_vect```
- CLANG Compiler
  - Equipped with -O3 optimization flag: ```make aocc_O3_compute_flowmap_vect```
  - Equipped with -Ofast optimization flag: ```make aocc_Of_compute_flowmap_vect```

Type ```make clean``` to empty the ```bin``` folder.

## Data files

The mesh data must be provided separated in files:
- ```coords.txt```: contains the mesh points coordinates (one value per line)
- ```faces.txt```: contains the mesh points indices (associated to the points order in coords.txt) that form each mesh face (one per line)
- ```times.txt```: contains the time instants for which the velocity is knwon (one per line)
- ```velocity.txt```: contains the knwon velocity values (one value per line) associated to each time in times.txt and each point in coords.txt

## Data generation

The ```UVaFlow_mesh-generation.py``` Python script automatically generates the four files described in the previous section, either for the 2D Double Gyre and 3D ABC flows. Please, note that the mesh files generation of big meshes might take some time to complete.

Remarks:

- In this script, the 2D Double Gyre Flow region is assumed to be [0:2, 0:1], this is, X axis values are taken in the [0, 2] interval, and Y axis values are taken in the [0, 1] interval. If one wants to change this, substitute ```0, 2``` in line 78 (```x_ = np.linspace(0, 2, int(args.x_steps_axis))```) and/or ```0, 1``` in line 79 (```y_ = np.linspace(0, 1, int(args.y_steps_axis))```).

- In this script, the 3D ABC Flow is assumed to be [0:1, 0:1, 0:1], this is, all the axis values are taken in the [0, 1] interval. If one wants to change this, substitute ```0, 1``` in line 158 (```x_ = np.linspace(0, 1, int(args.x_steps_axis))```), the ```0, 1``` in line 159 (```y_ = np.linspace(0, 1, int(args.y_steps_axis))```), and/or the ```0, 1``` in line 160 (```z_ = np.linspace(0, 1, int(args.x_steps_axis))```).

Usage (type ```python3 UVaFlow_mesh-generation.py -h``` to see the detailed help regarding usage): 

```UVaFlow_mesh-generation.py [-h] <nDim> <nt_knwon> <t_zero> <t_end> <coords_file> <faces_file> <times_file> <vel_file> <num_cores> <x_steps_axis> <y_steps_axis> [<z_steps_axis>]```

Where:
- ```nDim```: "2D_DGyre" or "3D_ABC"
- ```nt_knwown```: Number of time instants for which the velocity is knwon
- ```t_zero```: First time instant for which the velocity is knwon
- ```t_end```: Last time instant for which the velocity is knwon
- ```coords_file/faces_file/times_file/vel_file```: File path where the coords/faces/times/velocities information will be stored
- ```num_cores```: CPU threads to use in the mesh generation
- ```x_steps_axis/y_steps_axis/z_steps_axis```: Steps in X/Y/Z axis

Samples:

- 2D sample: Generation of a 2D Double Gyre Flow mesh knowing ${\color{red}100}$ time instants (between ${\color{green}t0=0,\ tend=10}$), composed of ${\color{blue}200*100}$ mesh points using ${\color{orange}1}$ thread and storing the generated data in the "coords.txt", "faces.txt", "times.txt", "velocity.txt" files.

  > __python3 UVaFlow_mesh-generation.py 2D_DGyre ${\color{red}100}$ ${\color{green}0\ 10}$ coords.txt faces.txt times.txt velocity.txt ${\color{orange}1}$ ${\color{blue}200\ 100\ 0}$__

- 3D sample: Generation of a 3D ABC Flow mesh knowing ${\color{red}50}$ time instants (between ${\color{green}t0=0,\ tend=1}$), composed of ${\color{blue}100\*100\*100}$ mesh points using ${\color{orange}8}$ threads and storing the generated data in the "coords.txt", "faces.txt", "times.txt", "velocity.txt" files.

  > __python3 UVaFlow_mesh-generation.py 3D_ABC ${\color{red}50}$ ${\color{green}0\ 1}$ coords.txt faces.txt times.txt velocity.txt ${\color{orange}8}$ ${\color{blue}100\ 100\ 100}$__
  
## Flowmap computation

Usage:

```bin/<executable> <nDim> <t_eval> <coords_file> <faces_file> <times_file> <velocity_file> <nsteps_rk4> <sched_policy> <print>```

Where:
- ```executable```: gcc_O3_compute_flowmap, gcc_Of_compute_flowmap, aocc_O3_compute_flowmap, aocc_Of_compute_flowmap... the executable stored in the ```bin``` folder that corresponds to the compiler (gcc, clang, icc) and optimization flag (-O3, -Ofast) desired (see the "Compilation from scratch" section for more details).
- ```nDim```: either ```2``` (for 2D Flows) or ```3``` (for 3D flows)
- ```t_eval```: time instant when the user wants to compute the flowmap
- ```coords_file/faces_file/times_file/velocity_file```:  file where mesh coordinates/faces/times/velocities are stored (see the previous data files section for more information).
- ```nsteps_rk4```: number of iterations to perform in the RK4 call.
- ```sched_policy```: either perform a sequential execution (```1```) or a parallel one based on OpenMP static (```2```), dynamic (```3```) or guided (```4```) scheduling policy
- ```print```: indicate whether the final result must be stored in an output file or not (```0```-NO, ```1```-YES)

Remainder: Remember to properly set the number of CPU threads to use with any of the OpenMP based versions, for example with ```export OMP_NUM_THREADS=X``` being X the number desired.

Samples:

- 2D sample: Computation of the flowmap using the GCC compiler provided with the -O3 optimization flag, for a ${\color{blue}2D}$ mesh at ${\color{red}t=8}$ having the data files stored in the ```source/dgyre_input/``` folder, using the ${\color{green}OpenMP dynamic}$ based parallel version, ${\color{orange}1\ step}$ for the RK4 calls, and ${\color{purple}not\ printing}$ the output in any file.

  > __bin/gcc_O3_compute_flowmap ${\color{blue}2}$ ${\color{red}8}$ source/dgyre_input/coords.txt source/dgyre_input/faces.txt source/dgyre_input/times.txt source/dgyre_input/velocity.txt ${\color{orange}1}$ ${\color{green}3}$ ${\color{purple}0}$__
  
- 3D sample: Computation of the flowmap using the AOCC compiler provided with the -Ofast optimization flag, for a ${\color{blue}3D}$ mesh at ${\color{red}t=0.5}$ having the data files stored in the ```source/abc_input/``` folder, using the ${\color{green}OpenMP static}$ based parallel version, ${\color{orange}1\ step}$ for the RK4 calls, and ${\color{purple}not\ printing}$ the output in any file.

  > __bin/aocc_Of_compute_flowmap ${\color{blue}3}$ ${\color{red}0.5}$ source/abc_input/coords.txt source/abc_input/faces.txt source/abc_input/times.txt source/abc_input/velocity.txt ${\color{orange}1}$ ${\color{green}2}$ ${\color{purple}0}$__

## How to cite

<a id="1">[1]</a> 
Carratalá-Sáez, R. and Sierra-Pallares, J. and Llanos, D. R. and Gonzalez-Escribano, A. (2023). 
UVaFlow: Lagrangian flowmap computation for fluid dynamic applications.
Submitted to the Journal Of Computational Science [major review revision in process], -(-), pp. -.

### Reproduction of the results in that paper

We provide two different paths to reproduce the results shown in the paper.

- Option A implies preparing everything from scratch, as a regular UVaFlow software user.
- Option B is a suggestion we have prepared in order to make the reproduction faster. 

The only difference between them is the fact that Option A uses ```src/compute_flowmap.c```, while B uses ```2023_Reproduction_JoCS_Paper_Results/compute_flowmap.c``` when compiling the code. This second C file does not read the coordinates file, neither the velocity file, but instead automatically generate that data. This is why this option is faster than the previous one.

#### Option A (slower)

1. Make sure that you have created the following folders in your working directory:

    ```
    mkdir bin
    mkdir source
    mkdir source/2D_500K
    mkdir source/3D_200K
    ```

2. Generate the data files associated to the test cases illustrated in the paper (modify the 10 after "velocity.txt" to specify the number of threads to use based on your system).

    2D Double Gyre Flow (500K): 
    ```
    python3 UVaFlow_mesh-generation.py 2D_DGyre 100 0 10 source/2D_500K/coords.txt source/2D_500K/faces.txt source/2D_500K/times.txt source/2D_500K/velocity.txt 10 1000 500 0
    ```
    3D ABC Flow (200K):
    ```
    python3 UVaFlow_mesh-generation.py 3D_ABC 500 0 10 source/3D_200K/coords.txt source/3D_200K/faces.txt source/3D_200K/times.txt source/3D_200K/velocity.txt 10 58 58 58
    ```

3. Compile the desired version of the code, following the instructions detailed in the "Compilation from scratch" section.
    ```
    # For example, both GCC versions (-O3 and -Ofast)
    make gcc
    ```
   
4. Test. Either use ```test_2D.sh``` and ```test_3D.sh``` or modified versions of them.

#### Option B (faster)

1. Enter the ```2023_Reproduction_JoCS_Paper_Results``` folder.
    ```
    cd 2023_Reproduction_JoCS_Paper_Results
    ```

2. Make sure that you have the ```bin``` folder created.
    ```
    mkdir bin
    ```

3. Compile the code available in that folder.
    ```
    # Compile the desired combination of compiler and optimization flag
    # Follow the instructions detailed in the "Compilation from scratch" section
    
    # For example, both GCC versions (-O3 and -Ofast)
    make gcc
    ```
    
 4. Test. Either use ```test_2D.sh``` and ```test_3D.sh``` or modified versions of them.
