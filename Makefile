# Compilers
CC=gcc
CLANG=clang
ICC=icc

# Flags
FLAGS=-lm
FLAG_OMP=-fopenmp
FLAG_ZEN3=-march=znver3

# Directories
DIR=.
DIR_src=${DIR}/src
DIR_bin=${DIR}/bin

# Complementary files
SRC=${DIR_src}/kdtree.c ${DIR_src}/search.c ${DIR_src}/interpolation.c ${DIR_src}/location.c ${DIR_src}/rk4.c ${DIR_src}/distance.c ${DIR_src}/read_python_files.c

# Include header files
INC=./include

# Compile lists

GCC_OBJS=gcc_O3_compute_flowmap gcc_Of_compute_flowmap
ICC_OBJS=icc_O3_compute_flowmap icc_Of_compute_flowmap

GCC_OBJS_Z3=gcc_O3_compute_flowmap_z3 gcc_Of_compute_flowmap_z3
AOCC_OBJS=aocc_O3_compute_flowmap aocc_Of_compute_flowmap

OBJS=${GCC_OBJS} ${GCC_OBJS_Z3} ${AOCC_OBJS} ${ICC_OBJS}

# Make lists
all:     ${OBJS}
gcc:     ${GCC_OBJS}
icc:     ${ICC_OBJS}
gccz3:   ${GCC_OBJS_Z3}
aocc:    ${AOCC_OBJS}

# -------------------------- #
# ---------- GCC ----------- #
# -------------------------- #

gcc_O3_compute_flowmap:
	${CC} ${FLAG_OMP} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/gcc_O3_compute_flowmap ${FLAGS} -O3

gcc_Of_compute_flowmap:
	${CC} ${FLAG_OMP} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/gcc_Of_compute_flowmap ${FLAGS} -Ofast

gcc_O3_compute_flowmap_z3:
	${CC} ${FLAG_OMP} ${FLAG_ZEN3} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/gcc_O3_compute_flowmap_z3 ${FLAGS} -O3

gcc_Of_compute_flowmap_z3:
	${CC} ${FLAG_OMP} ${FLAG_ZEN3} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/gcc_Of_compute_flowmap_z3 ${FLAGS} -Ofast

# -------------------------- #
# ---------- ICC ----------- #
# -------------------------- #

icc_O3_compute_flowmap:
	${ICC} ${FLAG_OMP} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/icc_O3_compute_flowmap ${FLAGS} -O3

icc_Of_compute_flowmap:
	${ICC} ${FLAG_OMP} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/icc_Of_compute_flowmap ${FLAGS} -Ofast

# -------------------------- #
# ---------- AOCC----------- #
# -------------------------- #

aocc_O3_compute_flowmap:
	${CLANG} ${FLAG_OMP} ${FLAG_ZEN3} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/aocc_O3_compute_flowmap_z3 ${FLAGS} -O3

aocc_Of_compute_flowmap:
	${CLANG} ${FLAG_OMP} ${FLAG_ZEN3} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/aocc_Of_compute_flowmap_z3 ${FLAGS} -Ofast

# -------------------------- #
# ---------- Clean---------- #
# -------------------------- #

clean:
	cd ${DIR_bin} && rm ${OBJS} && cd ..
