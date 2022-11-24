# Compilers
CC=gcc
CLANG=clang

# Flags
FLAGS=-lm
FLAG_OMP=-fopenmp

# Directories
DIR=.
DIR_src=${DIR}/src
DIR_bin=${DIR}/bin

# Complementary files
SRC=${DIR_src}/kdtree.c ${DIR_src}/interpolation.c ${DIR_src}/location.c ${DIR_src}/rk4.c ${DIR_src}/distance.c ${DIR_src}/read_python_files.c

# Include header files
INC=./include

# Compile lists

GCC_OBJS=gcc_O3_compute_flowmap gcc_Of_compute_flowmap
AOCC_OBJS=aocc_O3_compute_flowmap aocc_Of_compute_flowmap

OBJS=${GCC_OBJS} ${AOCC_OBJS}

# Make lists
all:     ${OBJS}
gcc:     ${GCC_OBJS}
aocc:    ${AOCC_OBJS}

# -------------------------- #
# ---------- GCC ----------- #
# -------------------------- #

gcc_O3_compute_flowmap:
	${CC} ${FLAG_OMP} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/gcc_O3_compute_flowmap ${FLAGS} -O3

gcc_Of_compute_flowmap:
	${CC} ${FLAG_OMP} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/gcc_Of_compute_flowmap ${FLAGS} -Ofast

# -------------------------- #
# ---------- AOCC----------- #
# -------------------------- #

aocc_O3_compute_flowmap:
	${CLANG} ${FLAG_OMP} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/aocc_O3_compute_flowmap ${FLAGS} -O3

aocc_Of_compute_flowmap:
	${CLANG} ${FLAG_OMP} ${DIR_src}/compute_flowmap.c ${SRC} -I ${INC} -o ${DIR_bin}/aocc_Of_compute_flowmap ${FLAGS} -Ofast

clean:
	cd ${DIR_bin} && rm ${OBJS} && cd ..
